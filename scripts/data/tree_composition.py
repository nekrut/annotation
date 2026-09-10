#!/usr/bin/env python3
"""Summarise the composition of a multiple-alignment tree: how many leaves
each clade contributes and how far they sit from the reference.

Standard library only.  Input is a Newick file with branch lengths (the
``.nh`` that UCSC serves next to each multiz or Cactus track) and,
optionally, a TSV that maps leaf names to a clade label.  Output is one
row per clade with the number of leaves, the number of distinct species
names (UCSC alignments can carry two assemblies of one species), and the
minimum, median and maximum patristic distance from the reference leaf in
substitutions per site, plus how many leaves fall within and beyond the
non-coding alignability horizon measured in ``docs/data-sources.md``
section 6 (mostly aligned below 0.5, essentially unaligned past 1.0).

The clade TSV has a header and at least the columns ``leaf`` and
``clade``; extra columns are ignored.  ``--lookup-gbif`` builds that TSV
from the GBIF species-match API (order and family per leaf; one request
per leaf, so about ten minutes for 470 leaves) and writes it to the path
given, after which the run is offline.  Leaves absent from the TSV are
reported under ``(unassigned)``.

Usage::

    python3 scripts/data/tree_composition.py hg38.470way.scientificNames.nh \
        --reference Homo_sapiens --clades scripts/data/hg38.470way.orders.tsv --markdown
    python3 scripts/data/tree_composition.py hg38.470way.scientificNames.nh \
        --reference Homo_sapiens --lookup-gbif scripts/data/hg38.470way.orders.tsv

With no ``--clades`` the script prints one row for the whole tree and a
histogram of distances at the bin edges given by ``--bins``.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict


# GBIF splits and lumps a few mammalian orders differently from the NCBI
# taxonomy that the rest of the project uses.  Fold them so that whales sit
# with the even-toed ungulates and shrews with hedgehogs.
CLADE_ALIASES = {
    "Cetacea": "Artiodactyla",
    "Soricomorpha": "Eulipotyphla",
    "Erinaceomorpha": "Eulipotyphla",
}


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None


def parse_newick(text: str) -> Node:
    """Minimal Newick parser: names, branch lengths, nested parentheses."""
    s = re.sub(r"\s+", "", text)
    pos = 0

    def read_label():
        nonlocal pos
        m = re.match(r"[^:,();]*", s[pos:])
        pos += len(m.group(0))
        return m.group(0)

    def read_length():
        nonlocal pos
        if pos < len(s) and s[pos] == ":":
            m = re.match(r":([-+0-9.eE]+)", s[pos:])
            pos += len(m.group(0))
            return float(m.group(1))
        return 0.0

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                child = node()
                child.parent = n
                n.children.append(child)
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError("bad Newick at offset %d" % pos)
        n.name = read_label()
        n.length = read_length()
        return n

    root = node()
    return root


def leaves(root: Node):
    out = []
    stack = [root]
    while stack:
        n = stack.pop()
        if n.children:
            stack.extend(reversed(n.children))
        else:
            out.append(n)
    return out


def path_to_root(n: Node):
    p = []
    while n is not None:
        p.append(n)
        n = n.parent
    return p


def patristic(a: Node, b: Node) -> float:
    pa = path_to_root(a)
    seen = {id(n): i for i, n in enumerate(pa)}
    d = 0.0
    n = b
    while id(n) not in seen:
        d += n.length
        n = n.parent
    for m in pa[: seen[id(n)]]:
        d += m.length
    return d


def species_name(leaf: str) -> str:
    """Collapse subspecies and assembly suffixes to a binomial for counting."""
    parts = leaf.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else leaf


def gbif_lookup(name: str, retries: int = 3):
    q = name.replace("_", " ")
    for attempt in range(retries):
        try:
            url = "https://api.gbif.org/v1/species/match?name=" + urllib.parse.quote(q)
            with urllib.request.urlopen(url, timeout=30) as r:
                j = json.load(r)
            if "order" not in j and len(q.split()) > 2:
                url = "https://api.gbif.org/v1/species/match?name=" + urllib.parse.quote(
                    " ".join(q.split()[:2])
                )
                with urllib.request.urlopen(url, timeout=30) as r:
                    j = json.load(r)
            return j
        except Exception as exc:  # network hiccup: back off and retry
            err = str(exc)
            time.sleep(2 * (attempt + 1))
    return {"error": err}


def write_clades_from_gbif(leaf_nodes, path):
    rows = []
    for n in leaf_nodes:
        j = gbif_lookup(n.name)
        rows.append(
            {
                "leaf": n.name,
                "clade": j.get("order", ""),
                "family": j.get("family", ""),
                "class": j.get("class", ""),
                "match": j.get("matchType", j.get("error", "")),
                "source": "gbif",
            }
        )
        time.sleep(0.05)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["leaf", "clade", "family", "class", "match", "source"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    missing = [r["leaf"] for r in rows if not r["clade"]]
    print("wrote %s (%d leaves, %d without a clade: %s)" % (path, len(rows), len(missing), ", ".join(missing)), file=sys.stderr)


def read_clades(path):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            c = row.get("clade", "").strip()
            out[row["leaf"]] = CLADE_ALIASES.get(c, c) or "(unassigned)"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("newick")
    ap.add_argument("--reference", required=True, help="leaf name of the reference species")
    ap.add_argument("--clades", help="TSV with columns leaf, clade")
    ap.add_argument("--lookup-gbif", metavar="TSV", help="query GBIF for every leaf and write this TSV, then use it")
    ap.add_argument("--bins", default="0.1,0.25,0.5,1,2", help="distance bin edges for the histogram")
    ap.add_argument("--near", type=float, default=0.5, help="'within' threshold, substitutions per site")
    ap.add_argument("--far", type=float, default=1.0, help="'beyond' threshold, substitutions per site")
    ap.add_argument("--markdown", action="store_true", help="print a Markdown table instead of TSV")
    ap.add_argument("--tsv", help="also write the per-leaf distance table here")
    args = ap.parse_args(argv)

    root = parse_newick(open(args.newick).read())
    lv = leaves(root)
    by_name = {n.name: n for n in lv}
    if args.reference not in by_name:
        sys.exit("reference %r is not a leaf; leaves look like %s" % (args.reference, [n.name for n in lv[:3]]))
    ref = by_name[args.reference]

    if args.lookup_gbif:
        write_clades_from_gbif(lv, args.lookup_gbif)
        args.clades = args.lookup_gbif
    clade_of = read_clades(args.clades) if args.clades else {}

    dist = {n.name: patristic(ref, n) for n in lv}
    informants = [n for n in lv if n is not ref]

    if args.tsv:
        with open(args.tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["leaf", "species", "clade", "distance_from_reference"])
            for n in sorted(lv, key=lambda n: dist[n.name]):
                w.writerow([n.name, species_name(n.name), clade_of.get(n.name, "(unassigned)"), "%.4f" % dist[n.name]])

    groups = defaultdict(list)
    for n in informants:
        groups[clade_of.get(n.name, "(unassigned)") if clade_of else "all"].append(n)

    header = ["clade", "leaves", "species", "min", "median", "max", "within %.2g" % args.near, "beyond %.2g" % args.far]
    rows = []
    for clade, ns in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        ds = [dist[n.name] for n in ns]
        rows.append(
            [
                clade,
                len(ns),
                len({species_name(n.name) for n in ns}),
                "%.3f" % min(ds),
                "%.3f" % statistics.median(ds),
                "%.3f" % max(ds),
                sum(d <= args.near for d in ds),
                sum(d > args.far for d in ds),
            ]
        )
    ds = [dist[n.name] for n in informants]
    rows.append(
        [
            "total",
            len(informants),
            len({species_name(n.name) for n in informants}),
            "%.3f" % min(ds),
            "%.3f" % statistics.median(ds),
            "%.3f" % max(ds),
            sum(d <= args.near for d in ds),
            sum(d > args.far for d in ds),
        ]
    )

    edges = [float(x) for x in args.bins.split(",")]
    hist = []
    lo = 0.0
    for hi in edges + [float("inf")]:
        hist.append(("%g-%g" % (lo, hi) if hi != float("inf") else ">%g" % lo, sum(lo < d <= hi for d in ds)))
        lo = hi

    if args.markdown:
        print("| " + " | ".join(header) + " |")
        print("|" + "---|" * len(header))
        for r in rows:
            print("| " + " | ".join(str(x) for x in r) + " |")
        print()
        print("| distance (subst/site) | informants |")
        print("|---|---|")
        for label, c in hist:
            print("| %s | %d |" % (label, c))
    else:
        print("\t".join(header))
        for r in rows:
            print("\t".join(str(x) for x in r))
        print()
        for label, c in hist:
            print("%s\t%d" % (label, c))
    dup = [s for s, c in _counts(species_name(n.name) for n in lv).items() if c > 1]
    print(
        "%d leaves, %d distinct species; %d species with more than one assembly%s"
        % (len(lv), len({species_name(n.name) for n in lv}), len(dup), (": " + ", ".join(sorted(dup))) if dup else ""),
        file=sys.stderr,
    )


def _counts(it):
    out = defaultdict(int)
    for x in it:
        out[x] += 1
    return out


if __name__ == "__main__":
    main()
