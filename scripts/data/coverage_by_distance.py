#!/usr/bin/env python3
"""Aggregate the per-species alignment coverage tables written by
``fetch_window.py`` over many windows, binned by tree distance from the
reference.

Standard library only.  Input is one or more directories (or manifest
files) produced by ``fetch_window.py``; each ``<name>.manifest.json`` carries
``coverage.species[<informant>][cds|utr|intron|intergenic]`` (fraction of the
window's reference bases of that class at which the informant has an
aligned, non-gap base) and ``coverage.window_bases_by_class``; the matching
``<name>.nh`` is the track's Newick tree with branch lengths in
substitutions per site.  The distance used is the patristic distance from
the reference leaf to the informant leaf on that tree.

For each distance bin and annotation class the script reports the
base-weighted coverage (sum over windows and informants of covered bases,
divided by the sum of class bases), so long windows count more than short
ones and every informant in the bin counts once per window.  It also
reports the share of informant-window pairs with coverage above 0.5, which
is the number a model designer cares about: how often is there anything to
look at.

Usage::

    python3 scripts/data/coverage_by_distance.py /tmp/cov/fly --reference dm6
    python3 scripts/data/coverage_by_distance.py /tmp/cov/human --reference hg38 \
        --bins 0.05,0.15,0.3,0.6,1.0,1.5 --tsv /tmp/cov/human.tsv

The reference name defaults to the ``assembly`` field of the first manifest.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

CLASSES = ("cds", "utr", "intron", "intergenic")


# ----------------------------------------------------------------- Newick

def parse_newick(text: str):
    """Return (children, parent, length, name) tables for a Newick string.
    Nodes are integers; leaves carry a name."""
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    children: list[list[int]] = []
    parent: list[int] = []
    length: list[float] = []
    name: list[str] = []

    def new(p: int) -> int:
        children.append([])
        parent.append(p)
        length.append(0.0)
        name.append("")
        if p >= 0:
            children[p].append(len(children) - 1)
        return len(children) - 1

    root = new(-1)
    cur = root
    i, n = 0, len(text)
    buf = ""

    def flush(node: int, token: str) -> None:
        token = token.strip()
        if not token:
            return
        if ":" in token:
            nm, ln = token.split(":", 1)
            try:
                length[node] = float(ln)
            except ValueError:
                length[node] = 0.0
        else:
            nm = token
        nm = nm.strip().strip("'")
        if nm and not name[node]:
            name[node] = nm

    while i < n:
        ch = text[i]
        if ch == "(":
            cur = new(cur)
        elif ch == ",":
            flush(cur, buf)
            buf = ""
            cur = new(parent[cur])
        elif ch == ")":
            flush(cur, buf)
            buf = ""
            cur = parent[cur]
        elif ch in " \t\r\n":
            pass
        else:
            buf += ch
        i += 1
    flush(cur, buf)
    return children, parent, length, name


def patristic_from(text: str, ref: str) -> dict[str, float]:
    """Distance from leaf ``ref`` to every other leaf, in branch-length units."""
    children, parent, length, name = parse_newick(text)
    idx = {nm: i for i, nm in enumerate(name) if nm}
    if ref not in idx:
        # tolerate a suffix such as ".chr2L" or a trailing version
        cands = [nm for nm in idx if nm.split(".")[0] == ref]
        if not cands:
            raise KeyError(f"reference {ref!r} not a leaf of the tree")
        ref = cands[0]
    # distance from every node to the reference leaf: walk up from ref, then down
    up: dict[int, float] = {}
    node, d = idx[ref], 0.0
    while node >= 0:
        up[node] = d
        d += length[node]
        node = parent[node]
    out: dict[str, float] = {}

    def down(node: int, d: float) -> None:
        for c in children[node]:
            if c in up:
                continue
            dc = d + length[c]
            if not children[c] and name[c]:
                out[name[c]] = dc
            down(c, dc)

    for node, d in up.items():
        down(node, d)
    out.pop(ref, None)
    return out


# ------------------------------------------------------------- aggregation

def load_windows(paths: list[str]) -> list[tuple[str, dict, str]]:
    files: list[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(glob.glob(os.path.join(p, "*.manifest.json"))))
        else:
            files.append(p)
    out = []
    for f in files:
        with open(f) as fh:
            m = json.load(fh)
        nh = f[: -len(".manifest.json")] + ".nh"
        if not os.path.exists(nh):
            print(f"warning: no tree next to {f}, skipped", file=sys.stderr)
            continue
        with open(nh) as fh:
            text = fh.read()
        # UCSC writes one multi-line tree; Ensembl writes one tree per block, take the first
        tree = next((t for t in text.split(";") if t.strip()), "").strip() + ";"
        out.append((f, m, tree))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("inputs", nargs="+", help="fetch_window output directories or manifest files")
    ap.add_argument("--reference", default=None, help="reference leaf name (default: manifest assembly)")
    ap.add_argument("--bins", default="0.1,0.25,0.5,1.0,2.0",
                    help="upper edges of the distance bins, substitutions per site")
    ap.add_argument("--tsv", default=None, help="also write the per-informant-window rows here")
    ap.add_argument("--markdown", action="store_true", help="print the summary as a Markdown table")
    args = ap.parse_args()

    windows = load_windows(args.inputs)
    if not windows:
        print("no manifests found", file=sys.stderr)
        return 1
    ref = args.reference or windows[0][1].get("assembly")
    edges = [float(x) for x in args.bins.split(",")]
    labels = []
    lo = 0.0
    for e in edges:
        labels.append(f"{lo:g}-{e:g}")
        lo = e
    labels.append(f">{lo:g}")

    def bin_of(d: float) -> int:
        for k, e in enumerate(edges):
            if d < e:
                return k
        return len(edges)

    # per bin, per class: covered bases and total bases; plus pairs above 0.5
    cov = [[0.0, 0.0] for _ in range(len(labels) * len(CLASSES))]
    above = [[0, 0] for _ in range(len(labels) * len(CLASSES))]
    informants: list[set[str]] = [set() for _ in labels]
    rows = []
    n_windows = 0
    totals_all = {c: 0 for c in CLASSES}
    missing_leaf: set[str] = set()
    for f, m, tree in windows:
        try:
            dist = patristic_from(tree, ref)
        except KeyError as e:
            print(f"warning: {f}: {e}", file=sys.stderr)
            continue
        n_windows += 1
        totals = m["coverage"]["window_bases_by_class"]
        for c in CLASSES:
            totals_all[c] += totals.get(c, 0)
        table = m["coverage"]["species"]
        # informants in the tree but absent from every block are unaligned everywhere
        for sp in dist:
            if sp not in table:
                table = dict(table)
                table[sp] = {c: (0.0 if totals.get(c) else None) for c in CLASSES}
        for sp, frac in table.items():
            if sp not in dist:
                missing_leaf.add(sp)
                continue
            b = bin_of(dist[sp])
            informants[b].add(sp)
            row = {"window": os.path.basename(f)[: -len(".manifest.json")], "informant": sp,
                   "distance": round(dist[sp], 4), "bin": labels[b]}
            for ci, c in enumerate(CLASSES):
                tot = totals.get(c, 0)
                fr = frac.get(c)
                row[c] = fr
                if not tot or fr is None:
                    continue
                k = b * len(CLASSES) + ci
                cov[k][0] += fr * tot
                cov[k][1] += tot
                above[k][1] += 1
                if fr > 0.5:
                    above[k][0] += 1
            rows.append(row)

    if missing_leaf:
        print(f"warning: {len(missing_leaf)} MAF sources not in the tree, ignored: "
              f"{sorted(missing_leaf)[:8]}", file=sys.stderr)

    print(f"reference {ref}; {n_windows} windows; reference bases by class: "
          + ", ".join(f"{c} {totals_all[c]:,}" for c in CLASSES))
    hdr = ["distance (subst/site)", "informants"] + [f"{c} cov" for c in CLASSES] + [f"{c} >0.5" for c in CLASSES]
    lines = []
    for b, lab in enumerate(labels):
        if not informants[b]:
            continue
        cells = [lab, str(len(informants[b]))]
        for ci in range(len(CLASSES)):
            k = b * len(CLASSES) + ci
            cells.append(f"{cov[k][0] / cov[k][1]:.3f}" if cov[k][1] else "-")
        for ci in range(len(CLASSES)):
            k = b * len(CLASSES) + ci
            cells.append(f"{above[k][0] / above[k][1]:.2f}" if above[k][1] else "-")
        lines.append(cells)
    if args.markdown:
        print("| " + " | ".join(hdr) + " |")
        print("|" + "---|" * len(hdr))
        for cells in lines:
            print("| " + " | ".join(cells) + " |")
    else:
        print("\t".join(hdr))
        for cells in lines:
            print("\t".join(cells))
    if args.tsv:
        with open(args.tsv, "w") as fh:
            keys = ["window", "informant", "distance", "bin", *CLASSES]
            fh.write("\t".join(keys) + "\n")
            for r in rows:
                fh.write("\t".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
        print(f"wrote {len(rows)} rows to {args.tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
