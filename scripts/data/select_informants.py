#!/usr/bin/env python3
"""Select up to K-1 informant species for one target from a multiz tree:
the fixed, annotation-free rule that proposal section 2 requires of B's
loader ("up to seven distinct informant species plus the target, K <= 8,
... alignment coverage and tree diversity without annotation, a fixed rule
frozen on training development chromosomes").

Rule (deterministic, one selection per track, reused for every chunk and
every compared tree arm):

1. Eligibility. A leaf is eligible when it is not the target, is not in an
   excluded relation of ``informant_membership.tsv`` (default: ``drop``, the
   held-out panel species that a training run must remove), and passes the
   coverage gate. The gate is ``coverage >= --min-coverage`` when a coverage
   table is given, otherwise the patristic distance to the target must be
   ``<= --horizon`` substitutions per site (default 1.0, the point past
   which non-coding sequence is essentially unaligned in
   ``docs/data-sources.md`` section 6). The horizon is a proxy for
   coverage until ``fetch_window.py`` coverage tables exist on the reserved
   development chromosomes; a run with ``--coverage`` supersedes it.
2. One assembly per taxon. Leaves sharing an NCBI taxid (two assemblies of
   one species) collapse to one: the higher coverage, else the later
   assembly name (``hg38`` over ``hg19``), so duplicate assemblies never
   count as independent taxa.
3. Greedy phylogenetic diversity. Starting from the target, repeatedly add
   the eligible leaf whose path to the already selected subtree adds the
   most branch length (Faith's PD of the induced subtree; the greedy choice
   is optimal for PD on a tree, Steel 2005 / Pardi and Goldman 2005). Ties
   go to the leaf closer to the target, then to the name. The selection
   for K is the first K-1 rows, so the K=16 sensitivity arm is the same
   list read further.

Coverage table: TSV with a header and columns ``leaf`` and ``coverage``
(fraction of target bases, or of target CDS bases, at which the informant
has an aligned non-gap base). ``coverage_by_distance.py`` and the
``fetch_window.py`` manifests supply the per-species numbers.

Standard library only, Python 3.11. Usage::

    python3 scripts/data/select_informants.py --nh mm39.35way.nh --reference mm39 \
        --membership scripts/data/informant_membership.tsv --track mm39/multiz35way
    python3 scripts/data/select_informants.py --nh hg38.100way.nh --reference hg38 \
        --membership ... --track hg38/multiz100way --exclude-relations '' --k 15

For the evaluation-only tracks (hg38 100-way, galGal6 77-way) pass
``--exclude-relations ''`` so held-out informants stay eligible; the
``relation`` column of the output is then the declaration the benchmark
(section 3.2) requires for such a run.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tree_composition import parse_newick, leaves, patristic  # noqa: E402


def read_membership(path: str, track: str) -> dict[str, dict]:
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["track"] == track]
    if not rows:
        raise SystemExit(f"no rows for track {track!r} in {path}")
    return {r["leaf"]: r for r in rows}


def read_coverage(path: str) -> dict[str, float]:
    with open(path, newline="") as fh:
        return {r["leaf"]: float(r["coverage"]) for r in csv.DictReader(fh, delimiter="\t")}


def version_key(name: str):
    m = re.search(r"(\d+)$", name)
    return (int(m.group(1)) if m else -1, name)


def select(root, reference: str, k: int, membership, coverage, min_coverage, horizon, exclude):
    by_name = {n.name: n for n in leaves(root)}
    if reference not in by_name:
        raise SystemExit(f"reference {reference!r} is not a leaf of the tree")
    ref = by_name[reference]
    dist = {name: patristic(ref, n) for name, n in by_name.items()}

    excluded = []  # (leaf, reason)
    eligible = []
    for name in sorted(by_name):
        if name == reference:
            continue
        rel = membership.get(name, {}).get("relation", "") if membership else ""
        if rel in exclude:
            excluded.append((name, f"relation={rel}"))
            continue
        if coverage is not None:
            c = coverage.get(name)
            if c is None:
                excluded.append((name, "no coverage row"))
                continue
            if c < min_coverage:
                excluded.append((name, f"coverage {c:.3f} < {min_coverage}"))
                continue
        elif dist[name] > horizon:
            excluded.append((name, f"distance {dist[name]:.3f} > horizon {horizon}"))
            continue
        eligible.append(name)

    # one assembly per taxid
    by_taxid: dict[str, list[str]] = {}
    for name in eligible:
        taxid = membership.get(name, {}).get("taxid", "") if membership else ""
        by_taxid.setdefault(taxid or f"leaf:{name}", []).append(name)
    kept = []
    for taxid, names in by_taxid.items():
        if len(names) == 1:
            kept.append(names[0])
            continue
        if coverage is not None:
            best = max(names, key=lambda n: (coverage[n], version_key(n)))
        else:
            best = max(names, key=version_key)
        kept.append(best)
        for n in names:
            if n != best:
                excluded.append((n, f"duplicate assembly of taxid {taxid}, kept {best}"))
    eligible = sorted(kept)

    # greedy PD: covered = node ids already inside the induced subtree
    covered = set()
    n = ref
    while n is not None:
        covered.add(id(n))
        n = n.parent

    def gain(name):
        g, n = 0.0, by_name[name]
        while id(n) not in covered:
            g += n.length
            n = n.parent
        return g

    chosen, cumulative = [], 0.0
    remaining = set(eligible)
    while remaining and len(chosen) < k:
        best = max(remaining, key=lambda nm: (gain(nm), -dist[nm], nm))
        # deterministic tie-break: larger gain, then closer, then name (reverse
        # so that max() picks the lexicographically first)
        best = min((nm for nm in remaining if gain(nm) == gain(best) and dist[nm] == dist[best]),
                   default=best)
        g = gain(best)
        cumulative += g
        chosen.append((best, g, cumulative))
        remaining.discard(best)
        n = by_name[best]
        while id(n) not in covered:
            covered.add(id(n))
            n = n.parent
    return chosen, dist, eligible, excluded


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--nh", required=True, help="Newick tree with branch lengths")
    ap.add_argument("--reference", required=True, help="target leaf name")
    ap.add_argument("--k", type=int, default=7, help="informants to select (K-1); default 7")
    ap.add_argument("--membership", help="informant_membership.tsv from informant_audit.py")
    ap.add_argument("--track", help="track name in the membership table, e.g. mm39/multiz35way")
    ap.add_argument("--exclude-relations", default="drop",
                    help="comma-separated relations to exclude (default: drop; '' keeps all)")
    ap.add_argument("--coverage", help="TSV with columns leaf, coverage")
    ap.add_argument("--min-coverage", type=float, default=0.5)
    ap.add_argument("--horizon", type=float, default=1.0,
                    help="max patristic distance without a coverage table")
    ap.add_argument("--out", default="-", help="output TSV ('-' for stdout)")
    ap.add_argument("--no-header", action="store_true")
    a = ap.parse_args(argv)

    if (a.membership is None) != (a.track is None):
        ap.error("--membership and --track go together")
    membership = read_membership(a.membership, a.track) if a.membership else {}
    coverage = read_coverage(a.coverage) if a.coverage else None
    exclude = {r for r in a.exclude_relations.split(",") if r}
    text = Path(a.nh).read_text()
    unit = ":" not in text
    if unit:
        # sacCer3/multiz7way ships a topology only, with whitespace between
        # sister labels: fall back to unit branch lengths and disable the
        # distance gate, which would be meaningless in those units.
        text = re.sub(r"(?<=[A-Za-z0-9_.\-)])\s+(?=[A-Za-z0-9_.(])", ",", text.strip())
        print(f"warning: {a.nh} has no branch lengths; using unit lengths and no "
              "distance gate (PD counts edges)", file=sys.stderr)
        a.horizon = float("inf")
    root = parse_newick(text)
    if unit:
        stack = [root]
        while stack:
            n = stack.pop(); n.length = 0.0 if n is root else 1.0; stack.extend(n.children)
    chosen, dist, eligible, excluded = select(
        root, a.reference, a.k, membership, coverage, a.min_coverage, a.horizon, exclude)

    gate = (f"coverage>={a.min_coverage}" if coverage is not None
            else "topology-only" if unit else f"distance<={a.horizon}")
    out = sys.stdout if a.out == "-" else open(a.out, "w", newline="")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    if not a.no_header:
        w.writerow(["track", "reference", "rank", "leaf", "scientific_name", "taxid", "relation",
                    "distance", "pd_gain", "pd_cumulative", "coverage", "gate"])
    for rank, (name, g, cum) in enumerate(chosen, 1):
        m = membership.get(name, {})
        cov = "" if coverage is None else f"{coverage[name]:.4f}"
        w.writerow([a.track or Path(a.nh).name, a.reference, rank, name,
                    m.get("scientific_name", ""), m.get("taxid", ""), m.get("relation", ""),
                    f"{dist[name]:.4f}", f"{g:.4f}", f"{cum:.4f}", cov, gate])
    if out is not sys.stdout:
        out.close()

    n_leaves = len(dist) - 1
    print(f"{a.track or a.nh}: {n_leaves} candidate leaves, {len(eligible)} eligible after "
          f"{gate} and taxid collapse, {len(excluded)} excluded, {len(chosen)} selected "
          f"(PD {chosen[-1][2]:.3f} subst/site)" if chosen else "nothing selected", file=sys.stderr)
    for name, why in excluded:
        if not why.startswith("distance"):
            print(f"  excluded {name}: {why}", file=sys.stderr)
    n_far = sum(1 for _, why in excluded if why.startswith("distance"))
    if n_far:
        print(f"  excluded {n_far} leaves beyond the horizon", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
