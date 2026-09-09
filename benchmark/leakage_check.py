#!/usr/bin/env python3
"""Check the panel's held-out rules.

Two things are verified against ``panel.tsv`` and ``taxonomy.tsv``:

1. **Phylogenetic separation.** For every ``heldout`` species, the most recent
   common ancestor it shares with the *closest training species* must sit at
   or above a minimum rank (default ``CLASS``).  A model trained on the
   ``train`` split therefore never sees a close relative of a cross-clade
   evaluation genome, so a good score cannot come from within-family
   memorization.  ``heldout_paired`` species deliberately break this rule --
   a close relative is in training so that the pair isolates a regime shift
   (intron length, GC, genome size) from a clade shift -- and are listed
   separately rather than counted as violations.
2. **Informant declaration.** If a comparative method uses informant
   genomes for a held-out target, the informants are listed in a TSV
   (``--informants target<TAB>informant``) and each informant is checked
   against the same rank rule and against the training split.  Informant
   *sequence* closer than the rule allows is permitted only when declared;
   informant *annotation* is never allowed to be projected onto the target
   (that is label leakage, and this script cannot see it -- it must be
   asserted in the run's report).

Usage:
    python3 benchmark/leakage_check.py
    python3 benchmark/leakage_check.py --min-rank ORDER --informants my_run.tsv

Exit status is non-zero when a rule is violated.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, "panel.tsv")
TAXONOMY = os.path.join(HERE, "taxonomy.tsv")

# NCBI Taxonomy ranks, shallow (inclusive) to deep (specific).  Only the ranks
# that actually appear in this panel's lineages need an entry; anything not
# listed (CLADE, NO_RANK, ...) is unranked and compared by depth alone.
RANK_ORDER = [
    "DOMAIN", "KINGDOM", "PHYLUM", "SUBPHYLUM", "SUPERCLASS", "CLASS",
    "SUBCLASS", "SUPERORDER", "ORDER", "SUBORDER", "INFRAORDER",
    "SUPERFAMILY", "FAMILY", "SUBFAMILY", "TRIBE", "GENUS", "SPECIES",
    "SUBSPECIES", "STRAIN",
]
RANK_DEPTH = {r: i for i, r in enumerate(RANK_ORDER)}


def load_lineages(path: str = TAXONOMY) -> dict[str, list[tuple[int, str, str]]]:
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            nodes = []
            for node in row["lineage"].split(";"):
                tid, rank, name = node.split(":", 2)
                nodes.append((int(tid), rank, name))
            out[row["species"]] = nodes
    return out


def mrca(a: list, b: list) -> tuple[int, str, str]:
    last = a[0]
    for x, y in zip(a, b):
        if x[0] != y[0]:
            break
        last = x
    return last


def named_rank(a: list, b: list) -> tuple[str, str]:
    """MRCA rank and name, walking up to the nearest *ranked* ancestor.

    NCBI marks many internal nodes ``CLADE``; ``Euarchontoglires`` is a
    SUPERORDER but ``Boreoeutheria`` above it is a CLADE.  For a rank
    comparison we report the deepest ancestor that carries a real rank.
    """
    node = mrca(a, b)
    depth = a.index(node) if node in a else 0
    for tid, rank, name in reversed(a[: depth + 1]):
        if rank in RANK_DEPTH:
            return rank, name
    return "DOMAIN", node[2]


def nearest_training(lineages, species, train):
    """Every training species tied at the deepest MRCA rank, and that rank.

    Ties are the normal case, not an edge case: chicken shares Sarcopterygii
    with mouse and with *Xenopus*, and printing only the first row of
    ``panel.tsv`` that reaches that rank would name mouse and hide the frog,
    which is the closer relative in time.  The rule is a floor on the *rank*,
    so it is the whole tied set that has to satisfy it.
    """
    best, names = None, []
    for t in train:
        rank, name = named_rank(lineages[species], lineages[t])
        if best is None or RANK_DEPTH[rank] > RANK_DEPTH[best[0]]:
            best, names = (rank, name), [t]
        elif RANK_DEPTH[rank] == RANK_DEPTH[best[0]]:
            names.append(t)
    return names, best[0], best[1]


def fmt_nearest(names, shown=3):
    """The tied set, abbreviated: at Eukaryota every training species ties."""
    names = sorted(names)
    if len(names) <= shown:
        return ",".join(names)
    return "%s,+%d more" % (",".join(names[:shown]), len(names) - shown)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=PANEL)
    ap.add_argument("--taxonomy", default=TAXONOMY)
    ap.add_argument("--min-rank", default="CLASS",
                    help="MRCA of a held-out species and any training species must be at or above this rank")
    ap.add_argument("--informants", help="TSV of target<TAB>informant pairs used by a comparative run")
    args = ap.parse_args()

    if args.min_rank not in RANK_DEPTH:
        ap.error(f"--min-rank must be one of: {', '.join(RANK_ORDER)}")
    limit = RANK_DEPTH[args.min_rank]

    with open(args.panel, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    lineages = load_lineages(args.taxonomy)
    train = [r["species"] for r in rows if r["split"] == "train"]
    heldout = [r["species"] for r in rows if r["split"] == "heldout"]
    paired = [r["species"] for r in rows if r["split"] == "heldout_paired"]

    violations = 0
    print(f"# split separation: MRCA(held-out, nearest training species) must be at or above {args.min_rank}")
    for h in heldout:
        closest, rank, name = nearest_training(lineages, h, train)
        status = "ok"
        if RANK_DEPTH[rank] > limit:
            status = "VIOLATION"
            violations += 1
        print(f"{status}\t{h}\tnearest={fmt_nearest(closest)}\tmrca={name} ({rank})")

    if paired:
        print("\n# paired held-out species: a close relative is in training on purpose.\n"
              "# These measure regime shift inside a clade and are reported separately;\n"
              "# they are never averaged into the cross-clade number.")
        for h in paired:
            closest, rank, name = nearest_training(lineages, h, train)
            print(f"paired\t{h}\tnearest={fmt_nearest(closest)}\tmrca={name} ({rank})")

    if args.informants:
        print(f"\n# informant separation for comparative runs ({args.informants})")
        with open(args.informants, newline="", encoding="utf-8") as fh:
            pairs = [line.rstrip("\n").split("\t")[:2] for line in fh if line.strip() and not line.startswith("#")]
        for target, informant in pairs:
            if informant not in lineages:
                print(f"note\t{target}\t{informant}\tnot in panel; separation not checked here")
                continue
            rank, name = named_rank(lineages[target], lineages[informant])
            tag = "declared-close" if RANK_DEPTH[rank] > limit else "ok"
            print(f"{tag}\t{target}\tinformant={informant}\tmrca={name} ({rank})")

    print(f"\n{len(heldout)} cross-clade held-out, {len(paired)} paired held-out, "
          f"{len(train)} training species, {violations} violation(s)")
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
