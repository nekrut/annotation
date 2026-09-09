#!/usr/bin/env python3
"""Assemble scored species into the section 4.8 table and its two aggregates.

Takes the JSON rows ``score.py`` writes, joins them to ``panel.tsv`` for the
clade and split of each species, and emits the primary result table plus the
cross-clade and regime-shift scores.  Cost columns (section 4.7) and BUSCO
(4.6) are not computed here; supply them with ``--cost`` as a TSV whose first
column is the species and whose remaining columns are merged into the row.

Refuses to report an aggregate over an incomplete panel unless
``--partial`` is given, because the aggregates in section 4.8 are unweighted
means over a fixed species set and a mean over a subset is a different
number wearing the same name.

Usage:
    python3 report.py results/*.json [--cost cost.tsv] [--markdown]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
import sys

PANEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "panel.tsv")

# Section 4.8: 20 species x these columns is the primary result.
COLUMNS = [
    ("nucleotide F1", ("nucleotide", "f1")),
    ("nucleotide MCC", ("nucleotide", "mcc")),
    ("exon F1", ("exon", "all", "f1")),
    ("donor F1", ("splice", "donor", "f1")),
    ("acceptor F1", ("splice", "acceptor", "f1")),
    ("transcript F1", ("transcript", "f1")),
    ("locus F1", ("locus", "f1")),
    ("fusion", ("locus", "fusion")),
    ("split", ("locus", "split")),
    ("start F1", ("codon", "start", "f1")),
    ("stop F1", ("codon", "stop", "f1")),
]

# Section 4.8: Tetrahymena is reported, never ranked (section 2.3).
NEVER_RANKED = {"Tetrahymena_thermophila"}


def dig(row, path):
    v = row
    for k in path:
        if v is None:
            return None
        v = v.get(k)
    return v


def load_panel():
    meta = {}
    with open(PANEL, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            meta[r["species"]] = r
    return meta


def aggregate(rows, meta, split, key, exclude=NEVER_RANKED):
    vals = [dig(r, key) for s, r in sorted(rows.items())
            if meta.get(s, {}).get("split") == split and s not in exclude]
    vals = [v for v in vals if isinstance(v, (int, float))]
    return round(statistics.fmean(vals), 5) if vals else None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", nargs="+", help="JSON rows from score.py")
    ap.add_argument("--cost", help="TSV of measured cost columns, keyed by species")
    ap.add_argument("--markdown", action="store_true", help="markdown instead of TSV")
    ap.add_argument("--partial", action="store_true",
                    help="report aggregates over the species present, not the panel")
    args = ap.parse_args()

    meta = load_panel()
    rows = {}
    seen = {}
    dupes = []
    for path in args.results:
        r = json.load(open(path))
        sp = r["species"]
        if sp in seen:
            dupes.append((sp, seen[sp], path))
        seen[sp] = path
        rows[sp] = r
    # Two results for one species is a submission with two runs in it, and
    # keeping the last silently reports one of them: an ablation scored
    # alongside the real run would vanish here without a word.
    if dupes:
        for sp, first, second in dupes:
            print("%s appears twice: %s and %s" % (sp, first, second),
                  file=sys.stderr)
        print("one result per species; drop or merge the duplicates",
              file=sys.stderr)
        return 2
    cost = {}
    cost_cols = []
    if args.cost:
        with open(args.cost, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            head = next(rd)
            cost_cols = head[1:]
            for line in rd:
                cost[line[0]] = dict(zip(cost_cols, line[1:]))

    # A row whose stop-codon convention was never established is 3 bp of
    # doubt at the 3' end of every chain in it (docs/benchmark.md 4.5), and
    # the aggregate would hide that behind nineteen settled ones.
    assumed = sorted(sp for sp, r in rows.items()
                     if r.get("codon", {}).get("stop_codon_convention_source")
                     == "assumed")
    if assumed:
        print("warning: stop-codon convention never established for %s. "
              "Their terminal-exon, single-exon, stop-codon and exact-"
              "transcript scores rest on the GFF3 default; re-score with "
              "--genome." % ", ".join(assumed), file=sys.stderr)

    unknown = sorted(set(rows) - set(meta))
    if unknown:
        print("not in panel.tsv: %s" % ", ".join(unknown), file=sys.stderr)
        return 2
    missing = sorted(set(meta) - set(rows))
    if missing and not args.partial:
        print("panel is incomplete, %d species missing: %s"
              % (len(missing), ", ".join(missing)), file=sys.stderr)
        print("re-run with --partial to report over what is present", file=sys.stderr)
        return 2

    header = ["species", "clade", "split"] + [c for c, _ in COLUMNS] + cost_cols
    table = [header]
    for species in sorted(rows, key=lambda s: (meta[s]["clade"], s)):
        r = rows[species]
        line = [species, meta[species]["clade"], meta[species]["split"]]
        for _, path in COLUMNS:
            v = dig(r, path)
            line.append("" if v is None else str(v))
        line += [cost.get(species, {}).get(c, "") for c in cost_cols]
        table.append(line)

    if args.markdown:
        print("| " + " | ".join(table[0]) + " |")
        print("|" + "---|" * len(table[0]))
        for line in table[1:]:
            print("| " + " | ".join(line) + " |")
    else:
        w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        w.writerows(table)

    print()
    print("# aggregates (unweighted means; %s excluded, section 2.3)"
          % ", ".join(sorted(NEVER_RANKED)))
    if missing:
        print("# PARTIAL: %d of %d panel species scored" % (len(rows), len(meta)))
    for label, split in (("cross-clade score", "heldout"),
                         ("regime-shift score", "heldout_paired")):
        parts = ["%s=%s" % (name, aggregate(rows, meta, split, path))
                 for name, path in COLUMNS if name.endswith("F1")]
        print("%s (%s): %s" % (label, split, "  ".join(parts)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
