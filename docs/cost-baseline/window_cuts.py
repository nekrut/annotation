#!/usr/bin/env python3
"""Count the reference loci that a windowed AUGUSTUS run cuts.

``run_augustus_windows.sh`` splits one sequence of length L into K windows of
ceil(L/K) bp and runs AUGUSTUS separately in each with --genemodel=complete,
so a gene whose CDS spans a window boundary cannot be predicted whole.  This
script reproduces the K-1 cut coordinates from L and K, applies the same
reference filter as ``benchmark/score.py`` (complete protein-coding
transcripts; pseudogenes and V/D/J/C segments dropped; loci are the GFF3
gene grouping with CDS merged over isoforms) and counts the distinct loci
whose merged-CDS span contains a cut, plus, for context, the loci whose CDS
span lies within a given distance of one.  The crossing count is what the
run cannot get right; the near count is not a bound on anything, because
--genemodel=complete can change a prediction anywhere in a window.

usage: window_cuts.py --reference GFF[.gz] --seqid NC_000021.9 --length 46709983 --windows 8 [--near 10000,50000]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "benchmark"))
import score  # noqa: E402


def cuts(length, k):
    step = (length + k - 1) // k
    if k < 1 or length < k or (k - 1) * step >= length:
        # Same rule as run_augustus_windows.sh: ceil(length/k) must leave no
        # empty window (length=10, k=6 gives step 2 and an empty window 5).
        raise ValueError("windows=%d leaves an empty window for %d bp (step %d)" % (k, length, step))
    out = []
    for i in range(k - 1):
        e = min((i + 1) * step, length)
        out.append((i, e, e + 1))  # window i ends at e; window i+1 starts at e+1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--seqid", required=True)
    ap.add_argument("--length", type=int, required=True, help="sequence length in bp, as the runner measures it")
    ap.add_argument("--windows", type=int, required=True)
    ap.add_argument("--near", default="10000,50000", help="comma-separated distances in bp for the context counts")
    args = ap.parse_args()

    ann = score.load_gff(args.reference)
    chains = score.select_transcripts(ann.chains(), ann)
    loci = score.loci_of(ann, chains, {args.seqid})
    spans = {gid: (blocks[0][0], blocks[-1][1]) for gid, (seqid, strand, blocks, tids) in loci.items()}
    boundaries = cuts(args.length, args.windows)
    near = [int(x) for x in args.near.split(",") if x]

    print("seqid=%s length=%d windows=%d step=%d reference_loci=%d" % (
        args.seqid, args.length, args.windows, (args.length + args.windows - 1) // args.windows, len(loci)))
    print("cut\tlast_bp_of_window\tfirst_bp_of_next\tloci_crossing\tgene_ids" + "".join("\tloci_within_%d" % d for d in near))
    crossing_all, near_all = set(), {d: set() for d in near}
    for i, e, s in boundaries:
        crossing = sorted(g for g, (a, b) in spans.items() if a <= e and b >= s)
        crossing_all.update(crossing)
        # score.py namespaces ids as seqid\0strand\0id; show the id alone.
        names = [g.split("\0")[-1] for g in crossing]
        row = [str(i), str(e), str(s), str(len(crossing)), ",".join(names) or "-"]
        for d in near:
            within = {g for g, (a, b) in spans.items() if a <= e + d and b >= s - d}
            near_all[d].update(within)
            row.append(str(len(within)))
        print("\t".join(row))
    print("distinct_loci_crossing_any_cut=%d" % len(crossing_all)
          + "".join(" distinct_loci_within_%d=%d" % (d, len(near_all[d])) for d in near))
    return 0


if __name__ == "__main__":
    sys.exit(main())
