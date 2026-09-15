#!/usr/bin/env python3
"""Decoder-based check that admitted chains have a finite legal numerator.

The admission audit (`admission.py`) admits a chain only when its rows are
in frame, every positive gap is at least m, the spliced CDS starts with an
initiator, ends with a stop and has no in-frame stop or ambiguous base. By
construction such a chain is a legal path of candidate A's grammar. This
script checks that claim with the grammar itself on a sample: it cuts each
sampled admitted gene's window from the checksummed FASTA, orients it,
rewards exactly the annotated CDS and intron bases (`Scores.favour`), and
runs the delayed-entry decoder with chunk seams. The chain has a finite
numerator when the partition and the Viterbi score are finite and the
Viterbi chain reproduces the annotated spliced CDS and intron intervals.

Pure-Python decoding is slow, so the sample is bounded by count and span;
the summary lists what was sampled and any failure with its reason.

    python3 -m model.labels.numerator_check --species X --gff G --fasta F \\
        --manifest model/labels/manifests/X.manifest.tsv.gz --per-species 5 --max-span 4000
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import sys
import time
from math import isfinite

from model.grammar import DelayedEntryDecoder, DurationMixture, Scores, TABLES
from model.labels.admission import load_gff_rows, revcomp, _iter_fasta

FLANK = 10
DUR = DurationMixture(m=20, pi=((0.5, 0.3, 0.2),) * 3, q=((0.9, 0.99, 0.999),) * 3)


def oriented_chain(t, seq):
    """(window, cds_ranges, intron_ranges) in oriented half-open coordinates
    with FLANK unrewarded bases on each side; adjacent rows are merged."""
    s, e = t.span
    ws, we = max(1, s - FLANK), min(len(seq), e + FLANK)
    window = seq[ws - 1:we].upper()
    rows = sorted(t.rows, key=lambda r: r.start)
    ivs = [(r.start - ws, r.end - ws + 1) for r in rows]
    if t.strand == "-":
        window = revcomp(window)
        n = len(window)
        ivs = [(n - b, n - a) for a, b in ivs][::-1]
    merged = []
    for a, b in ivs:
        if merged and merged[-1][1] == a:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    introns = [(a[1], b[0]) for a, b in zip(merged, merged[1:])]
    return window, merged, introns


def check(t, seq, table, seam=1024):
    window, cds, introns = oriented_chain(t, seq)
    n = len(window)
    sc = Scores.zeros(n).favour(cds, introns)
    dec = DelayedEntryDecoder(TABLES[table], DUR)
    seams = list(range(seam, n, seam))
    t0 = time.time()
    z = dec.partition(window, sc, seams=seams)
    v, chains = dec.viterbi(window, sc, seams=seams)
    dt = time.time() - t0
    out = {"transcript": t.tid, "seqid": t.seqid, "strand": t.strand, "window_bp": n,
           "introns": len(introns), "log_partition": z, "viterbi": v, "seconds": round(dt, 2)}
    want = "".join(window[a:b] for a, b in cds)
    if not (isfinite(z) and isfinite(v)):
        out["failure"] = "non-finite score"
    elif len(chains) != 1:
        out["failure"] = "%d chains decoded" % len(chains)
    else:
        c = chains[0]
        got_cds = [(s.start, s.end) for s in c.cds()]
        got_introns = [(s.start, s.end) for s in c.introns()]
        if c.spliced(window) != want or got_introns != introns or c.partial_5 or c.partial_3:
            out["failure"] = "chain differs: cds %s introns %s partial %s/%s" % (got_cds, got_introns, c.partial_5, c.partial_3)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--per-species", type=int, default=5)
    ap.add_argument("--max-span", type=int, default=4000)
    ap.add_argument("--min-introns", type=int, default=0)
    ap.add_argument("--table", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20260915)
    a = ap.parse_args(argv)
    admitted = []
    with gzip.open(a.manifest, "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            if r["status"] == "admitted" and int(r["cds_end"]) - int(r["cds_start"]) + 1 <= a.max_span \
                    and int(r["rows"]) - 1 >= a.min_introns:
                admitted.append((r["seqid"], r["strand"], r["transcript"]))
    rng = random.Random(a.seed)
    sample = set(rng.sample(admitted, min(a.per_species, len(admitted))))
    transcripts, _, _ = load_gff_rows(a.gff)
    wanted = {(t.seqid, t.strand, t.tid): t for t in transcripts if (t.seqid, t.strand, t.tid) in sample}
    results = []
    need = {k[0] for k in wanted}
    for name, seq in _iter_fasta(a.fasta):
        if name not in need:
            continue
        for k, t in wanted.items():
            if k[0] == name:
                results.append(check(t, seq, a.table))
        need.discard(name)
        if not need:
            break
    summary = {"species": a.species, "eligible_admitted": len(admitted), "sampled": len(results),
               "max_span": a.max_span, "failures": [r for r in results if "failure" in r],
               "results": results}
    print(json.dumps(summary, indent=1))
    return 1 if summary["failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
