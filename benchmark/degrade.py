#!/usr/bin/env python3
"""Make a synthetically degraded copy of a reference annotation.

An identity run (scoring a reference against itself) exercises every code
path in ``score.py`` but pins only the fixed points: every metric is 1.0 by
construction, so a scorer that silently dropped half the annotation from
*both* sides would still pass.  A *degraded* copy is the cheap control that
identity cannot be: the perturbation is known exactly, so the direction and
rough size of every metric's response is predictable before the run, and any
metric that does not move -- or moves the wrong way -- is a defect in the
scorer, not in the input.

Two independent, seeded perturbations, both applied to the reference GFF3:

``--drop-transcripts P``
    With probability ``P`` per transcript, delete the transcript feature and
    every child of it.  Expected response: transcript and locus sensitivity
    fall by about ``P``, precision stays at 1.0 (nothing wrong is added), and
    the locus level falls less than the transcript level in any species whose
    loci carry several isoforms, because a locus survives while one of its
    transcripts does.

``--shift-cds P``
    With probability ``P`` per CDS segment, move that segment's **downstream
    boundary in transcription order** 3 bp further downstream: the ``end``
    field on the ``+`` strand, the ``start`` field on the ``-`` strand.  Both
    sensitivity and precision fall.  On a terminal segment this is a wrong
    stop codon and a wrong terminal exon; on an internal one it is a wrong
    donor site (``+``) or acceptor site (``-``) and two wrong exons.  Because
    the draw is per segment and not per transcript, a transcript with ``n``
    CDS segments survives intact with probability ``(1-P)**n``, so the
    transcript level falls much faster than the nucleotide level on any
    species with many exons per transcript.  That gap is the point: it is
    what makes a per-transcript metric worth reporting next to a per-base one.

The 3 bp shift is deliberately in frame.  An out-of-frame shift would also
break the CDS phase, and a scorer that reads phase could catch it for the
wrong reason; 3 bp keeps the chain translatable and forces the boundary
metrics to do the work.

Strand matters and is the reason this script exists in the repository rather
than as a shell one-liner.  A shift applied to the ``end`` field regardless
of strand moves the 3' boundary of ``+`` strand genes and the *5'* boundary
of ``-`` strand ones, which perturbs start codons on roughly half the
annotation while claiming to perturb stop codons.  In a panel that spans
19.5% to 48.5% GC and 12 Mb to 2.2 Gb the strand split is close to even
everywhere, so the mistake is invisible in the aggregate and changes what the
run means.

Sequence bounds are read from the ``##sequence-region`` pragmas and from
RefSeq ``region`` features, so a shift is never allowed to run off the end of
a sequence.  Ensembl GFF3 has the pragmas; a file with neither is degraded
with the bound check skipped and a warning on stderr.

Nothing is written except the output GFF3 and a JSON summary of what was
actually done, which is what a run's declaration should cite.  Standard
library only.

Usage:
    python3 benchmark/degrade.py --reference REF.gff.gz --out DEGRADED.gff3
    python3 benchmark/degrade.py --reference REF.gff.gz --out D.gff3 \
        --drop-transcripts 0.10 --shift-cds 0.10 --seed 20260909
    python3 benchmark/degrade.py --self-test
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import sys

TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "ncRNA", "rRNA", "tRNA", "snRNA", "snoRNA",
    "lnc_RNA", "miRNA", "primary_transcript", "guide_RNA", "RNase_P_RNA",
    "SRP_RNA", "RNase_MRP_RNA", "telomerase_RNA", "antisense_RNA",
    "V_gene_segment", "C_gene_segment", "D_gene_segment", "J_gene_segment",
}


def attr(attrs: str, key: str) -> str:
    """Value of one GFF3 attribute, or ``''``."""
    for kv in attrs.split(";"):
        kv = kv.strip()
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return ""


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def degrade(lines, drop_p: float, shift_p: float, seed: int):
    """Yield degraded GFF3 lines; return a counts dict via ``.summary``.

    ``lines`` is any iterable of GFF3 text lines.  Two independent random
    streams are used, one per perturbation, so changing one rate does not
    reshuffle the draws of the other and two runs that differ only in
    ``--shift-cds`` drop exactly the same transcripts.
    """
    drop_rnd = random.Random(f"{seed}:drop")
    shift_rnd = random.Random(f"{seed}:shift")
    seq_len: dict[str, int] = {}
    dropped: dict[str, bool] = {}
    out: list[str] = []
    counts = {
        "transcripts_seen": 0, "transcripts_dropped": 0,
        "child_lines_dropped": 0, "cds_seen": 0, "cds_shifted": 0,
        "cds_shift_refused_at_sequence_end": 0,
        "cds_shifted_plus": 0, "cds_shifted_minus": 0,
        "sequence_lengths_known": 0,
    }

    for line in lines:
        if line.startswith("#"):
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        seq_len[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
            out.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            out.append(line)
            continue
        seqid, _, ftype, start, end, _, strand, _, attrs = fields[:9]

        if ftype == "region":
            # RefSeq writes one ``region`` feature per sequence, spanning it.
            try:
                seq_len.setdefault(seqid, int(end))
            except ValueError:
                pass
            out.append(line)
            continue

        if ftype in TRANSCRIPT_TYPES:
            tid = attr(attrs, "ID")
            counts["transcripts_seen"] += 1
            gone = drop_rnd.random() < drop_p
            if tid:
                dropped[tid] = gone
            if gone:
                counts["transcripts_dropped"] += 1
                continue
            out.append(line)
            continue

        # The shift is drawn for *every* CDS line, including one that is about
        # to be deleted with its transcript.  Drawing only for the survivors
        # would make the shift stream a function of the drop rate, so a
        # drop-only run and a combined run at the same seed would not shift the
        # same segments and the two could not be compared -- which is the whole
        # point of running them.
        shift_this = ftype == "CDS" and shift_rnd.random() < shift_p

        parent = attr(attrs, "Parent").split(",")[0]
        if parent and dropped.get(parent):
            counts["child_lines_dropped"] += 1
            continue

        if ftype == "CDS":
            counts["cds_seen"] += 1
            if shift_this:
                limit = seq_len.get(seqid)
                if strand == "-":
                    new_start = int(start) - 3
                    if new_start < 1:
                        counts["cds_shift_refused_at_sequence_end"] += 1
                    else:
                        fields[3] = str(new_start)
                        counts["cds_shifted"] += 1
                        counts["cds_shifted_minus"] += 1
                        line = "\t".join(fields) + "\n"
                else:
                    new_end = int(end) + 3
                    if limit is not None and new_end > limit:
                        counts["cds_shift_refused_at_sequence_end"] += 1
                    else:
                        fields[4] = str(new_end)
                        counts["cds_shifted"] += 1
                        counts["cds_shifted_plus"] += 1
                        line = "\t".join(fields) + "\n"
        out.append(line)

    counts["sequence_lengths_known"] = len(seq_len)
    degrade.summary = counts  # type: ignore[attr-defined]
    return out, counts


FIXTURE = """\
##gff-version 3
##sequence-region ctg1 1 1000
ctg1\t.\tgene\t100\t400\t.\t+\t.\tID=gene1
ctg1\t.\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene1
ctg1\t.\texon\t100\t200\t.\t+\t.\tID=e1;Parent=tx1
ctg1\t.\tCDS\t100\t200\t.\t+\t0\tID=c1;Parent=tx1
ctg1\t.\tCDS\t300\t400\t.\t+\t1\tID=c2;Parent=tx1
ctg1\t.\tgene\t500\t800\t.\t-\t.\tID=gene2
ctg1\t.\tmRNA\t500\t800\t.\t-\t.\tID=tx2;Parent=gene2
ctg1\t.\tCDS\t500\t600\t.\t-\t0\tID=c3;Parent=tx2
ctg1\t.\tCDS\t700\t800\t.\t-\t2\tID=c4;Parent=tx2
ctg1\t.\tgene\t1\t50\t.\t-\t.\tID=gene3
ctg1\t.\tmRNA\t1\t50\t.\t-\t.\tID=tx3;Parent=gene3
ctg1\t.\tCDS\t2\t50\t.\t-\t0\tID=c5;Parent=tx3
ctg1\t.\tgene\t900\t1000\t.\t+\t.\tID=gene4
ctg1\t.\tmRNA\t900\t1000\t.\t+\t.\tID=tx4;Parent=gene4
ctg1\t.\tCDS\t900\t1000\t.\t+\t0\tID=c6;Parent=tx4
"""


def self_test() -> int:
    """Check the perturbations on a fixture with both strands and both ends."""
    checks = 0

    def check(got, want, what):
        nonlocal checks
        checks += 1
        if got != want:
            print(f"FAIL {what}: got {got!r}, want {want!r}", file=sys.stderr)
            raise SystemExit(1)

    lines = FIXTURE.splitlines(keepends=True)

    # 1. Rates of zero are the identity.
    out, counts = degrade(lines, 0.0, 0.0, 1)
    check("".join(out), FIXTURE, "zero rates are the identity")
    check(counts["transcripts_seen"], 4, "transcripts seen")
    check(counts["cds_seen"], 6, "CDS seen")
    check(counts["sequence_lengths_known"], 1, "sequence lengths from pragma")

    # 2. Every CDS shifted: strand decides which field moves, and a shift is
    #    refused at either end of the sequence rather than running off it.
    out, counts = degrade(lines, 0.0, 1.0, 1)
    rows = [ln.split("\t") for ln in out if not ln.startswith("#")]
    cds = {r[8].split(";")[0][3:]: (int(r[3]), int(r[4])) for r in rows if r[2] == "CDS"}
    check(cds["c1"], (100, 203), "+ strand internal CDS end moves downstream")
    check(cds["c2"], (300, 403), "+ strand terminal CDS end moves downstream")
    check(cds["c3"], (497, 600), "- strand terminal CDS start moves downstream")
    check(cds["c4"], (697, 800), "- strand internal CDS start moves downstream")
    check(cds["c5"], (2, 50), "- strand CDS at coordinate 2 is not shifted off the start")
    check(cds["c6"], (900, 1000), "+ strand CDS at the sequence end is not shifted off it")
    check(counts["cds_shifted"], 4, "four of six CDS shifted")
    check(counts["cds_shift_refused_at_sequence_end"], 2, "two shifts refused")
    check(counts["cds_shifted_plus"], 2, "plus strand shifts")
    check(counts["cds_shifted_minus"], 2, "minus strand shifts")

    # 3. Dropping a transcript takes its children with it and nothing else.
    out, counts = degrade(lines, 1.0, 0.0, 1)
    kept = [ln for ln in out if not ln.startswith("#")]
    check(counts["transcripts_dropped"], 4, "all transcripts dropped")
    check(counts["child_lines_dropped"], 7, "all transcript children dropped")
    check([r.split("\t")[2] for r in kept], ["gene"] * 4, "only genes survive")

    # 4. The two streams are independent in *both* directions, which is what
    #    makes a drop-only and a shift-only run a decomposition of the combined
    #    one rather than three unrelated files.  Changing the shift rate must
    #    not change which transcripts are dropped, and changing the drop rate
    #    must not change which of the surviving CDS segments are shifted.
    a, ca = degrade(lines, 0.5, 0.0, 7)
    b, cb = degrade(lines, 0.5, 1.0, 7)
    ids = lambda o: [attr(ln.split("\t")[8], "ID") for ln in o
                     if not ln.startswith("#") and ln.split("\t")[2] == "mRNA"]
    check(ids(a), ids(b), "drop draw is independent of the shift rate")
    check(ca["transcripts_dropped"], cb["transcripts_dropped"], "same transcripts dropped")

    cds_coords = lambda o: {ln.split("\t")[8].split(";")[0][3:]:
                            (ln.split("\t")[3], ln.split("\t")[4])
                            for ln in o if not ln.startswith("#")
                            and ln.split("\t")[2] == "CDS"}
    for seed in (3, 7, 11, 23):
        combined, _ = degrade(lines, 0.5, 0.5, seed)
        shift_only, _ = degrade(lines, 0.0, 0.5, seed)
        surviving = cds_coords(combined)
        check(surviving,
              {k: v for k, v in cds_coords(shift_only).items() if k in surviving},
              f"shift draw is independent of the drop rate (seed {seed})")

    # 5. Same seed, same output; different seed, different output.
    c, _ = degrade(lines, 0.5, 0.5, 7)
    d, _ = degrade(lines, 0.5, 0.5, 7)
    check("".join(c), "".join(d), "seed is deterministic")
    differs = any("".join(degrade(lines, 0.5, 0.5, s)[0]) != "".join(c) for s in (8, 9, 10))
    check(differs, True, "a different seed gives a different file")

    # 6. A file with no sequence bounds still degrades, minus the guard.
    nobounds = [ln for ln in lines if not ln.startswith("##sequence-region")]
    out, counts = degrade(nobounds, 0.0, 1.0, 1)
    check(counts["sequence_lengths_known"], 0, "no sequence lengths")
    check(counts["cds_shifted"], 5, "only the coordinate-2 refusal survives")

    print(f"self-test: {checks} checks passed")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", help="reference GFF3 to degrade (.gz ok)")
    ap.add_argument("--out", help="write the degraded GFF3 here")
    ap.add_argument("--summary", help="write the JSON summary here as well as to stderr")
    ap.add_argument("--drop-transcripts", type=float, default=0.10,
                    metavar="P", help="probability of deleting a transcript (default 0.10)")
    ap.add_argument("--shift-cds", type=float, default=0.10, metavar="P",
                    help="probability of moving a CDS segment's downstream boundary "
                         "3 bp further downstream (default 0.10)")
    ap.add_argument("--seed", type=int, default=20260909, help="random seed")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if not args.reference or not args.out:
        ap.error("--reference and --out are required unless --self-test")
    for name, p in (("--drop-transcripts", args.drop_transcripts), ("--shift-cds", args.shift_cds)):
        if not 0.0 <= p <= 1.0:
            ap.error(f"{name} must be between 0 and 1")

    with opener(args.reference) as fh:
        out, counts = degrade(fh, args.drop_transcripts, args.shift_cds, args.seed)
    with open(args.out, "w") as fh:
        fh.writelines(out)

    counts["reference"] = args.reference
    counts["out"] = args.out
    counts["drop_transcripts"] = args.drop_transcripts
    counts["shift_cds"] = args.shift_cds
    counts["seed"] = args.seed
    if not counts["sequence_lengths_known"]:
        print("warning: no ##sequence-region pragmas and no region features; "
              "a shift near a sequence end was not bounds-checked", file=sys.stderr)
    text = json.dumps(counts, indent=2, sort_keys=True)
    print(text, file=sys.stderr)
    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
