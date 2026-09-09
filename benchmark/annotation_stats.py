#!/usr/bin/env python3
"""Summarize a RefSeq/Ensembl-style GFF3 annotation.

Reads one ``*_genomic.gff.gz`` (or plain GFF3) and reports the numbers the
benchmark panel is justified by: intron length quantiles, exons per coding
transcript, CDS length, and coding density.  Standard library only, streams
the file, and never writes anything.

Introns are derived from the exons of ``mRNA`` features and de-duplicated by
(seqid, start, end, strand) so that alternative isoforms sharing an intron
count once.  Quantiles are nearest-rank on the sorted unique-intron list.

Usage:
    python3 annotation_stats.py PATH.gff.gz [--json]
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
from collections import defaultdict

MRNA_TYPES = {"mRNA"}


def _open(path: str) -> io.TextIOBase:
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _attr(field: str, key: str) -> str | None:
    # GFF3 attributes: key=value;key=value.  Values here never contain '='.
    for part in field.split(";"):
        if part.startswith(key) and part[len(key):len(key) + 1] == "=":
            return part[len(key) + 1:]
    return None


def parse(path: str) -> dict:
    mrna_ids: set[str] = set()
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cds_len: dict[str, int] = defaultdict(int)
    mrna_seq: dict[str, tuple[str, str]] = {}
    genome_bp = 0
    genes = 0

    with _open(path) as fh:
        for line in fh:
            if line[0] == "#":
                if line.startswith("##sequence-region"):
                    f = line.split()
                    if len(f) == 4:
                        genome_bp += int(f[3]) - int(f[2]) + 1
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            ftype = f[2]
            if ftype in MRNA_TYPES:
                fid = _attr(f[8], "ID")
                if fid:
                    mrna_ids.add(fid)
                    mrna_seq[fid] = (f[0], f[6])
            elif ftype == "exon":
                parent = _attr(f[8], "Parent")
                if parent in mrna_ids:
                    exons[parent].append((int(f[3]), int(f[4])))
            elif ftype == "CDS":
                parent = _attr(f[8], "Parent")
                if parent in mrna_ids:
                    cds_len[parent] += int(f[4]) - int(f[3]) + 1
            elif ftype == "gene":
                if _attr(f[8], "gene_biotype") in (None, "protein_coding"):
                    genes += 1

    introns: set[tuple[str, int, int, str]] = set()
    exon_counts: list[int] = []
    intronless = 0
    for mid, ex in exons.items():
        ex.sort()
        exon_counts.append(len(ex))
        if len(ex) == 1:
            intronless += 1
        seqid, strand = mrna_seq.get(mid, ("?", "?"))
        for i in range(len(ex) - 1):
            start, end = ex[i][1] + 1, ex[i + 1][0] - 1
            if end >= start:
                introns.add((seqid, start, end, strand))

    lens = sorted(e - s + 1 for _, s, e, _ in introns)
    cds = sorted(cds_len.values())
    exon_counts.sort()

    def q(xs: list[int], p: float) -> int | None:
        if not xs:
            return None
        return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))]

    return {
        "protein_coding_genes": genes,
        "mrna_transcripts": len(exons),
        "intronless_transcripts": intronless,
        "unique_introns": len(lens),
        "intron_len_p10": q(lens, 0.10),
        "intron_len_median": q(lens, 0.50),
        "intron_len_p90": q(lens, 0.90),
        "intron_len_p99": q(lens, 0.99),
        "intron_len_max": lens[-1] if lens else None,
        "intron_len_mean": round(sum(lens) / len(lens), 1) if lens else None,
        "exons_per_transcript_median": q(exon_counts, 0.50),
        "exons_per_transcript_mean": round(sum(exon_counts) / len(exon_counts), 2) if exon_counts else None,
        "cds_len_median": q(cds, 0.50),
        "cds_total_bp": sum(cds),
        "sequence_region_bp": genome_bp or None,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gff")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of key\\tvalue lines")
    args = ap.parse_args()
    stats = parse(args.gff)
    if args.json:
        print(json.dumps(stats))
    else:
        for k, v in stats.items():
            print(f"{k}\t{v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
