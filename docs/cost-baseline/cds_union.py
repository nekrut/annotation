#!/usr/bin/env python3
"""Genomic coding fraction: the union of CDS intervals over the genome.

``benchmark/panel.tsv``'s ``cds_fraction_pct`` sums CDS length over every
mRNA, so a gene with seven isoforms counts its shared exons seven times.
This script reports the fraction of genome bases covered by at least one
CDS feature (either strand), which is the number a candidate-region stage
must keep.  It also recomputes the per-isoform sum so the two can be read
side by side.

Numerator and denominator must cover the same sequences.  ``genome_bp`` in
the panel is the primary assembly without organelles, but a RefSeq GFF
annotates every sequence in the assembly: for GRCh38.p14 that adds 515
alt-scaffold, fix-patch and novel-patch sequences (199 Mb, largely
duplicates of primary-assembly genes such as the MHC) plus chrM; mouse,
maize and yeast have no alt or patch units, only organelles.  By default
the script therefore keeps only sequences whose role in the NCBI assembly
report is ``assembled-molecule``, ``unlocalized-scaffold`` or
``unplaced-scaffold`` and whose assembly unit is not ``non-nuclear``, and
checks that their lengths sum to the panel's ``genome_bp``.  ``--all``
restores the unfiltered union (every seqid in the GFF) for comparison.

Usage:
    python3 docs/cost-baseline/cds_union.py [--all] Homo_sapiens [Mus_musculus ...]

Reads ``benchmark/panel.tsv`` for the accession, FTP directory, genome size
and GFF md5; downloads ``*_genomic.gff.gz`` and ``*_assembly_report.txt``
into ``$CDS_UNION_CACHE`` (default ``/tmp/cds_union``) if absent; verifies
the GFF md5; prints one TSV row per species.  Standard library only.
"""
import csv
import gzip
import hashlib
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PANEL = os.path.join(ROOT, "benchmark", "panel.tsv")
FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
UA = "relay-benchmark/0.1 (https://github.com/nekrut/annotation)"
CACHE = os.environ.get("CDS_UNION_CACHE", "/tmp/cds_union")
PRIMARY_ROLES = {"assembled-molecule", "unlocalized-scaffold", "unplaced-scaffold"}


def ftp_dir_url(accession: str, ftp_dir: str) -> str:
    prefix, digits = accession.split("_", 1)
    num = digits.split(".", 1)[0]
    parts = [num[i:i + 3] for i in range(0, 9, 3)]
    return f"{FTP_ROOT}/{prefix}/{parts[0]}/{parts[1]}/{parts[2]}/{ftp_dir}/{ftp_dir}"


def fetch(url: str, dest: str, md5: str = "") -> None:
    if not os.path.exists(dest):
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r, open(dest + ".part", "wb") as w:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                w.write(chunk)
        os.rename(dest + ".part", dest)
    if not md5:
        return
    h = hashlib.md5()
    with open(dest, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != md5:
        sys.exit(f"md5 mismatch for {dest}: {h.hexdigest()} != {md5}")


def primary_seqids(report_path: str) -> dict:
    """RefSeq accession -> length for primary-assembly, nuclear sequences."""
    keep = {}
    with open(report_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fld = line.rstrip("\n").split("\t")
            role, unit, refseq, length = fld[1], fld[7], fld[6], fld[8]
            if role in PRIMARY_ROLES and unit != "non-nuclear" and refseq != "na":
                keep[refseq] = int(length)
    return keep


def union_length(intervals: list) -> int:
    intervals.sort()
    total, cur_s, cur_e = 0, None, None
    for s, e in intervals:
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        elif e > cur_e:
            cur_e = e
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def measure(path: str, keep: dict = None) -> tuple:
    """Return (union, per-isoform sum, CDS count) over kept seqids, and the
    union over excluded seqids.  ``keep=None`` keeps every seqid."""
    by_seq: dict = {}
    per_isoform = 0
    n_cds = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fld = line.rstrip("\n").split("\t")
            if len(fld) != 9 or fld[2] != "CDS":
                continue
            s, e = int(fld[3]), int(fld[4])
            by_seq.setdefault(fld[0], []).append((s, e))
    excluded = 0
    for seqid, iv in by_seq.items():
        if keep is not None and seqid not in keep:
            excluded += union_length(iv)
            continue
        per_isoform += sum(e - s + 1 for s, e in iv)
        n_cds += len(iv)
    union = sum(union_length(iv) for seqid, iv in by_seq.items() if keep is None or seqid in keep)
    return union, per_isoform, n_cds, excluded


def main(argv: list) -> None:
    use_all = "--all" in argv
    species = [a for a in argv if not a.startswith("--")]
    rows = {r["species"]: r for r in csv.DictReader(open(PANEL), delimiter="\t")}
    print("species\tgenome_bp\tsequences\tcds_union_bp\tcds_union_pct\tcds_isoform_sum_bp\tcds_isoform_sum_pct\tcds_features\tisoform_inflation\texcluded_cds_union_bp")
    for sp in species:
        r = rows[sp]
        base = ftp_dir_url(r["accession"], r["ftp_dir"])
        gff = os.path.join(CACHE, os.path.basename(base) + "_genomic.gff.gz")
        fetch(base + "_genomic.gff.gz", gff, r["md5_gff"])
        g = int(r["genome_bp"])
        keep = None
        if not use_all:
            report = os.path.join(CACHE, os.path.basename(base) + "_assembly_report.txt")
            fetch(base + "_assembly_report.txt", report)
            keep = primary_seqids(report)
            kept_bp = sum(keep.values())
            if kept_bp != g:
                print(f"warning: {sp}: primary nuclear sequences sum to {kept_bp} bp, panel genome_bp is {g}", file=sys.stderr)
        union, iso, n, excluded = measure(gff, keep)
        print(f"{sp}\t{g}\t{'all' if use_all else 'primary'}\t{union}\t{100*union/g:.2f}\t{iso}\t{100*iso/g:.1f}\t{n}\t{iso/union:.2f}\t{excluded}")
        sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
