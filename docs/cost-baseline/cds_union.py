#!/usr/bin/env python3
"""Genomic coding fraction: the union of CDS intervals over the genome.

``benchmark/panel.tsv``'s ``cds_fraction_pct`` sums CDS length over every
mRNA, so a gene with seven isoforms counts its shared exons seven times.
This script reports the fraction of genome bases covered by at least one
CDS feature (either strand), which is the number a candidate-region stage
must keep.  It also recomputes the per-isoform sum so the two can be read
side by side.

Usage:
    python3 docs/cost-baseline/cds_union.py Homo_sapiens [Mus_musculus ...]

Reads ``benchmark/panel.tsv`` for the accession, FTP directory, genome size
and GFF md5; downloads ``*_genomic.gff.gz`` into ``$CDS_UNION_CACHE`` (default
``/tmp/cds_union``) if absent; verifies the md5; prints one TSV row per
species.  Standard library only.
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


def gff_url(accession: str, ftp_dir: str) -> str:
    prefix, digits = accession.split("_", 1)
    num = digits.split(".", 1)[0]
    parts = [num[i:i + 3] for i in range(0, 9, 3)]
    return f"{FTP_ROOT}/{prefix}/{parts[0]}/{parts[1]}/{parts[2]}/{ftp_dir}/{ftp_dir}_genomic.gff.gz"


def fetch(url: str, dest: str, md5: str) -> None:
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
    h = hashlib.md5()
    with open(dest, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    if md5 and h.hexdigest() != md5:
        sys.exit(f"md5 mismatch for {dest}: {h.hexdigest()} != {md5}")


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


def measure(path: str) -> tuple:
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
            per_isoform += e - s + 1
            n_cds += 1
            by_seq.setdefault(fld[0], []).append((s, e))
    union = sum(union_length(v) for v in by_seq.values())
    return union, per_isoform, n_cds


def main(argv: list) -> None:
    rows = {r["species"]: r for r in csv.DictReader(open(PANEL), delimiter="\t")}
    print("species\tgenome_bp\tcds_union_bp\tcds_union_pct\tcds_isoform_sum_bp\tcds_isoform_sum_pct\tcds_features\tisoform_inflation")
    for sp in argv:
        r = rows[sp]
        url = gff_url(r["accession"], r["ftp_dir"])
        dest = os.path.join(CACHE, os.path.basename(url))
        fetch(url, dest, r["md5_gff"])
        union, iso, n = measure(dest)
        g = int(r["genome_bp"])
        print(f"{sp}\t{g}\t{union}\t{100*union/g:.2f}\t{iso}\t{100*iso/g:.1f}\t{n}\t{iso/union:.2f}")
        sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
