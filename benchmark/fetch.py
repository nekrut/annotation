#!/usr/bin/env python3
"""Fetch benchmark panel genomes and annotations from the NCBI FTP mirror.

Reads ``benchmark/panel.tsv``, resolves each assembly accession to its
directory under ``https://ftp.ncbi.nlm.nih.gov/genomes/all/``, downloads the
requested files into a destination directory *outside* the repository, and
verifies them against the MD5 checksums NCBI publishes alongside them.

Nothing is written into the repository: the panel manifest records what to
fetch and the checksum to expect, the data itself stays on the machine that
runs the benchmark.  Standard library only.

Examples:
    python3 benchmark/fetch.py --species Saccharomyces_cerevisiae --dest /tmp/panel
    python3 benchmark/fetch.py --all --what gff --dest /data/panel
    python3 benchmark/fetch.py --species Homo_sapiens --dry-run

Exit status is non-zero if any requested file failed to download or failed
its checksum.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import sys
import urllib.error
import urllib.request

FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
UA = "relay-benchmark/0.1 (https://github.com/nekrut/annotation)"
SUFFIX = {
    "gff": "_genomic.gff.gz",
    "fasta": "_genomic.fna.gz",
    "protein": "_protein.faa.gz",
    "cds": "_cds_from_genomic.fna.gz",
}
PANEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "panel.tsv")


def asm_parent(accession: str) -> str:
    """https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405"""
    prefix, digits = accession.split("_", 1)
    num = digits.split(".", 1)[0]
    parts = [num[i:i + 3] for i in range(0, 9, 3)]
    return f"{FTP_ROOT}/{prefix}/{parts[0]}/{parts[1]}/{parts[2]}"


def asm_dir(accession: str, assembly_name: str, ftp_dir: str = "") -> tuple[str, str]:
    """Resolve the assembly directory and the file-name stem inside it.

    The directory is named after the *submitter's* assembly name, which is not
    always the ``assembly_name`` the Datasets API reports: GCF_034140825.1 is
    ``AGIS1.0`` in the API and ``ASM3414082v1`` on the FTP site.  ``panel.tsv``
    therefore records the resolved directory name in ``ftp_dir``; when it is
    absent or stale we fall back to listing the parent directory.
    """
    parent = asm_parent(accession)
    stem = ftp_dir.strip() or f"{accession}_{assembly_name.replace(' ', '_')}"
    if ftp_dir.strip():
        return f"{parent}/{stem}", stem
    try:
        listing = get(f"{parent}/").decode("utf-8", "replace")
    except (urllib.error.URLError, OSError):
        return f"{parent}/{stem}", stem
    found = re.findall(re.escape(accession) + r"_[^\"/<>\s]+", listing)
    if found:
        stem = sorted(found, key=len)[0]
    return f"{parent}/{stem}", stem


def read_panel(path: str = PANEL) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if not r["species"].startswith("#")]
    return rows


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as resp:
        return resp.read()


def download(url: str, dest: str) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    tmp = dest + ".part"
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    os.replace(tmp, dest)


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ncbi_md5sums(base: str) -> dict[str, str]:
    """Parse the md5checksums.txt NCBI ships in every assembly directory."""
    out = {}
    for line in get(f"{base}/md5checksums.txt").decode().splitlines():
        digest, _, name = line.partition("  ")
        out[os.path.basename(name.strip("./ "))] = digest
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=PANEL)
    ap.add_argument("--species", action="append", default=[], help="species column value; repeatable")
    ap.add_argument("--clade", action="append", default=[], help="clade column value; repeatable")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--what", default="gff", help="comma-separated: " + ",".join(SUFFIX))
    ap.add_argument("--dest", default="./panel-data", help="download directory (keep it out of the repo)")
    ap.add_argument("--dry-run", action="store_true", help="print URLs, download nothing")
    args = ap.parse_args()

    rows = read_panel(args.panel)
    if not args.all:
        rows = [r for r in rows if r["species"] in args.species or r["clade"] in args.clade]
    if not rows:
        ap.error("no species selected: use --all, --species NAME, or --clade NAME")
    kinds = [k.strip() for k in args.what.split(",")]
    bad = [k for k in kinds if k not in SUFFIX]
    if bad:
        ap.error(f"unknown --what value(s): {', '.join(bad)}")

    failures = 0
    for row in rows:
        base, stem = asm_dir(row["accession"], row["assembly_name"], row.get("ftp_dir", ""))
        sums = None
        for kind in kinds:
            name = stem + SUFFIX[kind]
            url = f"{base}/{name}"
            if args.dry_run:
                print(url)
                continue
            os.makedirs(os.path.join(args.dest, row["species"]), exist_ok=True)
            out = os.path.join(args.dest, row["species"], name)
            try:
                if not os.path.exists(out):
                    download(url, out)
                if sums is None:
                    sums = ncbi_md5sums(base)
                want = sums.get(name)
                got = md5(out)
                if want and want != got:
                    print(f"CHECKSUM MISMATCH {out}: want {want}, got {got}", file=sys.stderr)
                    failures += 1
                    continue
                expected = row.get(f"md5_{kind}", "").strip()
                if expected and expected not in ("", "-", "?") and expected != got:
                    print(f"PANEL CHECKSUM MISMATCH {out}: panel.tsv says {expected}, got {got}",
                          file=sys.stderr)
                    failures += 1
                    continue
                print(f"ok\t{row['species']}\t{kind}\t{got}\t{out}")
            except (urllib.error.URLError, OSError) as exc:
                print(f"FAILED {url}: {exc}", file=sys.stderr)
                failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
