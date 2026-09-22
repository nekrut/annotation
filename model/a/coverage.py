"""Turnkey train-set coverage for candidate A (a-pilot.md section 7).

The section-3.6 loader (`model.a.dataset`) reports, per species, how much of the
admitted set the current complete-target/clean-window scope actually trains on
(`coverage_report`/`format_coverage`). Running it needs the *checksummed* source
GFF3/FASTA that the committed manifests were built from, which are too large to
commit (charter: never commit files over 5 MB). Both incremental reviews
(engels-0065, stalin-0068) asked for the panel-wide numbers before any
train-panel conclusion is drawn from A, so producing them must be one command,
not a manual checkout.

This module is that command. It

- constructs the deterministic NCBI FTP URL of each pinned source from the
  filename recorded in the species `*.summary.json` (`ncbi_url`), so no URL is
  guessed or scraped;
- downloads the exact `_genomic.gff.gz`/`_genomic.fna.gz` and verifies each
  against the summary's MD5 (`fetch`) -- the same digest `dataset.verify_source`
  re-checks before yielding a window, so a corrupted or wrong-assembly download
  can never reach the loader;
- runs the existing `coverage_report` over the checked-out sources and renders
  the section-7 TSV (`format_coverage`).

The report itself re-derives nothing: it delegates admission to
`model.labels.admission.audit_species` and the checksum gate to
`dataset.verify_source`, exactly as training would, so the printed numbers are
the fitting population under the current scope. `coverage_report` runs species
sequentially and retains only counts, so only one genome is resident at a time,
not the panel sum; all ten train species were reported on a laptop this way. The
manifests' `peak_rss_mb` is the *standalone* `audit_species` high-water mark
(0.08-3.6 GB across the ten, 1.9-3.6 GB for the four large ones; MiB despite the
name) -- the full coverage-loader peak, which keeps the audit result live while
re-parsing the GFF and building windows, is higher by an unmeasured margin and
is not reported here. A whole-panel `--fetch` additionally needs scratch disk
for the ten compressed genomes at once, since the fetch retains each download.
Standard library only (urllib for the fetch); Python 3.11.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
import urllib.request
from typing import List, Optional, Tuple

NCBI_ALL = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
Source = Tuple[str, str, str, str]  # (species, summary, gff, fasta)


def _accession_and_dir(filename: str) -> Tuple[str, str]:
    """``(accession, assembly_dir)`` for an NCBI ``<acc>_<asm>_genomic.<ext>.gz``
    filename. The accession is the first two ``_``-delimited tokens (for example
    ``GCF_000146045.2``); the assembly directory is the filename up to
    ``_genomic.`` (``GCF_000146045.2_R64``), so an assembly name that itself
    contains underscores (``Release_6_plus_ISO1_MT``) is kept whole."""
    marker = "_genomic."
    if marker not in filename:
        raise ValueError(f"not an NCBI genomic filename: {filename!r}")
    asm_dir = filename[: filename.index(marker)]
    parts = filename.split("_")
    if len(parts) < 3 or "." not in parts[1] or parts[0] not in ("GCF", "GCA"):
        raise ValueError(f"cannot parse accession from {filename!r}")
    return parts[0] + "_" + parts[1], asm_dir


def ncbi_url(filename: str) -> str:
    """The deterministic NCBI ``genomes/all`` FTP URL for a pinned genomic file.

    The nine accession digits split into three-character path components, for
    example ``GCF_000146045.2_R64_genomic.fna.gz`` ->
    ``.../GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz``.
    """
    accession, asm_dir = _accession_and_dir(filename)
    prefix, digits = accession.split("_", 1)
    number = digits.split(".")[0]
    if len(number) != 9 or not number.isdigit():
        raise ValueError(f"unexpected accession number in {accession!r}")
    triplets = "/".join(number[i:i + 3] for i in range(0, 9, 3))
    return f"{NCBI_ALL}/{prefix}/{triplets}/{asm_dir}/{filename}"


def species_sources(manifests_dir: str, sources_dir: str) -> List[Source]:
    """``(species, summary, gff, fasta)`` for every ``*.summary.json`` under
    ``manifests_dir``, with the GFF/FASTA resolved to their expected paths under
    ``sources_dir/<species>/`` using the filenames the summary pins. The files
    need not exist yet -- :func:`fetch` puts them there and
    :func:`dataset.coverage_report` checks their digests."""
    out: List[Source] = []
    for name in sorted(os.listdir(manifests_dir)):
        if not name.endswith(".summary.json"):
            continue
        summary = os.path.join(manifests_dir, name)
        with open(summary, "r", encoding="utf-8") as fh:
            meta = json.load(fh)
        species = meta["species"]
        d = os.path.join(sources_dir, species)
        out.append((species, summary,
                    os.path.join(d, meta["gff"]), os.path.join(d, meta["fasta"])))
    return out


def _md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: str) -> None:
    """Download ``url`` to ``dest`` atomically (temp file in the same directory,
    renamed on success), so an interrupted fetch never leaves a half file that
    would fail the digest check on a rerun."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(dest), suffix=".part")
    os.close(fd)
    try:
        with urllib.request.urlopen(url) as resp, open(tmp, "wb") as out:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def fetch(manifests_dir: str, sources_dir: str,
          species: Optional[str] = None, log=print) -> List[Source]:
    """Download and MD5-verify every pinned source (or one ``species``) into
    ``sources_dir/<species>/``. A file already present and matching its manifest
    digest is left alone, so the fetch is resumable; a present file with the
    wrong digest is re-downloaded. Returns the resolved sources."""
    sources = species_sources(manifests_dir, sources_dir)
    if species is not None:
        sources = [s for s in sources if s[0] == species]
        if not sources:
            raise ValueError(f"no summary for species {species!r}")
    for sp, summary, gff, fasta in sources:
        with open(summary, "r", encoding="utf-8") as fh:
            meta = json.load(fh)
        for path, filename, want in (
            (gff, meta["gff"], meta["gff_md5"]),
            (fasta, meta["fasta"], meta["fasta_md5"]),
        ):
            if os.path.exists(path) and _md5(path) == want:
                log(f"{sp}: {filename} present, md5 ok")
                continue
            url = ncbi_url(filename)
            log(f"{sp}: fetching {url}")
            _download(url, path)
            got = _md5(path)
            if got != want:
                raise RuntimeError(
                    f"{sp}: {filename} md5 {got} != manifest {want}")
            log(f"{sp}: {filename} md5 ok")
    return sources


def _report(sources: List[Source], complete_only: bool,
            max_window: Optional[int]) -> str:
    # Imported here so `fetch`/URL construction stay usable on a host that
    # cannot import the loader's dependencies.
    from model.a.dataset import coverage_report, format_coverage

    return format_coverage(coverage_report(
        sources, complete_only=complete_only, max_window=max_window))


def main(argv: Optional[List[str]] = None) -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    default_manifests = os.path.normpath(
        os.path.join(here, "..", "labels", "manifests"))
    ap = argparse.ArgumentParser(
        description="Fetch pinned train sources and report loader coverage "
                    "(a-pilot.md section 7).")
    ap.add_argument("--manifests", default=default_manifests,
                    help="directory of <species>.summary.json (default: "
                         "model/labels/manifests)")
    ap.add_argument("--sources", required=True,
                    help="scratch directory for the checked-out GFF3/FASTA")
    ap.add_argument("--species", help="restrict to one species")
    ap.add_argument("--fetch", action="store_true",
                    help="download and MD5-verify the sources before reporting")
    ap.add_argument("--all", dest="complete_only", action="store_false",
                    help="include edge-partial chains (complete_only=False)")
    ap.add_argument("--max-window", type=int, default=None,
                    help="skip windows longer than this many bases")
    ap.add_argument("--out", help="write the TSV here (default: stdout)")
    ap.add_argument("--self-test", action="store_true",
                    help="offline checks of URL construction, then exit")
    args = ap.parse_args(argv)

    if args.self_test:
        _self_test()
        print("coverage self-test ok")
        return 0

    # Resolve and validate the selection up front, in both modes, so an empty
    # manifest directory or an unknown --species is a clean argument error
    # rather than a valid-looking zero-row table. A known species with zero
    # admitted labels still resolves here and reports as its own row.
    sources = species_sources(args.manifests, args.sources)
    if not sources:
        ap.error(f"no *.summary.json manifests found in {args.manifests}")
    if args.species is not None:
        selected = [s for s in sources if s[0] == args.species]
        if not selected:
            ap.error(f"unknown species {args.species!r}; manifests offer: "
                     + ", ".join(s[0] for s in sources))
        sources = selected

    if args.fetch:
        # Route fetch progress to stderr so a redirected/piped stdout is a clean
        # TSV even on a fully cached run.
        sources = fetch(args.manifests, args.sources, species=args.species,
                        log=lambda m: print(m, file=sys.stderr))
    else:
        missing = [p for _, _, g, f in sources
                   for p in (g, f) if not os.path.exists(p)]
        if missing:
            ap.error("missing source files (run with --fetch):\n  "
                     + "\n  ".join(missing))

    table = _report(sources, args.complete_only, args.max_window)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(table)
    else:
        sys.stdout.write(table)
    return 0


def _self_test() -> None:
    cases = {
        "GCF_000146045.2_R64_genomic.fna.gz":
            NCBI_ALL + "/GCF/000/146/045/GCF_000146045.2_R64/"
            "GCF_000146045.2_R64_genomic.fna.gz",
        "GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz":
            NCBI_ALL + "/GCF/902/167/145/"
            "GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/"
            "GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz",
        "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz":
            NCBI_ALL + "/GCF/000/001/215/"
            "GCF_000001215.4_Release_6_plus_ISO1_MT/"
            "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
    }
    for filename, want in cases.items():
        got = ncbi_url(filename)
        assert got == want, (filename, got, want)
    for bad in ("random.txt", "GCF_000146045.2_R64.fna.gz",
                "notacc_000146045.2_R64_genomic.fna.gz"):
        try:
            ncbi_url(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted malformed filename {bad!r}")


if __name__ == "__main__":
    raise SystemExit(main())
