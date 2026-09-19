"""Tests for the turnkey coverage command (`model.a.coverage`).

Standard library only. URL construction and source resolution are checked
offline; the fetch path is exercised with a stubbed downloader (no network);
the CLI is run end to end over a tiny synthetic species (one admitted two-exon
gene) written to temp files, so it produces the section-7 TSV on a host without
torch and without touching NCBI.
"""
import gzip
import hashlib
import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stdout

from model.a import coverage as cov

# One admitted complete gene on chr1 (same construction as test_a_dataset):
# exon1 11..25 = ATG + 12 A, 25-base intron, exon2 51..62 = 9 A + TAA.
CONTIG_LEN = 10000
EXON1 = "ATGAAAAAAAAAAAA"
INTRON = "GT" + "C" * 21 + "AG"
EXON2 = "AAAAAAAAATAA"
GFF = "\n".join([
    "##gff-version 3",
    "##sequence-region chr1 1 %d" % CONTIG_LEN,
    "chr1\ttest\tgene\t11\t62\t.\t+\t.\tID=gene1;gene_biotype=protein_coding",
    "chr1\ttest\tmRNA\t11\t62\t.\t+\t.\tID=tx1;Parent=gene1",
    "chr1\ttest\tCDS\t11\t25\t.\t+\t0\tID=cds1;Parent=tx1",
    "chr1\ttest\tCDS\t51\t62\t.\t+\t0\tID=cds1;Parent=tx1",
    "",
])
GFF_NAME = "GCF_000146045.2_R64_genomic.gff.gz"
FASTA_NAME = "GCF_000146045.2_R64_genomic.fna.gz"


def _contig():
    head = "C" * 10 + EXON1 + INTRON + EXON2
    return head + "A" * (CONTIG_LEN - len(head))


def _fasta():
    seq = _contig()
    lines = [">chr1 test contig"] + [seq[i:i + 70] for i in range(0, len(seq), 70)]
    return "\n".join(lines) + "\n"


# The loader gunzips sources by their ``.gz`` suffix, so the synthetic sources
# are real gzip bytes (deterministic ``mtime=0``); the summary pins the md5 of
# those exact bytes, matching how the loader's checksum gate hashes the file.
GFF_GZ = gzip.compress(GFF.encode(), mtime=0)
FASTA_GZ = gzip.compress(_fasta().encode(), mtime=0)
GFF_MD5 = hashlib.md5(GFF_GZ).hexdigest()
FASTA_MD5 = hashlib.md5(FASTA_GZ).hexdigest()


class UrlTest(unittest.TestCase):
    def test_self_test_passes(self):
        cov._self_test()  # raises on any mismatch

    def test_known_urls(self):
        base = cov.NCBI_ALL
        self.assertEqual(
            cov.ncbi_url("GCF_000146045.2_R64_genomic.fna.gz"),
            base + "/GCF/000/146/045/GCF_000146045.2_R64/"
                   "GCF_000146045.2_R64_genomic.fna.gz")
        # an assembly name with underscores stays whole in the directory
        self.assertEqual(
            cov.ncbi_url("GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz"),
            base + "/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/"
                   "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz")

    def test_accession_and_dir(self):
        acc, d = cov._accession_and_dir(
            "GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz")
        self.assertEqual(acc, "GCA_902167145.1")
        self.assertEqual(d, "GCA_902167145.1_Zm-B73-REFERENCE-NAM-5.0")

    def test_malformed_filenames_rejected(self):
        for bad in ("random.txt", "GCF_000146045.2_R64.fna.gz",
                    "notacc_000146045.2_R64_genomic.fna.gz"):
            with self.assertRaises(ValueError):
                cov.ncbi_url(bad)

    def test_real_manifests_construct_urls(self):
        # every committed summary yields a well-formed URL, if the tree is present
        here = os.path.dirname(os.path.abspath(__file__))
        manifests = os.path.normpath(
            os.path.join(here, "..", "model", "labels", "manifests"))
        if not os.path.isdir(manifests):
            self.skipTest("manifests not checked out")
        for _sp, summary, gff, fasta in cov.species_sources(manifests, "/scratch"):
            for path in (gff, fasta):
                url = cov.ncbi_url(os.path.basename(path))
                self.assertTrue(url.startswith(cov.NCBI_ALL + "/GC"))
                self.assertTrue(url.endswith(os.path.basename(path)))


class SourcesAndReportTest(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.manifests = os.path.join(self.root, "manifests")
        self.sources = os.path.join(self.root, "sources")
        os.makedirs(self.manifests)
        with open(os.path.join(self.manifests, "Yeast.summary.json"), "w") as fh:
            json.dump({"species": "Yeast", "gff": GFF_NAME, "fasta": FASTA_NAME,
                       "gff_md5": GFF_MD5, "fasta_md5": FASTA_MD5}, fh)

    def _write_sources(self):
        d = os.path.join(self.sources, "Yeast")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, GFF_NAME), "wb") as fh:
            fh.write(GFF_GZ)
        with open(os.path.join(d, FASTA_NAME), "wb") as fh:
            fh.write(FASTA_GZ)

    def test_species_sources_resolves_paths(self):
        srcs = cov.species_sources(self.manifests, self.sources)
        self.assertEqual(len(srcs), 1)
        sp, summary, gff, fasta = srcs[0]
        self.assertEqual(sp, "Yeast")
        self.assertEqual(gff, os.path.join(self.sources, "Yeast", GFF_NAME))
        self.assertEqual(fasta, os.path.join(self.sources, "Yeast", FASTA_NAME))

    def test_fetch_verifies_and_is_resumable(self):
        # stub the network: "download" writes the correct synthetic bytes
        payload = {GFF_NAME: GFF_GZ, FASTA_NAME: FASTA_GZ}
        calls = []

        def fake_download(url, dest):
            calls.append(url)
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "wb") as fh:
                fh.write(payload[os.path.basename(dest)])

        orig = cov._download
        cov._download = fake_download
        try:
            cov.fetch(self.manifests, self.sources, log=lambda *a: None)
            self.assertEqual(len(calls), 2)          # gff + fasta downloaded
            # both files present and md5-verified
            gff = os.path.join(self.sources, "Yeast", GFF_NAME)
            self.assertEqual(cov._md5(gff), GFF_MD5)
            # a second fetch downloads nothing (present and verified)
            cov.fetch(self.manifests, self.sources, log=lambda *a: None)
            self.assertEqual(len(calls), 2)
        finally:
            cov._download = orig

    def test_fetch_raises_on_digest_mismatch(self):
        def bad_download(url, dest):
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "w") as fh:
                fh.write("WRONG")

        orig = cov._download
        cov._download = bad_download
        try:
            with self.assertRaises(RuntimeError):
                cov.fetch(self.manifests, self.sources, log=lambda *a: None)
        finally:
            cov._download = orig

    def test_main_reports_coverage_tsv(self):
        self._write_sources()
        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = cov.main(["--manifests", self.manifests,
                           "--sources", self.sources])
        self.assertEqual(rc, 0)
        lines = buf.getvalue().splitlines()
        self.assertEqual(lines[0].split("\t")[0], "species")
        # one admitted, one yielded, 100%
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        self.assertEqual(row["species"], "Yeast")
        self.assertEqual(row["admitted"], "1")
        self.assertEqual(row["yielded"], "1")
        self.assertEqual(row["yielded_pct"], "100.0")
        self.assertEqual(lines[-1].split("\t")[0], "TOTAL")

    def test_main_errors_when_sources_missing(self):
        with self.assertRaises(SystemExit):
            cov.main(["--manifests", self.manifests, "--sources", self.sources])


if __name__ == "__main__":
    unittest.main()
