"""Tests for the section-3.6 structured training-window loader (`model.a.dataset`).

Standard library only: a tiny synthetic species (one admitted two-exon gene on
a >=10 kb contig) is written to temp GFF3/FASTA/summary files, so the loader is
exercised end to end -- checksum gate, admission delegation, oriented window,
and the round-trip through the chain-loss oracle -- on a host without torch.
"""
import hashlib
import json
import math
import os
import tempfile
import unittest

from model.a.dataset import (
    LoaderStats,
    SourceMismatch,
    WindowExample,
    iter_windows,
    verify_source,
)
from model.a.loss import chain_nll
from model.grammar import ReferenceDecoder
from model.labels.admission import revcomp

# One admitted complete gene on chr1: exon1 = 11..25 (ATG + 12 A, 5 codons),
# a 25-base intron 26..50, exon2 = 51..62 (9 A + TAA stop, 4 codons). Spliced
# CDS ATG(AAA)x7 TAA: in frame, initiator, terminal stop, no internal stop.
CONTIG_LEN = 10000
EXON1 = "ATGAAAAAAAAAAAA"          # 15 bases, positions 11..25
INTRON = "GT" + "C" * 21 + "AG"     # 25 bases, positions 26..50
EXON2 = "AAAAAAAAATAA"             # 12 bases, positions 51..62


def _contig():
    head = "C" * 10 + EXON1 + INTRON + EXON2      # bases 1..62
    return head + "A" * (CONTIG_LEN - len(head))  # pad to 10 kb


GFF = "\n".join([
    "##gff-version 3",
    "##sequence-region chr1 1 %d" % CONTIG_LEN,
    "chr1\ttest\tgene\t11\t62\t.\t+\t.\tID=gene1;gene_biotype=protein_coding",
    "chr1\ttest\tmRNA\t11\t62\t.\t+\t.\tID=tx1;Parent=gene1",
    "chr1\ttest\tCDS\t11\t25\t.\t+\t0\tID=cds1;Parent=tx1",
    "chr1\ttest\tCDS\t51\t62\t.\t+\t0\tID=cds1;Parent=tx1",
    "",
])


def _fasta():
    seq = _contig()
    lines = [">chr1 test contig"]
    lines += [seq[i:i + 70] for i in range(0, len(seq), 70)]
    return "\n".join(lines) + "\n"


def _md5(text):
    return hashlib.md5(text.encode()).hexdigest()


class DatasetTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.gff = os.path.join(self.dir, "x.gff3")
        self.fasta = os.path.join(self.dir, "x.fna")
        self.summary = os.path.join(self.dir, "x.summary.json")
        with open(self.gff, "w") as fh:
            fh.write(GFF)
        with open(self.fasta, "w") as fh:
            fh.write(_fasta())
        with open(self.summary, "w") as fh:
            json.dump({"gff_md5": _md5(GFF), "fasta_md5": _md5(_fasta())}, fh)

    # -- checksum gate -------------------------------------------------------
    def test_verify_source_ok(self):
        got = verify_source(self.summary, self.gff, self.fasta)
        self.assertEqual(got["gff_md5"], _md5(GFF))
        self.assertEqual(got["fasta_md5"], _md5(_fasta()))

    def test_verify_source_detects_tampered_fasta(self):
        with open(self.fasta, "a") as fh:
            fh.write("A\n")            # one base changed after the manifest
        with self.assertRaises(SourceMismatch):
            verify_source(self.summary, self.gff, self.fasta)

    def test_verify_source_fasta_presence_must_match(self):
        # summary recorded a fasta digest; verifying without a fasta is a mismatch
        with self.assertRaises(SourceMismatch):
            verify_source(self.summary, self.gff, None)

    # -- windowing -----------------------------------------------------------
    def test_iter_windows_yields_admitted_complete_chain(self):
        stats = LoaderStats()
        windows = list(iter_windows(self.summary, self.gff, self.fasta, stats=stats))
        self.assertEqual(len(windows), 1)
        ex = windows[0]
        self.assertIsInstance(ex, WindowExample)
        self.assertEqual(ex.key, ("chr1", "+", "tx1"))
        # FLANK=10 on each side: window is bases 1..72, cds/introns rebased to ws=1
        self.assertEqual(ex.n, 72)
        self.assertEqual(ex.cds_ranges, [(10, 25), (50, 62)])
        self.assertEqual(ex.intron_ranges, [(25, 50)])
        self.assertFalse(ex.partial_5 or ex.partial_3)
        # the window really carries the initiator at the first CDS base
        self.assertEqual(ex.window[10:13], "ATG")
        self.assertEqual(stats.admitted, 1)
        self.assertEqual(stats.yielded, 1)
        self.assertEqual(stats.skipped_partial, 0)

    def test_window_round_trips_through_chain_loss(self):
        ex = next(iter_windows(self.summary, self.gff, self.fasta))
        loss, log_z, log_z_num = chain_nll(
            ReferenceDecoder(), ex.window, None, ex.cds_ranges, ex.intron_ranges)
        # a legal admitted chain has a finite numerator and nonnegative loss
        self.assertTrue(math.isfinite(log_z) and math.isfinite(log_z_num))
        self.assertGreaterEqual(loss, -1e-9)
        # the support-masked numerator Viterbi-decodes back to the gold chain
        num = ex.support(None)
        score, chains = ReferenceDecoder().viterbi(ex.window, num)
        self.assertTrue(math.isfinite(score))
        self.assertEqual(len(chains), 1)
        c = chains[0]
        self.assertFalse(c.partial_5 or c.partial_3)
        self.assertEqual([(s.start, s.end) for s in c.cds()], ex.cds_ranges)
        self.assertEqual([(s.start, s.end) for s in c.introns()], ex.intron_ranges)

    def test_max_window_skips_and_counts(self):
        stats = LoaderStats()
        windows = list(iter_windows(self.summary, self.gff, self.fasta,
                                    max_window=50, stats=stats))
        self.assertEqual(windows, [])
        self.assertEqual(stats.skipped_too_long, 1)
        self.assertEqual(stats.yielded, 0)

    # -- source pins are enforced by the loading interface (finding 3) --------
    def test_iter_windows_enforces_source_pin(self):
        # A flank base changed after pinning: verify_source raises, and because
        # iter_windows gates on it, the loader yields nothing rather than a
        # window built from unpinned sequence.
        with open(self.fasta, "w") as fh:
            fh.write(_fasta().replace("CCCCCCCCCC", "GCCCCCCCCC", 1))
        with self.assertRaises(SourceMismatch):
            list(iter_windows(self.summary, self.gff, self.fasta))

    def test_table_is_bound_from_summary(self):
        # The summary, not a loader default, supplies the genetic-code id.
        with open(self.summary, "w") as fh:
            json.dump({"gff_md5": _md5(GFF), "fasta_md5": _md5(_fasta()),
                       "m": 20, "table": 1}, fh)
        ex = next(iter_windows(self.summary, self.gff, self.fasta))
        self.assertEqual(ex.table, 1)

    # -- soft masking survives orientation (finding 2) -----------------------
    def test_soft_mask_channel_preserved_plus(self):
        # An all-lower-case (soft-masked) contig must yield an all-lower-case
        # window, so the featurizer's soft-mask channel is fully populated;
        # oriented_chain's .upper() would have zeroed it.
        from model.a.features import _classify
        lower = _contig().lower()
        fasta = ">chr1\n" + "\n".join(lower[i:i + 70]
                                      for i in range(0, len(lower), 70)) + "\n"
        with open(self.fasta, "w") as fh:
            fh.write(fasta)
        with open(self.summary, "w") as fh:
            json.dump({"gff_md5": _md5(GFF), "fasta_md5": _md5(fasta)}, fh)
        ex = next(iter_windows(self.summary, self.gff, self.fasta))
        self.assertTrue(ex.window.islower())
        soft = sum(c[2] for c in _classify(ex.window))
        self.assertEqual(soft, ex.n)

    def test_oriented_window_preserves_case_both_strands(self):
        # _oriented_window keeps case on both strands (revcomp preserves case).
        from model.a.dataset import _oriented_window

        class _Stub:
            def __init__(self, strand):
                self.strand, self.span, self.rows = strand, (11, 20), []
        seq = _contig().lower()
        plus = _oriented_window(_Stub("+"), seq)
        minus = _oriented_window(_Stub("-"), seq)
        self.assertTrue(plus.islower() and minus.islower())
        self.assertEqual(minus, revcomp(plus))

    # -- a close neighbouring gene makes the window unclean (finding 1) -------
    def test_neighbouring_gene_skips_window(self):
        # Two nonoverlapping complete genes 11..19 and 21..29: each sits in the
        # other's flank, so neither window is clean.
        seq = "C" * 10 + "ATGAAATAA" + "C" + "ATGCCCTAA" + "A" * (CONTIG_LEN - 29)
        rows = ["##gff-version 3", "##sequence-region chr1 1 %d" % CONTIG_LEN]
        for i, (a, b) in enumerate(((11, 19), (21, 29)), 1):
            rows += [
                "chr1\tt\tgene\t%d\t%d\t.\t+\t.\tID=g%d;gene_biotype=protein_coding" % (a, b, i),
                "chr1\tt\tmRNA\t%d\t%d\t.\t+\t.\tID=t%d;Parent=g%d" % (a, b, i, i),
                "chr1\tt\tCDS\t%d\t%d\t.\t+\t0\tID=c%d;Parent=t%d" % (a, b, i, i)]
        gff = "\n".join(rows) + "\n"
        fasta = ">chr1\n" + seq + "\n"
        for path, text in ((self.gff, gff), (self.fasta, fasta)):
            with open(path, "w") as fh:
                fh.write(text)
        with open(self.summary, "w") as fh:
            json.dump({"gff_md5": _md5(gff), "fasta_md5": _md5(fasta)}, fh)
        stats = LoaderStats()
        windows = list(iter_windows(self.summary, self.gff, self.fasta, stats=stats))
        self.assertEqual(windows, [])
        self.assertEqual(stats.admitted, 2)
        self.assertEqual(stats.skipped_neighbor, 2)


if __name__ == "__main__":
    unittest.main()
