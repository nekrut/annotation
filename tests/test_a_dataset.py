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
    CoverageRow,
    LoaderStats,
    SourceMismatch,
    WindowExample,
    coverage_report,
    format_coverage,
    gene_free_intervals,
    iter_background_windows,
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

    def test_gene_free_intervals_widen_spans_by_margin(self):
        # gene 11..62 (1-based inclusive) widened by 10: blocks [0, 72)
        self.assertEqual(gene_free_intervals(100, [(11, 62)], margin=10), [(72, 100)])
        # overlapping/adjacent spans merge; a leading gap is kept
        self.assertEqual(gene_free_intervals(100, [(30, 40), (35, 50), (60, 61)], margin=0),
                         [(0, 29), (50, 59), (61, 100)])
        self.assertEqual(gene_free_intervals(100, [], margin=5), [(0, 100)])
        self.assertEqual(gene_free_intervals(0, [], margin=5), [])

    def test_background_windows_are_gene_free_and_seeded(self):
        stats = LoaderStats()
        bg = list(iter_background_windows(self.summary, self.gff, self.fasta,
                                          length=1000, count=3, seed=1, stats=stats))
        # gene-free tail is [72, 10000): nine 1,000-base tiles from 72
        self.assertEqual(stats.background_candidates, 9)
        self.assertEqual(stats.background_yielded, 3)
        self.assertEqual(len(bg), 3)
        for ex in bg:
            self.assertEqual(ex.n, 1000)
            self.assertEqual((ex.cds_ranges, ex.intron_ranges), ([], []))
            seqid, strand, tag = ex.key
            self.assertEqual(seqid, "chr1")
            self.assertIn(strand, "+-")
            a, b = (int(v) for v in tag.split(":")[1].split("-"))
            self.assertGreaterEqual(a, 72)
            self.assertEqual(b - a, 1000)
            self.assertEqual(ex.window, "A" * 1000 if strand == "+" else "T" * 1000)
            # the empty chain is the all-intergenic support: finite nll
            loss, log_z, log_z_num = chain_nll(
                ReferenceDecoder(), ex.window, None, ex.cds_ranges, ex.intron_ranges)
            self.assertTrue(math.isfinite(log_z) and math.isfinite(log_z_num))
            self.assertGreaterEqual(loss, -1e-9)
        # deterministic under the seed; a different seed draws differently
        again = list(iter_background_windows(self.summary, self.gff, self.fasta,
                                             length=1000, count=3, seed=1))
        self.assertEqual([e.key for e in again], [e.key for e in bg])
        # count above the candidate pool yields every tile once
        every = list(iter_background_windows(self.summary, self.gff, self.fasta,
                                             length=1000, count=50))
        self.assertEqual(len(every), 9)
        self.assertEqual(len({e.key[2] for e in every}), 9)
        # chain-window accounting untouched (dev reservations key off it)
        self.assertEqual(dict(stats.windows_by_seqid), {})

    def test_background_windows_enforce_source_pin(self):
        with open(self.fasta, "a") as fh:
            fh.write("\n")
        with self.assertRaises(SourceMismatch):
            list(iter_background_windows(self.summary, self.gff, self.fasta,
                                         length=100, count=1))

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


    # -- coverage accounting across a panel (engels-0065, stalin-0068) --------
    def _write_species(self, name, gff, fasta):
        d = tempfile.mkdtemp()
        g = os.path.join(d, "%s.gff3" % name)
        f = os.path.join(d, "%s.fna" % name)
        s = os.path.join(d, "%s.summary.json" % name)
        with open(g, "w") as fh:
            fh.write(gff)
        with open(f, "w") as fh:
            fh.write(fasta)
        with open(s, "w") as fh:
            json.dump({"gff_md5": _md5(gff), "fasta_md5": _md5(fasta)}, fh)
        return (name, s, g, f)

    def _neighbour_species(self):
        # Two nonoverlapping complete genes each in the other's flank: admitted
        # 2, yielded 0, skipped_neighbor 2 (mirrors test_neighbouring_gene).
        seq = "C" * 10 + "ATGAAATAA" + "C" + "ATGCCCTAA" + "A" * (CONTIG_LEN - 29)
        rows = ["##gff-version 3", "##sequence-region chr1 1 %d" % CONTIG_LEN]
        for i, (a, b) in enumerate(((11, 19), (21, 29)), 1):
            rows += [
                "chr1\tt\tgene\t%d\t%d\t.\t+\t.\tID=g%d;gene_biotype=protein_coding" % (a, b, i),
                "chr1\tt\tmRNA\t%d\t%d\t.\t+\t.\tID=t%d;Parent=g%d" % (a, b, i, i),
                "chr1\tt\tCDS\t%d\t%d\t.\t+\t0\tID=c%d;Parent=t%d" % (a, b, i, i)]
        return self._write_species("neigh", "\n".join(rows) + "\n", ">chr1\n" + seq + "\n")

    def test_coverage_report_aggregates_panel(self):
        clean = self._write_species("clean", GFF, _fasta())     # admitted 1, yielded 1
        neigh = self._neighbour_species()                        # admitted 2, yielded 0
        rows = coverage_report([clean, neigh])
        self.assertEqual([r.species for r in rows], ["clean", "neigh"])
        by = {r.species: r for r in rows}
        self.assertEqual((by["clean"].admitted, by["clean"].yielded), (1, 1))
        self.assertEqual(by["clean"].yielded_fraction, 1.0)
        self.assertEqual((by["neigh"].admitted, by["neigh"].yielded), (2, 0))
        self.assertEqual(by["neigh"].skipped_neighbor, 2)
        self.assertEqual(by["neigh"].yielded_fraction, 0.0)

    def test_coverage_report_honours_max_window(self):
        clean = self._write_species("clean", GFF, _fasta())
        rows = coverage_report([clean], max_window=50)
        self.assertEqual(rows[0].yielded, 0)
        self.assertEqual(rows[0].skipped_too_long, 1)

    def test_format_coverage_has_header_and_total(self):
        rows = coverage_report([
            self._write_species("clean", GFF, _fasta()),
            self._neighbour_species(),
        ])
        table = format_coverage(rows)
        lines = table.rstrip("\n").split("\n")
        self.assertEqual(lines[0].split("\t")[0], "species")
        self.assertEqual(len(lines), 1 + len(rows) + 1)   # header + rows + TOTAL
        total = lines[-1].split("\t")
        self.assertEqual(total[0], "TOTAL")
        self.assertEqual(total[1], "3")                   # 1 + 2 admitted
        self.assertEqual(total[2], "1")                   # 1 + 0 yielded
        self.assertEqual(total[3], "33.3")                # 1/3 yielded
        self.assertEqual(total[5], "2")                   # skip_neighbor column

    def test_coverage_report_enforces_source_pin(self):
        name, s, g, f = self._write_species("clean", GFF, _fasta())
        with open(f, "a") as fh:
            fh.write("A\n")               # tamper after pinning
        with self.assertRaises(SourceMismatch):
            coverage_report([(name, s, g, f)])


if __name__ == "__main__":
    unittest.main()
