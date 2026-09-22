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
    annotated_gene_spans,
    coverage_report,
    format_coverage,
    gene_free_intervals,
    is_gene_feature,
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

    def test_is_gene_feature_types(self):
        for t in ("gene", "pseudogene", "ncRNA_gene", "V_gene_segment", "mRNA",
                  "lnc_RNA", "tRNA", "transcript", "exon", "CDS",
                  "five_prime_UTR", "three_prime_UTR", "pseudogenic_transcript"):
            self.assertTrue(is_gene_feature(t), t)
        for t in ("region", "chromosome", "centromere", "telomere",
                  "origin_of_replication", "long_terminal_repeat",
                  "mobile_genetic_element", "sequence_feature",
                  "regulatory_region", "biological_region"):
            self.assertFalse(is_gene_feature(t), t)

    def test_background_excludes_utrs_and_noncoding_genes(self):
        # engels-0096: the admission parser keeps only transcripts with CDS
        # rows and its span covers the CDS endpoints, so a UTR (mRNA 11..2000
        # around a 9-base CDS) and a separate lncRNA (3000..5000) must be
        # blocked from the raw GFF3 rows, not from the parsed transcripts.
        cds = "C" * 10 + "ATGAAATAA"
        seq = cds + ("aCgT" * 2500)[:CONTIG_LEN - len(cds)]
        fasta = ">chr1 utr fixture\n" + "\n".join(
            seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n"
        gff = "\n".join([
            "##gff-version 3",
            "##sequence-region chr1 1 %d" % CONTIG_LEN,
            "chr1\tt\tgene\t11\t2000\t.\t+\t.\tID=g1;gene_biotype=protein_coding",
            "chr1\tt\tmRNA\t11\t2000\t.\t+\t.\tID=t1;Parent=g1",
            "chr1\tt\texon\t11\t2000\t.\t+\t.\tParent=t1",
            "chr1\tt\tCDS\t11\t19\t.\t+\t0\tParent=t1",
            "chr1\tt\tgene\t3000\t5000\t.\t-\t.\tID=g2;gene_biotype=lncRNA",
            "chr1\tt\tlnc_RNA\t3000\t5000\t.\t-\t.\tID=t2;Parent=g2",
            "chr1\tt\texon\t3000\t5000\t.\t-\t.\tParent=t2",
            "chr1\tt\tcentromere\t6000\t6100\t.\t+\t.\tID=cen1",
            "",
        ])
        with open(self.gff, "w") as fh:
            fh.write(gff)
        with open(self.fasta, "w") as fh:
            fh.write(fasta)
        with open(self.summary, "w") as fh:
            json.dump({"gff_md5": _md5(gff), "fasta_md5": _md5(fasta),
                       "m": 20, "table": 1}, fh)
        spans = annotated_gene_spans(self.gff)
        self.assertEqual(sorted(set(spans["chr1"])), [(11, 19), (11, 2000), (3000, 5000)])
        # gene-free: [2010, 2989) and [5010, 10000) with FLANK=10
        self.assertEqual(gene_free_intervals(CONTIG_LEN, spans["chr1"]),
                         [(2010, 2989), (5010, 10000)])
        stats = LoaderStats()
        bg = list(iter_background_windows(self.summary, self.gff, self.fasta,
                                          length=1000, count=99, seed=1, stats=stats))
        # no tile fits [2010, 2989); [5010, 10000) holds four 1,000-base tiles
        self.assertEqual(stats.background_candidates, 4)
        self.assertEqual(sorted(e.key[2] for e in bg),
                         ["background:5010-6010", "background:6010-7010",
                          "background:7010-8010", "background:8010-9010"])
        for e in bg:
            a, b = (int(v) for v in e.key[2].split(":")[1].split("-"))
            self.assertFalse(a < 2000 or (a < 5000 and b > 2999), e.key)
            # the centromere is not a gene: tiles may cover it
        self.assertTrue(any(e.key[2] == "background:6010-7010" for e in bg) or
                        any(e.key[2] == "background:5010-6010" for e in bg))
        # source case is preserved on the forward strand, complemented on minus
        for e in bg:
            a = int(e.key[2].split(":")[1].split("-")[0])
            raw = seq[a:a + 1000]
            self.assertEqual(e.window, raw if e.key[1] == "+" else revcomp(raw))

    def test_load_all_windows_bounds_background_by_max_window(self):
        # engels-0096 P2, through the training loader itself: max_window=100
        # with a 1,000-base background window is refused; with max_window=1000
        # every loaded window (the 72-base chain and the background tile) fits.
        from model.a.train import TrainConfig, _load_all_windows
        src = {"name": "sp", "summary": self.summary, "gff": self.gff,
               "fasta": self.fasta, "dev_seqids": []}
        base = {"sources": [src], "out_dir": self.dir, "background_windows": 1,
                "background_length": 1000}
        with self.assertRaises(ValueError):
            TrainConfig.from_dict(dict(base, max_window=100))
        config = TrainConfig.from_dict(dict(base, max_window=1000))
        examples, stats = _load_all_windows(config)
        self.assertEqual(sorted(ex.n for ex in examples), [72, 1000])
        self.assertTrue(all(ex.n <= config.max_window for ex in examples))
        self.assertEqual(stats["sp"].background_yielded, 1)
        # the bound is re-checked on a config edited after parsing
        config.background_length = 4096
        with self.assertRaises(ValueError):
            _load_all_windows(config)

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

    # -- intergenic context clipped at neighbours (fit v2 loader increment) --
    def test_context_widens_window_and_shifts_ranges(self):
        # Gene at 11..62 on a 10 kb contig: 10 bases of sequence on the left,
        # so the 5' context is clipped to 0 by the sequence start (ws stays 1)
        # and the 3' side gets the full 100 bases.
        ex = next(iter_windows(self.summary, self.gff, self.fasta, context=100))
        self.assertEqual((ex.context_5, ex.context_3), (0, 100))
        self.assertEqual(ex.n, 172)
        self.assertEqual(ex.cds_ranges, [(10, 25), (50, 62)])
        self.assertEqual(ex.intron_ranges, [(25, 50)])
        self.assertEqual(ex.window[10:13], "ATG")
        self.assertEqual(ex.intergenic_bases, 172 - 27 - 25)
        # every added base is supervised U and the gold chain still decodes
        num = ex.support(None)
        score, chains = ReferenceDecoder().viterbi(ex.window, num)
        self.assertTrue(math.isfinite(score))
        self.assertEqual([(a.start, a.end) for a in chains[0].cds()], ex.cds_ranges)

    def test_context_stops_at_neighbouring_gene_both_strands(self):
        # Three complete genes: g1 at 101..109 (+), focal g2 at 301..309 (-),
        # g3 at 401..409 (+). Context 500 for g2 must stop at 110 on the left
        # and 400 on the right, i.e. add 181 genomic bases on the left and
        # 81 on the right; on the minus strand these are the 3' and 5' sides.
        seq = list("A" * CONTIG_LEN)
        seq[100:109] = "ATGAAATAA"
        seq[300:309] = revcomp("ATGCCCTAA")
        seq[400:409] = "ATGGGGTAA"
        seq = "".join(seq)
        rows = ["##gff-version 3", "##sequence-region chr1 1 %d" % CONTIG_LEN]
        for i, (a, b, st) in enumerate(((101, 109, "+"), (301, 309, "-"), (401, 409, "+")), 1):
            rows += [
                "chr1\tt\tgene\t%d\t%d\t.\t%s\t.\tID=g%d;gene_biotype=protein_coding" % (a, b, st, i),
                "chr1\tt\tmRNA\t%d\t%d\t.\t%s\t.\tID=t%d;Parent=g%d" % (a, b, st, i, i),
                "chr1\tt\tCDS\t%d\t%d\t.\t%s\t0\tID=c%d;Parent=t%d" % (a, b, st, i, i)]
        gff = "\n".join(rows) + "\n"
        fasta = ">chr1\n" + seq + "\n"
        for path, text in ((self.gff, gff), (self.fasta, fasta)):
            with open(path, "w") as fh:
                fh.write(text)
        with open(self.summary, "w") as fh:
            json.dump({"gff_md5": _md5(gff), "fasta_md5": _md5(fasta)}, fh)
        stats = LoaderStats()
        by_key = {ex.key: ex for ex in iter_windows(
            self.summary, self.gff, self.fasta, stats=stats, context=500)}
        self.assertEqual(stats.yielded, 3)
        g2 = by_key[("chr1", "-", "t2")]
        # flank window 291..319; context left to 110, right to 400
        self.assertEqual((g2.context_5, g2.context_3), (81, 181))
        self.assertEqual(g2.n, 400 - 110 + 1)
        self.assertEqual(g2.cds_ranges, [(81 + 10, 81 + 19)])
        self.assertEqual(g2.window[91:94], "ATG")
        self.assertEqual(g2.window[91:100], "ATGCCCTAA")
        g1 = by_key[("chr1", "+", "t1")]     # left: sequence start, right: g2
        self.assertEqual((g1.context_5, g1.context_3), (90, 181))
        self.assertEqual(g1.window[100:109], "ATGAAATAA")
        g3 = by_key[("chr1", "+", "t3")]     # left: g2, right: full 500
        self.assertEqual((g3.context_5, g3.context_3), (81, 500))
        self.assertEqual(stats.context_bases, 81 + 181 + 90 + 181 + 81 + 500)
        self.assertEqual(stats.context_clipped, 3)
        # no context: unchanged flank-only windows
        plain = {ex.key: ex for ex in iter_windows(self.summary, self.gff, self.fasta)}
        for k, ex in plain.items():
            self.assertEqual((ex.context_5, ex.context_3), (0, 0))
            self.assertEqual(ex.n, 29)
            self.assertEqual(ex.cds_ranges, [(10, 19)])

    def test_context_counts_toward_max_window(self):
        stats = LoaderStats()
        got = list(iter_windows(self.summary, self.gff, self.fasta, stats=stats,
                                context=100, max_window=100))
        self.assertEqual(got, [])
        self.assertEqual(stats.skipped_too_long, 1)
        with self.assertRaises(ValueError):
            list(iter_windows(self.summary, self.gff, self.fasta, context=-1))


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
