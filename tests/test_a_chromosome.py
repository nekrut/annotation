"""Tests for the end-to-end chromosome path (``model.a.chromosome``): the
tiling partitions the sequence exactly with the declared overlap, core
de-duplication reports each chain once, chains map to genomic GFF3 on both
strands, and (with torch) a synthetic two-gene sequence decodes to the same
genes whether or not a window seam falls inside the overlap halo."""
import importlib.util
import os
import tempfile
import unittest

from model.grammar import Chain, Segment
from model.a import chromosome as C

try:
    import torch  # noqa: F401
    _HAS_TORCH = importlib.util.find_spec("model.a.pooled") is not None
except ModuleNotFoundError:
    _HAS_TORCH = False


class Tiles(unittest.TestCase):
    def test_cores_partition_and_windows_overlap(self):
        for n, window, overlap in [(1, 10, 0), (10, 10, 4), (11, 10, 4), (100, 10, 4),
                                   (100, 12, 5), (1000, 64, 0), (12289, 12288, 2048)]:
            layout = C.tiles(n, window, overlap)
            self.assertEqual(layout[0].start, 0)
            self.assertEqual(layout[-1].end, n)
            self.assertEqual(layout[0].core_start, 0)
            self.assertEqual(layout[-1].core_end, n)
            for t in layout:
                self.assertLessEqual(t.end - t.start, window)
                self.assertLessEqual(t.start, t.core_start)
                self.assertLessEqual(t.core_end, t.end)
            for a, b in zip(layout, layout[1:]):
                self.assertEqual(a.end, b.start + overlap)  # a is full length
                self.assertEqual(a.core_end, b.core_start)  # cores partition [0, n)
                self.assertEqual(b.core_start, b.start + overlap // 2)

    def test_short_gene_has_exactly_one_owner_that_contains_it(self):
        for n, window, overlap in [(500, 100, 40), (503, 100, 41), (12289, 12288, 4096)]:
            layout = C.tiles(n, window, overlap)
            half = overlap // 2
            longest = overlap - half
            for s in range(0, n, 7):
                for length in (1, longest // 2, longest):
                    e = s + length
                    if e > n:
                        continue
                    owners = [t for t in layout if t.core_start <= s < t.core_end]
                    self.assertEqual(len(owners), 1, (s, e))
                    t = owners[0]
                    self.assertTrue(t.start <= s and e <= t.end, (s, e, t))
                    if t.start > 0:
                        self.assertGreaterEqual(s - t.start, half)

    def test_rejects_bad_geometry(self):
        with self.assertRaises(ValueError):
            C.tiles(10, 0, 0)
        with self.assertRaises(ValueError):
            C.tiles(10, 10, 10)
        self.assertEqual(C.tiles(0, 10, 2), [])


def _gene(start, end):
    return Chain(segments=[Segment("cds", start, end)])


class Claiming(unittest.TestCase):
    def test_claimed_by_core_and_shifted(self):
        tile = C.Tile(100, 200, 120, 180)
        chains = [_gene(0, 9), _gene(20, 29), _gene(79, 90), _gene(80, 95)]
        got = C.claimed(chains, tile)
        self.assertEqual([(c.segments[0].start, c.segments[0].end) for c in got],
                         [(120, 129), (179, 190)])

    def test_gff3_rows_both_strands(self):
        n = 50
        # oriented [10, 19) on '-' is genomic [31, 40) one-based 32..40
        rows = C.predict_sequence.__globals__["gff3_rows"](_gene(10, 19), "chr", n, "-", "g1")
        cds = [r for r in rows if "\tCDS\t" in r][0].split("\t")
        self.assertEqual((cds[3], cds[4], cds[6]), ("32", "40", "-"))
        rows = C.predict_sequence.__globals__["gff3_rows"](_gene(10, 19), "chr", n, "+", "g1")
        cds = [r for r in rows if "\tCDS\t" in r][0].split("\t")
        self.assertEqual((cds[3], cds[4], cds[6]), ("11", "19", "+"))


@unittest.skipUnless(_HAS_TORCH, "needs torch")
class EndToEnd(unittest.TestCase):
    """A ``model`` whose encoder emits a strong hand-built signal for two
    intronless genes (one per strand) decodes both regardless of seam."""

    def _model_for(self, seq, genes):
        import torch
        from types import SimpleNamespace
        from model.a.encoder import DecoderParams
        from model.a.torch_loss import support_mask
        from model.grammar import reverse_complement

        n = len(seq)
        # Per-strand emission tables built from the gold support: +6 on the
        # gene's channels, everything else 0. The fake encoder looks up the
        # window by content (it sees oriented windows), so it is a dict.
        table = {}
        for strand in "+-":
            oriented = seq if strand == "+" else reverse_complement(seq)
            e = torch.zeros(11, n, dtype=torch.float64)
            for (a, b, gs) in genes:
                if gs != strand:
                    continue
                oa, ob = (a, b) if strand == "+" else (n - b, n - a)
                m = support_mask(n, [(oa, ob)], [])
                m = torch.where(torch.isinf(m), torch.zeros_like(m), torch.full_like(m, 6.0))
                # only inside the gene: keep U free elsewhere
                e[:, oa:ob] += m[:, oa:ob]
            table[oriented] = e

        def encoder(feats):
            # The encoder input is the stride-anchored slice ending at the
            # tile (``predict_sequence``); locate it by content in either
            # oriented strand and return the globally indexed scores.
            L = feats.shape[-1]
            for oriented, e in table.items():
                pos = oriented.find(self._current)
                if pos >= 0 and oriented.count(self._current) == 1:
                    out = torch.zeros(1, 11, L, dtype=torch.float64)
                    out[0, :, :len(self._current)] = e[:, pos:pos + len(self._current)]
                    return out
            raise AssertionError("unknown window")

        dec = DecoderParams().double()
        with torch.no_grad():
            dec.donor_dinuc.zero_(); dec.acceptor_dinuc.zero_()
        return SimpleNamespace(encoder=encoder, decoder=dec)

    def test_two_genes_both_strands_across_seams(self):
        import random
        import torch
        from model.a import pooled
        from model.grammar import TABLES, reverse_complement

        random.seed(7)
        n = 300
        bases = [random.choice("ACGT") for _ in range(n)]
        # plus gene at [40, 58): ATG + 4 codons + TAA ; minus gene occupying genomic [200, 218)
        plus = "ATG" + "GCA" * 4 + "TAA"
        bases[40:58] = list(plus)
        bases[200:218] = list(reverse_complement(plus))
        seq = "".join(bases)
        genes = [(40, 58, "+"), (200, 218, "-")]
        model = self._model_for(seq, genes)
        dec = model.decoder
        expected = {("+", 41, 58), ("-", 201, 218)}
        # The last three are the carried-seam mode: one segment scanned in
        # 50-base tiles (both genes cross a tile seam: 40-58 crosses 50, the
        # minus gene at oriented 82-100 crosses 100), two segments overlapping
        # by 60 (segment seam near 180, tile seams every 50 inside), and one
        # segment in tiles of 7 (every ``m``-window of a pending donor
        # straddles a seam).
        for window, overlap, segments in [(300, 0, None), (120, 60, None), (100, 70, None),
                                          (50, 0, 1), (50, 60, 2), (7, 0, 1)]:
            self.window, self.overlap = window, overlap
            # patch encode_sequence path: predict_sequence calls model.encoder(feats)
            # with the featurized window; remember the window text via a hook.
            from model.a import chromosome as CC
            orig = CC.predict_sequence

            def _run():
                rows_all = []
                counts_all = None
                import model.a.features as Fe
                real = Fe.encode_sequence

                def spy(padded, available=None):
                    self._current = padded[:sum(available)] if available else padded
                    return real(padded, available=available)
                Fe.encode_sequence = spy
                try:
                    rows, counts = orig(model, "chr", seq, code=TABLES[1],
                                        structure=pooled.structure(dec, 20),
                                        tables=pooled.duration_tables(dec),
                                        window=window, overlap=overlap, decode_batch=4,
                                        device=torch.device("cpu"), dtype=torch.float64,
                                        segments=segments)
                finally:
                    Fe.encode_sequence = real
                return rows, counts
            rows, counts = _run()
            found = {(r.split("\t")[6], int(r.split("\t")[3]), int(r.split("\t")[4]))
                     for r in rows if "\tCDS\t" in r}
            self.assertEqual(found, expected, (window, overlap, segments, rows))
            self.assertEqual(counts["chains"], 2, (window, overlap, segments))
            self.assertEqual(counts["strand_bases"], 2 * n)
            if segments is not None:
                self.assertEqual(counts["windows"], 2 * segments)
                if segments == 1:
                    self.assertEqual(counts["tiles"], 2 * -(-n // window))
                self.assertGreaterEqual(counts["tiles"], counts["windows"])

    @unittest.skipUnless(_HAS_TORCH, "torch")
    def test_real_encoder_pooling_grid_is_chromosome_anchored(self):
        """Two overlapping windows whose step is not a stride multiple give
        the same emissions to every shared interior base on both strands
        (stalin-0085: the pooling grid is anchored to the oriented origin,
        not reset per window)."""
        import random
        from unittest.mock import patch
        import torch
        from model.a import fast_viterbi, pooled
        from model.a.encoder import CandidateA, POOL_STRIDE
        from model.grammar import TABLES

        torch.set_num_threads(1)
        rng = random.Random(85)
        seq = "".join(rng.choices("ACGT", k=7400))
        torch.manual_seed(85)
        model = CandidateA().eval()
        window, overlap = 4800, 2200
        layout = C.tiles(len(seq), window, overlap)
        self.assertEqual(len(layout), 2)
        self.assertNotEqual(layout[1].start % POOL_STRIDE, 0)  # misaligned step
        seen = []

        def inspect(windows, emissions, **kw):
            start = layout[1].start
            lo, hi = start + 600, layout[0].end - 600
            delta = (emissions[0][:, lo:hi] - emissions[1][:, lo - start:hi - start]).abs()
            seen.append((hi - lo, float(delta.max())))
            return [(0.0, []) for _ in windows]

        with patch.object(fast_viterbi, "viterbi_windows", side_effect=inspect):
            _rows, counts = C.predict_sequence(
                model, "synthetic", seq, code=TABLES[1],
                structure=pooled.structure(model.decoder),
                tables=pooled.duration_tables(model.decoder),
                window=window, overlap=overlap, decode_batch=8,
                device=torch.device("cpu"), dtype=torch.float64)
        self.assertEqual(counts["windows"], 4)
        self.assertEqual(len(seen), 2)
        for compared, diff in seen:
            self.assertEqual(compared, 1000)
            self.assertLess(diff, 1e-6, seen)

    @unittest.skipUnless(_HAS_TORCH, "torch")
    def test_decode_batch_bounds_live_emissions(self):
        """Windows are flushed through the decoder in groups of
        ``decode_batch``: no call receives more, and no more emission tensors
        than that are alive at any call (engels-0084 P2). Encoder, featurizer
        and bias are stubbed; only the buffering is under test."""
        from types import SimpleNamespace
        from unittest.mock import patch
        import torch
        from model.a import fast_viterbi, features, pooled
        from model.a.encoder import DecoderParams
        from model.grammar import TABLES

        dec = DecoderParams().double()
        n, window, overlap, batch = 100_000, 12288, 4096, 2
        layout = C.tiles(n, window, overlap)
        seen = []

        def spy(windows, emissions, **kw):
            seen.append((len(windows), [e.shape[-1] for e in emissions],
                         len({e.untyped_storage().data_ptr() for e in emissions})))
            return [(0.0, []) for _ in windows]

        with patch.object(features, "encode_sequence",
                          side_effect=lambda x, available=None: torch.zeros(8, len(x))), \
             patch.object(pooled, "motif_bias",
                          side_effect=lambda x, dec, **kw: torch.zeros(11, len(x), dtype=kw["dtype"])), \
             patch.object(fast_viterbi, "viterbi_windows", side_effect=spy):
            _rows, counts = C.predict_sequence(
                SimpleNamespace(encoder=lambda f: torch.zeros(1, 11, f.shape[-1]), decoder=dec),
                "synthetic", "A" * n, code=TABLES[1],
                structure=pooled.structure(dec), tables=pooled.duration_tables(dec),
                window=window, overlap=overlap, decode_batch=batch,
                device=torch.device("cpu"), dtype=torch.float64)
        self.assertEqual(counts["windows"], 2 * len(layout))
        self.assertEqual(len(seen), 2 * -(-len(layout) // batch))
        self.assertTrue(all(k <= batch and alive <= batch for k, _l, alive in seen), seen)
        # every window decoded once, in tile order, on each strand
        lengths = [t.end - t.start for t in layout]
        got = [l for _k, ls, _a in seen for l in ls]
        self.assertEqual(got, lengths * 2)

    def test_write_gff3(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.gff3")
            nbytes = C.write_gff3(p, ["a\tb", "c\td"])
            with open(p) as fh:
                text = fh.read()
            self.assertEqual(text, "##gff-version 3\na\tb\nc\td\n")
            self.assertEqual(nbytes, len(text))


if __name__ == "__main__":
    unittest.main()
