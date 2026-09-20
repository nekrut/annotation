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
            L = feats.shape[-1]
            for oriented, e in table.items():
                for t in C.tiles(len(oriented), self.window, self.overlap):
                    x = oriented[t.start:t.end]
                    if self._current == x:
                        out = torch.zeros(1, 11, L, dtype=torch.float64)
                        out[0, :, :len(x)] = e[:, t.start:t.end]
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
        for window, overlap in [(300, 0), (120, 60), (100, 70)]:
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
                                        device=torch.device("cpu"), dtype=torch.float64)
                finally:
                    Fe.encode_sequence = real
                return rows, counts
            rows, counts = _run()
            found = {(r.split("\t")[6], int(r.split("\t")[3]), int(r.split("\t")[4]))
                     for r in rows if "\tCDS\t" in r}
            self.assertEqual(found, expected, (window, overlap, rows))
            self.assertEqual(counts["chains"], 2, (window, overlap))
            self.assertEqual(counts["strand_bases"], 2 * n)

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
