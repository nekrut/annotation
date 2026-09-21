"""Parity tests for the tensor Viterbi (``model.a.fast_viterbi``).

The max-product scan is held to the standard-library delayed-entry decoder
(``model.grammar.delayed.DelayedEntryDecoder.viterbi``): identical best score
and identical ``Chain`` objects on the section-3.4 fixtures and on random
small lattices with IUPAC ambiguity, both genetic-code tables, alternative
initiators, ``m = 1..4``, ``R = 1..3`` and ``-inf`` masks; infeasible windows
are ``-inf`` in both; length-bucketed batched decoding equals per-window
decoding; the sparse-predecessor scan equals the dense-transition reference
bit for bit on scores, tail back-pointers and every traceback; the
carried-seam scan over segments (``scan_segments`` / ``viterbi_segments``)
equals the unbroken scan at every tile length, in both modes. Skips without
torch, like the other candidate-A tensor suites.
"""
import importlib.util
import random
import unittest
from math import isinf, log

from model.grammar import DurationMixture, EdgePrior, TABLES
from model.grammar.delayed import DelayedEntryDecoder
from model.grammar.scores import Scores

try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ModuleNotFoundError:
    _HAS_TORCH = False

if _HAS_TORCH and importlib.util.find_spec("model.a.fast_viterbi") is not None:
    from model.a import fast_viterbi
    HAS_FAST = True
else:
    HAS_FAST = False

SHORT = DurationMixture(m=2, pi=((0.5, 0.3, 0.2),) * 3, q=((0.9, 0.99, 0.999),) * 3)
SINGLE_X = "TT" + "ATGAAAAAA" + "TAA" + "TT"
TWO_X = "TT" + "ATGAA" + "GTAC" + "ATAA" + "TT"


def _scores(e):
    e = e.tolist()
    return Scores(n=len(e[0]), u=e[0], cds=[e[1], e[2], e[3]], intron=[e[4], e[5], e[6]],
                  start=e[7], stop=e[8], donor=e[9], acceptor=e[10])


def _random_mixture(rng):
    m = rng.choice([1, 2, 3, 4])
    R = rng.choice([1, 2, 3])
    pi = [rng.random() + 0.1 for _ in range(R)]
    pi = tuple(p / sum(pi) for p in pi)
    q = tuple(rng.uniform(0.3, 0.95) for _ in range(R))
    return DurationMixture(m=m, pi=(pi,) * 3, q=(q,) * 3)


def _random_case(rng):
    n = rng.randint(4, 30)
    x = "".join(rng.choice("ACGTACGTNRY") for _ in range(n))
    code = rng.choice([TABLES[1], TABLES[6], TABLES[1].with_alternative_initiators()])
    e = torch.randn(11, n, dtype=torch.float64) * 2
    if rng.random() < 0.5:
        e[torch.rand(11, n) < 0.15] = float("-inf")
    return x, code, e


@unittest.skipUnless(HAS_FAST, "tensor Viterbi not available")
class FixtureParity(unittest.TestCase):
    def _check(self, x, e, code=TABLES[1], duration=DurationMixture()):
        a, ca = DelayedEntryDecoder(code, duration).viterbi(x, _scores(e))
        b, cb = fast_viterbi.viterbi(x, e, code=code, duration=duration)
        if isinf(a):
            self.assertTrue(isinf(b))
            self.assertEqual(cb, [])
            return a, ca
        self.assertAlmostEqual(a, b, places=9)
        self.assertEqual(ca, cb)
        return a, ca

    def test_zero_emissions_prefer_no_gene(self):
        for x in (SINGLE_X, TWO_X):
            best, chains = self._check(x, torch.zeros(11, len(x), dtype=torch.float64), duration=SHORT)
            self.assertEqual(best, 0.0)
            self.assertEqual(chains, [])

    def test_single_exon_recovered(self):
        e = torch.zeros(11, len(SINGLE_X), dtype=torch.float64)
        e[0, 2:14] = -5.0                                 # u disfavoured over the CDS
        best, chains = self._check(SINGLE_X, e)
        self.assertEqual(len(chains), 1)
        self.assertEqual([(s.start, s.end) for s in chains[0].cds()], [(2, 14)])
        self.assertEqual(best, 0.0)

    def test_two_exon_recovered(self):
        e = torch.zeros(11, len(TWO_X), dtype=torch.float64)
        e[0, 2:15] = -5.0
        e[4:7, 7:11] = 2.0                                # intron favoured
        e[1:4, 7:11] = -2.0
        best, chains = self._check(TWO_X, e, duration=SHORT)
        self.assertEqual(len(chains), 1)
        self.assertEqual([(s.start, s.end) for s in chains[0].cds()], [(2, 7), (11, 15)])
        self.assertEqual([(s.start, s.end) for s in chains[0].introns()], [(7, 11)])

    def test_ambiguity_prior_is_per_base(self):
        # The gene path through one N costs -log 4 in Viterbi (one base of
        # four), while the all-U path is free; the sum-product prior would
        # have priced the N at 0 and tied them.
        x = "ATGAANTAA"
        e = torch.zeros(11, len(x), dtype=torch.float64)
        best, chains = self._check(x, e)
        self.assertEqual(best, 0.0)
        e[0, :] = -1.0
        best, chains = self._check(x, e)
        self.assertAlmostEqual(best, -log(4.0), places=12)
        self.assertTrue(chains[0].uncertain)

    def test_infeasible_window(self):
        x = "GGGGG"
        e = torch.zeros(11, len(x), dtype=torch.float64)
        e[0, :] = float("-inf")                            # no U, no initiator
        self._check(x, e)

    def test_empty_window(self):
        best, chains = fast_viterbi.viterbi("", torch.zeros(11, 0, dtype=torch.float64))
        self.assertEqual((best, chains), (0.0, []))


@unittest.skipUnless(HAS_FAST, "tensor Viterbi not available")
class RandomLatticeParity(unittest.TestCase):
    def test_random_lattices(self):
        rng = random.Random(7)
        torch.manual_seed(7)
        finite = infeasible = 0
        for _ in range(120):
            dur = _random_mixture(rng)
            x, code, e = _random_case(rng)
            a, ca = DelayedEntryDecoder(code, dur).viterbi(x, _scores(e))
            b, cb = fast_viterbi.viterbi(x, e, code=code, duration=dur)
            if isinf(a):
                self.assertTrue(isinf(b), (x, dur))
                infeasible += 1
                continue
            finite += 1
            self.assertAlmostEqual(a, b, places=9, msg=(x, dur, code))
            self.assertEqual(ca, cb, (x, dur, code))
        self.assertGreater(finite, 40)
        self.assertGreater(infeasible, 10)


@unittest.skipUnless(HAS_FAST, "tensor Viterbi not available")
class Batching(unittest.TestCase):
    def test_batched_equals_single(self):
        rng = random.Random(3)
        torch.manual_seed(3)
        cases = [_random_case(rng) for _ in range(40)]
        ws, cs, es = zip(*cases)
        single = [fast_viterbi.viterbi(x, e, code=c, duration=SHORT) for x, c, e in cases]
        batched = fast_viterbi.viterbi_windows(ws, es, codes=cs, duration=SHORT, batch_size=7)
        for (a, ca), (b, cb) in zip(single, batched):
            if isinf(a):
                self.assertTrue(isinf(b))
            else:
                self.assertAlmostEqual(a, b, places=9)
            self.assertEqual(ca, cb)

    def test_sparse_scan_equals_dense_reference(self):
        """Scores, ``exit_r``/``entered`` and the traced states of every finite
        window agree exactly between ``viterbi_batch`` (sparse predecessors)
        and ``viterbi_batch_reference`` (dense ``(B, K, K)`` block), across
        codes, dtypes, ``m``, ``-inf`` masks and padded batches."""
        from model.a.fast_loss import Grammar, symbol_index_tensor
        rng = random.Random(88)
        torch.manual_seed(88)
        finite = 0
        for code in (TABLES[1], TABLES[6], TABLES[1].with_alternative_initiators()):
            for dtype in (torch.float32, torch.float64):
                for m in (1, 2, 4, 20):
                    for mask in (0.0, 0.12):
                        B, L = 6, 61
                        xs = ["".join(rng.choices("ACGTacgtNRY", k=L)) for _ in range(B)]
                        lengths = torch.tensor([L, L - 7, L // 2, L, 1, 0])
                        e = 2 * torch.randn(B, 11, L, dtype=dtype)
                        if mask:
                            e[torch.rand(B, 11, L) < mask] = float("-inf")
                        for b in range(B):
                            e[b, :, lengths[b]:] = float("-inf")
                        sym = torch.stack([symbol_index_tensor(x) for x in xs])
                        dur = _random_mixture(rng)
                        dur = DurationMixture(m=m, pi=dur.pi, q=dur.q)
                        g = Grammar.get(code, dur, dtype=dtype, device=e.device)
                        a = fast_viterbi.viterbi_batch(e, sym, lengths, g)
                        b = fast_viterbi.viterbi_batch_reference(e, sym, lengths, g)
                        self.assertTrue(torch.equal(a[0], b[0]), (code.table, dtype, m, mask))
                        self.assertTrue(torch.equal(a[2], b[2]) and torch.equal(a[3], b[3]))
                        for i in range(B):
                            if isinf(float(a[0][i])):
                                continue
                            finite += 1
                            n = int(lengths[i])
                            sa = fast_viterbi.traceback(n, g, a[1][i], a[2][i], a[3][i])
                            sb = fast_viterbi.traceback(n, g, b[1][i], b[2][i], b[3][i])
                            self.assertEqual(sa, sb)
        self.assertGreater(finite, 100)

    def test_bad_inputs_rejected(self):
        e = torch.zeros(11, 4, dtype=torch.float64)
        with self.assertRaises(ValueError):
            fast_viterbi.viterbi("ACG", e)
        with self.assertRaises(ValueError):
            fast_viterbi.viterbi_windows(["ACGT"], [e, e])
        with self.assertRaises(ValueError):
            fast_viterbi.viterbi_windows(["ACGT"], [e], batch_size=0)
        with self.assertRaises(ValueError):
            fast_viterbi.viterbi_batch(e[None, :5], torch.zeros(1, 4, dtype=torch.long),
                                       torch.tensor([4]),
                                       fast_viterbi.Grammar.get(TABLES[1], DurationMixture(),
                                                                dtype=e.dtype, device=e.device))


@unittest.skipUnless(HAS_FAST, "tensor Viterbi not available")
class SymbolLookup(unittest.TestCase):
    def test_symbol_index_tensor_equals_reference(self):
        from model.a.fast_loss import symbol_index_tensor, symbol_indices

        every_byte = "".join(chr(b) for b in range(128))
        rng = random.Random(11)
        cases = ["", "A", every_byte, every_byte[::-1],
                 "".join(rng.choices("ACGTacgtNnRYSWKMBDHVryswkmbdhvXxUu-*.", k=3001))]
        for x in cases:
            got = symbol_index_tensor(x)
            self.assertEqual(got.dtype, torch.long)
            self.assertEqual(got.tolist(), symbol_indices(x))
        with self.assertRaises(ValueError):
            symbol_index_tensor("AC\u00e9")

    def test_phase_tables_match_state_tables(self):
        from model.a.fast_loss import Grammar, duration_by_phase, duration_by_state

        g = Grammar.get(TABLES[1], SHORT, dtype=torch.float64, device="cpu")
        for by_phase, by_state in zip(duration_by_phase(g), duration_by_state(g)):
            self.assertTrue(torch.equal(by_phase[g.phase], by_state))


if __name__ == "__main__":
    unittest.main()


@unittest.skipUnless(HAS_FAST, "tensor Viterbi not available")
class EdgeParity(unittest.TestCase):
    """Sequence-edge partials (proposal 3.1) under the four-family
    ``EdgePrior``: the edge-enabled tensor scan is held to
    ``DelayedEntryDecoder(code, duration, edges).viterbi`` on random small
    lattices and on fixtures whose best path is a partial chain, batched and
    single, with the four scalars distinct so a family swap would show."""

    EDGES = (EdgePrior(log(0.3), log(0.4)),
             EdgePrior(entry=-0.2, exit=-1.7, intron_entry=-0.9, intron_exit=-0.05),
             EdgePrior(entry=2.0, exit=1.0, intron_entry=3.0, intron_exit=2.5))

    def _check(self, x, e, code, dur, edges, places=9):
        a, ca = DelayedEntryDecoder(code, dur, edges).viterbi(x, _scores(e))
        b, cb = fast_viterbi.viterbi(x, e, code=code, duration=dur, edges=edges)
        if isinf(a):
            self.assertTrue(isinf(b), (x, dur, edges))
            self.assertEqual(cb, [])
            return a, ca
        self.assertAlmostEqual(a, b, places=places, msg=(x, dur, code, edges))
        self.assertEqual(ca, cb, (x, dur, code, edges))
        return a, ca

    def test_random_lattices_with_edges(self):
        rng = random.Random(11)
        torch.manual_seed(11)
        finite = partial5 = partial3 = both = 0
        for _ in range(240):
            dur = _random_mixture(rng)
            x, code, e = _random_case(rng)
            edges = rng.choice(self.EDGES)
            a, ca = self._check(x, e, code, dur, edges)
            if isinf(a):
                continue
            finite += 1
            partial5 += any(c.partial_5 for c in ca)
            partial3 += any(c.partial_3 for c in ca)
            both += any(c.partial_5 and c.partial_3 for c in ca)
        self.assertGreater(finite, 150)
        self.assertGreater(partial5, 30)
        self.assertGreater(partial3, 30)
        self.assertGreater(both, 10)

    def test_families_are_distinct(self):
        # CDS from the edge (E0) pays the coding entry; a residual intron (J)
        # the intron entry; an exit in S/E the coding exit; in T the intron
        # exit. Each fixture's best path uses exactly one family, and the
        # score moves with that scalar only.
        dur = SHORT
        x = "AAAAAAAAA"
        e = torch.zeros(11, len(x), dtype=torch.float64)
        e[0, :] = -3.0                                    # U disfavoured everywhere
        e[4:7, :] = -10.0                                 # and so is every intronic base
        base = EdgePrior(entry=0.0, exit=0.0, intron_entry=0.0, intron_exit=0.0)
        best0, chains = self._check(x, e, TABLES[1], dur, base)
        self.assertEqual(len(chains), 1)
        self.assertTrue(chains[0].partial_5 and chains[0].partial_3)
        self.assertEqual([s.kind for s in chains[0].segments], ["cds"])
        for name, delta in (("entry", -1.5), ("exit", -0.5)):
            best, _ = self._check(x, e, TABLES[1], dur, EdgePrior(**{**base.__dict__, name: delta}))
            self.assertAlmostEqual(best, best0 + delta, places=12)
        for name in ("intron_entry", "intron_exit"):
            best, _ = self._check(x, e, TABLES[1], dur, EdgePrior(**{**base.__dict__, name: -1.5}))
            self.assertAlmostEqual(best, best0, places=12)
        # residual intron from the edge closing into CDS, then a censored intron
        e = torch.zeros(11, len(x), dtype=torch.float64)
        e[0, :] = -3.0
        e[4:7, :3] = 2.0                                  # intronic at the start
        e[1:4, :3] = -10.0
        e[4:7, 3:6] = -10.0                               # coding in the middle
        e[4:7, 6:] = 2.0                                  # intronic at the end
        e[1:4, 6:] = -10.0
        best0, chains = self._check(x, e, TABLES[1], dur, base)
        self.assertEqual([s.kind for s in chains[0].segments], ["intron", "cds", "intron"])
        for name, delta in (("intron_entry", -1.5), ("intron_exit", -0.5)):
            best, _ = self._check(x, e, TABLES[1], dur, EdgePrior(**{**base.__dict__, name: delta}))
            self.assertAlmostEqual(best, best0 + delta, places=12)
        for name in ("entry", "exit"):
            best, _ = self._check(x, e, TABLES[1], dur, EdgePrior(**{**base.__dict__, name: -1.5}))
            self.assertAlmostEqual(best, best0, places=12)

    def test_terminal_tail_survives_residual_intron(self):
        # engels-0088 / stalin-0090: the best path is E0 at base 0, a donor at
        # 1 and the tail T through the end (score 6 - log 3 - 4 log 2, closed
        # form). With J folded into the tail layer, a residual-intron entry
        # that outscored the donor entry erased the tail's only terminal
        # candidate and the scan fell back to CDS [0, 4). The score cannot
        # depend on the J weight since the winning path never uses J.
        x = "AAAA"
        e = torch.zeros(11, 4, dtype=torch.float64)
        e[0] = -10.0
        e[4:7] = 2.0
        e[4:7, 0] = 10.0
        e[7:11] = -100.0
        e[9, 1] = 0.0
        dur = DurationMixture(m=1)
        want = 6 - log(3) - 4 * log(2)
        for j in (-20.0, -10.0, -5.0, log(0.5), 0.0, 5.0):
            edges = EdgePrior(intron_entry=j)
            best, chains = self._check(x, e, TABLES[1], dur, edges)
            self.assertAlmostEqual(best, want, places=12, msg=j)
            self.assertEqual([(s.kind, s.start, s.end) for s in chains[0].segments],
                             [("cds", 0, 1), ("intron", 1, 4)])
            self.assertTrue(chains[0].partial_5 and chains[0].partial_3)

    def test_terminal_tail_boundary_sweep(self):
        # stalin-0090's boundary sweep: the intron reaches the tail at n = m + 1
        # and beyond; n <= m is a censored pending donor or single-base CDS.
        # Both J weights, codes 1 and 6, both dtypes, R = 1 and 3, batched
        # equal to single.
        mixes = {1: lambda m: DurationMixture(m=m, pi=((1.0,),) * 3, q=((0.6,), (0.7,), (0.8,))),
                 3: lambda m: DurationMixture(m=m, pi=((.2, .3, .5), (.5, .2, .3), (.3, .5, .2)),
                                              q=((.2, .5, .9), (.3, .7, .8), (.4, .6, .95)))}
        tails = 0
        for code in (TABLES[1], TABLES[6]):
            for dtype in (torch.float64, torch.float32):
                for m in (2, 4, 20):
                    for R, mk in mixes.items():
                        dur = mk(m)
                        cases = []
                        for n in (m - 1, m, m + 1, m + 2):
                            x = "A" * n
                            e = torch.zeros(11, n, dtype=dtype)
                            e[0] = -1000.0
                            e[7:11] = -1000.0
                            e[1:4, 1:] = -1000.0
                            e[4:7] = 2.0
                            e[4:7, 0] = 100.0
                            if n > 1:
                                e[9, 1] = 0.0
                            cases.append((x, e))
                        for j in (-1000.0, log(0.5)):
                            edges = EdgePrior(intron_entry=j)
                            single = [self._check(x, e, code, dur, edges,
                                                  places=9 if dtype is torch.float64 else 4)
                                      for x, e in cases]
                            tails += sum(any(sg.kind == "intron" and sg.end == len(x)
                                             for c in ch for sg in c.segments)
                                         for (x, _), (_, ch) in zip(cases, single))
                            for bs in (1, 7, 64):
                                batched = fast_viterbi.viterbi_windows(
                                    [x for x, _ in cases], [e for _, e in cases], codes=[code] * 4,
                                    duration=dur, batch_size=bs, edges=edges)
                                for (x, e), (a, ca), (b, cb) in zip(cases, single, batched):
                                    tb = fast_viterbi.viterbi(x, e, code=code, duration=dur, edges=edges)
                                    self.assertEqual(tb[0], b, (n, m, R, j))
                                    self.assertEqual(ca, cb)
        self.assertGreaterEqual(tails, 48)

    def test_no_edges_unchanged_and_empty(self):
        e = torch.zeros(11, len(SINGLE_X), dtype=torch.float64)
        e[0, 2:14] = -5.0
        for edges in self.EDGES:
            best, chains = self._check(SINGLE_X, e, TABLES[1], DurationMixture(), edges)
            if edges is self.EDGES[0]:                    # priors below log 0.5: the complete gene wins
                self.assertEqual(best, 0.0)
                self.assertEqual([(s.start, s.end) for s in chains[0].cds()], [(2, 14)])
                self.assertFalse(chains[0].partial_5 or chains[0].partial_3)
            else:                                         # parity holds whichever path wins
                self.assertGreaterEqual(best, 0.0)
        best, chains = fast_viterbi.viterbi("", torch.zeros(11, 0, dtype=torch.float64), edges=self.EDGES[1])
        self.assertEqual((best, chains), (0.0, []))

    def test_batched_with_edges_equals_single(self):
        rng = random.Random(5)
        torch.manual_seed(5)
        cases = [_random_case(rng) for _ in range(48)] + [("", TABLES[1], torch.zeros(11, 0, dtype=torch.float64))]
        ws, cs, es = zip(*cases)
        for edges in self.EDGES[1:]:
            single = [fast_viterbi.viterbi(x, e, code=c, duration=SHORT, edges=edges) for x, c, e in cases]
            for bs in (1, 7, 64):
                batched = fast_viterbi.viterbi_windows(ws, es, codes=cs, duration=SHORT, batch_size=bs, edges=edges)
                for (a, ca), (b, cb) in zip(single, batched):
                    if isinf(a):
                        self.assertTrue(isinf(b))
                    else:
                        self.assertEqual(a, b)
                    self.assertEqual(ca, cb)

    def test_learned_tables_and_pooled_prior(self):
        from model.a.encoder import DecoderParams
        from model.a import pooled
        torch.manual_seed(2)
        dec = DecoderParams()
        with torch.no_grad():
            for prm in dec.parameters():
                prm.add_(torch.randn_like(prm))
        edges = pooled.edge_prior(dec)
        self.assertEqual((edges.entry, edges.exit, edges.intron_entry, edges.intron_exit),
                         tuple(dec.partial_families.tolist()))
        tables, mix = pooled.duration_tables(dec), pooled.as_mixture(dec, m=2)
        rng = random.Random(9)
        finite = 0
        for _ in range(60):
            x, code, e = _random_case(rng)
            a, ca = DelayedEntryDecoder(code, mix, edges).viterbi(x, _scores(e))
            b, cb = fast_viterbi.viterbi(x, e, code=code, duration=pooled.structure(dec, m=2),
                                         tables=tables, edges=edges)
            if isinf(a):
                self.assertTrue(isinf(b))
                continue
            finite += 1
            self.assertAlmostEqual(a, b, places=9)
            self.assertEqual(ca, cb)
        self.assertGreater(finite, 30)


@unittest.skipUnless(HAS_FAST, "tensor Viterbi not available")
class SeamCarry(unittest.TestCase):
    """``viterbi_segments`` (the scan tiled with the state carried across
    seams, proposal 3.3) equals ``viterbi`` on the whole segment: same score
    bit for bit and the same ``Chain`` objects, for tiles of 1..64 bases
    (shorter than ``m``, so pending donors straddle seams), segments of
    unequal length in one batch (rows that end early idle through later
    tiles), empty segments, both genetic-code families, ``-inf`` masks, and
    with and without the sequence-edge partials (a tail entered before a
    seam must still be offered at the segment's end)."""

    EDGES = (None, EdgeParity.EDGES[1], EdgeParity.EDGES[2])

    def _segment_case(self, rng, nmax=60):
        n = rng.randint(0, nmax)
        x = "".join(rng.choice("ACGTACGTNRY") for _ in range(n))
        e = torch.randn(11, n, dtype=torch.float64) * 2
        if rng.random() < 0.5:
            e[torch.rand(11, n) < 0.15] = float("-inf")
        return x, e

    def test_tiled_equals_whole(self):
        rng = random.Random(3)
        torch.manual_seed(3)
        finite = seamed = 0
        for _ in range(200):
            dur = _random_mixture(rng)
            code = rng.choice([TABLES[1], TABLES[6], TABLES[1].with_alternative_initiators()])
            edges = rng.choice(self.EDGES)
            cases = [self._segment_case(rng) for _ in range(rng.randint(1, 6))]
            xs, es = zip(*cases)
            tile = rng.choice([1, 2, 3, 5, 7, 11, 64])
            tiled = fast_viterbi.viterbi_segments(xs, es, tile=tile, code=code, duration=dur, edges=edges)
            for (x, e), (b, cb) in zip(cases, tiled):
                a, ca = fast_viterbi.viterbi(x, e, code=code, duration=dur, edges=edges)
                if isinf(a):
                    self.assertTrue(isinf(b), (x, dur, code, edges, tile))
                    self.assertEqual(cb, [], (x, dur, code, edges, tile))
                    continue
                finite += 1
                seamed += any(s.start < k < s.end for c in ca for s in c.segments
                              for k in range(tile, len(x), tile))
                self.assertEqual(a, b, (x, dur, code, edges, tile))
                self.assertEqual(ca, cb, (x, dur, code, edges, tile))
        self.assertGreater(finite, 400)
        self.assertGreater(seamed, 100)

    def test_back_pointers_equal_unbroken_scan(self):
        # Beyond the traceback: the concatenated per-tile stores agree with
        # the single scan wherever the single scan's values are specified
        # (reachable states and boundaries inside each row's own length).
        from model.a.fast_loss import Grammar, symbol_index_tensor
        rng = random.Random(4)
        torch.manual_seed(4)
        for _ in range(40):
            dur = _random_mixture(rng)
            code = rng.choice([TABLES[1], TABLES[6]])
            edges = rng.choice(self.EDGES)
            cases = [self._segment_case(rng) for _ in range(rng.randint(1, 4))]
            lengths = [len(x) for x, _ in cases]
            L = max(lengths)
            g = Grammar.get(code, dur, dtype=torch.float64, device=torch.device("cpu"))
            e = torch.full((len(cases), 11, L), float("-inf"), dtype=torch.float64)
            sym = torch.zeros((len(cases), L), dtype=torch.long)
            for b, (x, em) in enumerate(cases):
                e[b, :, :len(x)] = em
                if x:
                    sym[b, :len(x)] = symbol_index_tensor(x)
            lens = torch.tensor(lengths)
            whole = fast_viterbi._scan(e, sym, lens, g, None, edges)
            for tile in (1, 3, 8):
                score, packed, final = fast_viterbi.scan_segments(e, sym, lens, g, tile, edges)
                self.assertTrue(torch.equal(whole[0], score), (tile, dur, edges))
                self.assertIsInstance(packed, fast_viterbi.PackedBackPointers)
                # 2 (prev slots) + K / 2 (exit_r nibbles) + ceil(K R / 8) (entered bits)
                self.assertEqual(packed.nbytes(), len(cases) * (L + 1) * (2 + (g.K + 1) // 2 + (g.K * g.R + 7) // 8))
                for b, n in enumerate(lengths):
                    prev, exit_r, entered = packed.dense(b)
                    # exit_r and entered are written at every step of the
                    # scan, so they are fully specified inside the row; prev
                    # is read for t >= 1.
                    self.assertTrue(torch.equal(whole[1][b, 1:n + 1], prev[1:n + 1]), (tile, b))
                    self.assertTrue(torch.equal(whole[2][b, :n], exit_r[:n]), (tile, b))
                    self.assertTrue(torch.equal(whole[3][b, :n + 1], entered[:n + 1]), (tile, b))
                    if edges is not None:
                        self.assertTrue(torch.equal(whole[4][b], final[b]), (tile, b))
                    if not isinf(float(whole[0][b])):
                        st_w = fast_viterbi.traceback(n, g, whole[1][b], whole[2][b], whole[3][b],
                                                      None if edges is None else whole[4][b])
                        st_p = fast_viterbi.traceback(n, g, prev, exit_r, entered,
                                                      None if edges is None else final[b])
                        self.assertEqual(st_w, st_p, (tile, b))

    def test_bad_tile_rejected(self):
        with self.assertRaises(ValueError):
            fast_viterbi.viterbi_segments(["ATG"], [torch.zeros(11, 3, dtype=torch.float64)], tile=0)
