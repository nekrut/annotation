"""Parity tests for the tensor Viterbi (``model.a.fast_viterbi``).

The max-product scan is held to the standard-library delayed-entry decoder
(``model.grammar.delayed.DelayedEntryDecoder.viterbi``): identical best score
and identical ``Chain`` objects on the section-3.4 fixtures and on random
small lattices with IUPAC ambiguity, both genetic-code tables, alternative
initiators, ``m = 1..4``, ``R = 1..3`` and ``-inf`` masks; infeasible windows
are ``-inf`` in both; length-bucketed batched decoding equals per-window
decoding; the sparse-predecessor scan equals the dense-transition reference
bit for bit on scores, tail back-pointers and every traceback. Skips without
torch, like the other candidate-A tensor suites.
"""
import importlib.util
import random
import unittest
from math import isinf, log

from model.grammar import DurationMixture, TABLES
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
