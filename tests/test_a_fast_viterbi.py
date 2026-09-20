"""Parity tests for the tensor Viterbi (``model.a.fast_viterbi``).

The max-product scan is held to the standard-library delayed-entry decoder
(``model.grammar.delayed.DelayedEntryDecoder.viterbi``): identical best score
and identical ``Chain`` objects on the section-3.4 fixtures and on random
small lattices with IUPAC ambiguity, both genetic-code tables, alternative
initiators, ``m = 1..4``, ``R = 1..3`` and ``-inf`` masks; infeasible windows
are ``-inf`` in both; length-bucketed batched decoding equals per-window
decoding. Skips without torch, like the other candidate-A tensor suites.
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


if __name__ == "__main__":
    unittest.main()
