"""Parity tests for the vectorized delayed-entry chain loss (``model.a.fast_loss``).

The fast kernel is held to the *reference* torch forward (``model.a.torch_loss``,
itself pinned to the standard-library oracle): identical partition, identical
loss, identical emission gradient on the section-3.4 fixtures, on random small
lattices with ambiguous bases, ``m = 1``, multi-component mixtures, the ciliate
table and ``-inf`` masked channels, and batched execution equal to per-window
execution. Skips without torch, like the other candidate-A tensor suites.
"""
import importlib.util
import random
import unittest
from math import log

from model.grammar import DurationMixture, TABLES

try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ModuleNotFoundError:
    _HAS_TORCH = False

if _HAS_TORCH and importlib.util.find_spec("model.a.fast_loss") is not None:
    from model.a import fast_loss, torch_loss
    HAS_FAST = True
else:
    HAS_FAST = False

SHORT = DurationMixture(m=2, pi=((0.5, 0.3, 0.2),) * 3, q=((0.9, 0.99, 0.999),) * 3)
SINGLE_X = "TT" + "ATGAAAAAA" + "TAA" + "TT"
SINGLE_CDS = [(2, 14)]
SINGLE_INTRON: list = []
TWO_X = "TT" + "ATGAA" + "GTAC" + "ATAA" + "TT"
TWO_CDS = [(2, 7), (11, 15)]
TWO_INTRON = [(7, 11)]
FIXTURES = [
    (DurationMixture(), SINGLE_X, SINGLE_CDS, SINGLE_INTRON),
    (SHORT, SINGLE_X, SINGLE_CDS, SINGLE_INTRON),
    (SHORT, TWO_X, TWO_CDS, TWO_INTRON),
    # The gene-free (background) window: the empty chain, every base U.
    (SHORT, "TTATGAAGTAAGTAATT", [], []),
]


def _random_mixture(rng):
    m = rng.choice([1, 2, 3, 4])
    R = rng.choice([1, 2, 3])
    pi = [rng.random() + 0.1 for _ in range(R)]
    pi = tuple(p / sum(pi) for p in pi)
    q = tuple(rng.uniform(0.3, 0.95) for _ in range(R))
    return DurationMixture(m=m, pi=(pi,) * 3, q=(q,) * 3)


@unittest.skipUnless(HAS_FAST, "fast chain loss not available")
class FixtureParity(unittest.TestCase):
    def test_zero_emission_loss(self):
        for dur, x, cds, itr in FIXTURES:
            a = torch_loss.chain_nll(x, None, cds, itr, duration=dur)
            b = fast_loss.chain_nll(x, None, cds, itr, duration=dur)
            self.assertAlmostEqual(float(a), float(b), places=10)

    def test_single_exon_is_log_two(self):
        self.assertAlmostEqual(float(fast_loss.chain_nll(SINGLE_X, None, SINGLE_CDS, SINGLE_INTRON)),
                               log(2.0), places=10)

    def test_loss_and_gradient_with_emissions(self):
        torch.manual_seed(0)
        for dur, x, cds, itr in FIXTURES:
            e = torch.randn(11, len(x), dtype=torch.float64)
            e1 = e.clone().requires_grad_(True)
            e2 = e.clone().requires_grad_(True)
            a = torch_loss.chain_nll(x, e1, cds, itr, duration=dur)
            b = fast_loss.chain_nll(x, e2, cds, itr, duration=dur)
            self.assertAlmostEqual(float(a), float(b), places=9)
            a.backward()
            b.backward()
            self.assertLess(float((e1.grad - e2.grad).abs().max()), 1e-9)
            self.assertTrue(torch.isfinite(e2.grad).all())

    def test_partition_with_emissions(self):
        torch.manual_seed(1)
        for dur, x, cds, itr in FIXTURES:
            e = torch.randn(11, len(x), dtype=torch.float64)
            self.assertAlmostEqual(float(torch_loss.partition(x, e, duration=dur)),
                                   float(fast_loss.partition(x, e, duration=dur)), places=9)


@unittest.skipUnless(HAS_FAST, "fast chain loss not available")
class RandomLatticeParity(unittest.TestCase):
    """Random short sequences with IUPAC ambiguity, both genetic-code tables,
    every small ``m`` and ``R``, and sparse ``-inf`` masks: partition and
    gradient must match the reference forward, and an infeasible lattice must
    be ``-inf`` in both."""

    def test_random_lattices(self):
        rng = random.Random(7)
        torch.manual_seed(7)
        checked_infeasible = 0
        for _ in range(30):
            dur = _random_mixture(rng)
            n = rng.randint(4, 14)
            x = "".join(rng.choice("ACGTACGTNRY") for _ in range(n))
            code = TABLES[rng.choice([1, 6])]
            e = torch.randn(11, n, dtype=torch.float64) * 2
            if rng.random() < 0.5:
                e[torch.rand(11, n) < 0.15] = float("-inf")
            e1 = e.clone().requires_grad_(True)
            e2 = e.clone().requires_grad_(True)
            a = torch_loss.partition(x, e1, code=code, duration=dur)
            b = fast_loss.partition(x, e2, code=code, duration=dur)
            if not torch.isfinite(a):
                self.assertFalse(torch.isfinite(b))
                checked_infeasible += 1
                continue
            self.assertAlmostEqual(float(a), float(b), places=8, msg=(x, dur, code.table))
            a.backward()
            b.backward()
            self.assertLess(float((e1.grad - e2.grad).abs().max()), 1e-8)
            self.assertTrue(torch.isfinite(e2.grad).all())


@unittest.skipUnless(HAS_FAST, "fast chain loss not available")
class Batching(unittest.TestCase):
    WINDOWS = [SINGLE_X, TWO_X, "ATGAAATAA"]
    CHAINS = [(SINGLE_CDS, SINGLE_INTRON), (TWO_CDS, TWO_INTRON), ([(0, 9)], [])]

    def test_batched_equals_single(self):
        torch.manual_seed(2)
        L = max(len(w) for w in self.WINDOWS) + 3          # extra padding columns
        E = torch.randn(len(self.WINDOWS), 11, L, dtype=torch.float64, requires_grad=True)
        batched = fast_loss.batch_chain_nll(self.WINDOWS, E, self.CHAINS, duration=SHORT)
        batched.sum().backward()
        for b, (w, (cds, itr)) in enumerate(zip(self.WINDOWS, self.CHAINS)):
            e = E.detach()[b, :, : len(w)].clone().requires_grad_(True)
            single = torch_loss.chain_nll(w, e, cds, itr, duration=SHORT)
            single.backward()
            self.assertAlmostEqual(float(batched[b]), float(single), places=9)
            self.assertLess(float((E.grad[b, :, : len(w)] - e.grad).abs().max()), 1e-9)
            # padding columns receive no gradient
            self.assertEqual(float(E.grad[b, :, len(w):].abs().sum()), 0.0)

    def test_infeasible_support_is_inf_not_clamped(self):
        # A CDS interval that does not end on a stop admits no numerator path.
        v = fast_loss.chain_nll("TTATGAAATAATT", None, [(2, 10)], [])
        self.assertTrue(torch.isinf(v) and v > 0)

    def test_grammar_cache_keyed_by_genetic_code_value(self):
        # Standard table 1 and its alternative-initiator variant share a table
        # number; whichever is used first must not decide the other's grammar
        # (engels-0080 / stalin-0081, both call orders).
        x = "TTGAAATAA"
        e = torch.zeros(11, len(x), dtype=torch.float64)
        std, alt = TABLES[1], TABLES[1].with_alternative_initiators()
        for order in ((std, alt), (alt, std)):
            fast_loss.Grammar._cache.clear()
            for code in order:
                ref_z = float(torch_loss.partition(x, e, code=code))
                self.assertAlmostEqual(float(fast_loss.partition(x, e, code=code)), ref_z, places=12)
                ref_nll = float(torch_loss.chain_nll(x, e, [(0, 9)], [], code=code))
                got = float(fast_loss.chain_nll(x, e, [(0, 9)], [], code=code))
                if ref_nll == float("inf"):
                    self.assertEqual(got, float("inf"))
                else:
                    self.assertAlmostEqual(got, ref_nll, places=12)
            # the alternative code admits TTG..TAA (log Z = log 2), the standard does not
            self.assertAlmostEqual(float(fast_loss.partition(x, e, code=alt)), log(2), places=12)
            self.assertEqual(float(fast_loss.partition(x, e, code=std)), 0.0)

    def test_empty_window_matches_reference(self):
        e = torch.zeros(11, 0, dtype=torch.float64)
        self.assertEqual(float(fast_loss.partition("", e)), float(torch_loss.partition("", e)))
        self.assertEqual(float(fast_loss.partition("", e)), 0.0)

    def test_bad_inputs_rejected(self):
        with self.assertRaises(ValueError):
            fast_loss.chain_nll("ATG", torch.zeros(11, 5, dtype=torch.float64), [(0, 3)], [])
        e = torch.zeros(11, 3, dtype=torch.float64)
        e[0, 0] = float("nan")
        with self.assertRaises(ValueError):
            fast_loss.partition("ATG", e)
        with self.assertRaises(ValueError):
            fast_loss.batch_chain_nll(["ATGAAATAA"], torch.zeros(1, 11, 5, dtype=torch.float64),
                                      [([(0, 9)], [])])


if __name__ == "__main__":
    unittest.main()
