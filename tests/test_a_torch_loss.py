"""Parity and autograd tests for candidate A's differentiable torch chain loss.

Every check here is gated on PyTorch: it asserts that ``model.a.torch_loss``
reproduces the standard-library ``model.a.loss`` oracle to floating-point
tolerance on the section-3.4 fixtures, and that autograd through the loss yields
the CRF marginal difference a training loop applies to the emission head. On a
torch-less host the whole module skips, exactly as the encoder and loss suites
gate their torch checks. The same fixtures the oracle pins (single-exon,
two-exon, both duration mixtures) are reused so the two implementations are held
to one specification.
"""
import importlib.util
import unittest
from math import isfinite, log

from model.grammar import DurationMixture, EdgePrior, ReferenceDecoder, Scores, TABLES
from model.a.loss import chain_nll as oracle_chain_nll, numerator_scores

try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ModuleNotFoundError:
    _HAS_TORCH = False

if _HAS_TORCH and importlib.util.find_spec("model.a.torch_loss") is not None:
    from model.a import torch_loss
    HAS_TORCH_LOSS = True
else:
    HAS_TORCH_LOSS = False

# Section-3.4 fixtures, identical to tests/test_a_loss.py.
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


@unittest.skipUnless(HAS_TORCH_LOSS, "torch chain loss not available")
class SupportMask(unittest.TestCase):
    def test_mask_matches_numerator_scores(self):
        # The additive torch mask must reproduce the oracle's -inf/keep pattern:
        # a 0 in the mask marks an allowed channel, -inf a forbidden one.
        for _dur, x, cds, itr in FIXTURES:
            n = len(x)
            m = torch_loss.support_mask(n, cds, itr, dtype=torch.float64)
            sc = numerator_scores(None, n, cds, itr)  # base zeros -> 0.0 where allowed
            channels = [sc.u, sc.cds[0], sc.cds[1], sc.cds[2],
                        sc.intron[0], sc.intron[1], sc.intron[2],
                        sc.start, sc.stop, sc.donor, sc.acceptor]
            for ci, arr in enumerate(channels):
                for t in range(n):
                    allowed = isfinite(arr[t])
                    self.assertEqual(bool(torch.isfinite(m[ci, t])), allowed,
                                     f"channel {ci} base {t}")

    def test_incomplete_and_overlapping_ranges_rejected(self):
        with self.assertRaises(ValueError):
            torch_loss.support_mask(10, [], [(2, 5)])
        with self.assertRaises(ValueError):
            torch_loss.support_mask(10, [(0, 5)], [(3, 6)])

    def test_empty_chain_is_the_all_intergenic_window(self):
        m = torch_loss.support_mask(10, [], [])
        self.assertTrue(bool(torch.isfinite(m[0]).all()))
        self.assertTrue(bool(torch.isinf(m[1:]).all()))


@unittest.skipUnless(HAS_TORCH_LOSS, "torch chain loss not available")
class Parity(unittest.TestCase):
    def test_loss_matches_oracle_zero_emissions(self):
        for dur, x, cds, itr in FIXTURES:
            loss, _z, _zn = oracle_chain_nll(ReferenceDecoder(TABLES[1], dur), x, None, cds, itr)
            t = torch_loss.chain_nll(x, None, cds, itr, duration=dur)
            self.assertAlmostEqual(float(t), loss, places=6)

    def test_partition_matches_oracle(self):
        for dur, x, _cds, _itr in FIXTURES:
            dec = ReferenceDecoder(TABLES[1], dur)
            z = dec.partition(x, Scores.zeros(len(x)))
            zt = torch_loss.partition(x, torch.zeros((11, len(x)), dtype=torch.float64), duration=dur)
            self.assertAlmostEqual(float(zt), z, places=6)

    def test_loss_matches_oracle_with_emissions(self):
        # A non-trivial emission tensor: reward the gold start base and penalise a
        # U emission inside the CDS, the same perturbations the oracle's sign test
        # uses. The two implementations must agree channel for channel.
        dur = SHORT
        base = Scores.zeros(len(TWO_X))
        base.cds[0][2] = 1.3
        base.u[5] = 0.7
        loss, _z, _zn = oracle_chain_nll(ReferenceDecoder(TABLES[1], dur), TWO_X, base,
                                         TWO_CDS, TWO_INTRON)
        e = torch.zeros((11, len(TWO_X)), dtype=torch.float64)
        e[1, 2] = 1.3   # cds[0] at base 2
        e[0, 5] = 0.7   # U at base 5
        t = torch_loss.chain_nll(TWO_X, e, TWO_CDS, TWO_INTRON, duration=dur)
        self.assertAlmostEqual(float(t), loss, places=6)

    def test_loss_nonnegative(self):
        for dur, x, cds, itr in FIXTURES:
            t = torch_loss.chain_nll(x, None, cds, itr, duration=dur)
            self.assertGreaterEqual(float(t), -1e-9)

    def test_single_exon_is_log_two(self):
        t = torch_loss.chain_nll(SINGLE_X, None, SINGLE_CDS, SINGLE_INTRON)
        self.assertAlmostEqual(float(t), log(2), places=6)


@unittest.skipUnless(HAS_TORCH_LOSS, "torch chain loss not available")
class Autograd(unittest.TestCase):
    def test_grad_matches_finite_difference_of_oracle(self):
        # d loss / d e_c[t] = P_free(emit c at t) - P_num(emit c at t). Autograd
        # through the torch loss must equal the oracle's central difference at the
        # gold start base (bounded in (-1, 0], since that channel sits only on the
        # numerator support).
        dur = SHORT
        e = torch.zeros((11, len(TWO_X)), dtype=torch.float64, requires_grad=True)
        loss = torch_loss.chain_nll(TWO_X, e, TWO_CDS, TWO_INTRON, duration=dur)
        loss.backward()
        g = float(e.grad[1, 2])   # cds[0] at gold start base 2

        eps = 1e-4
        hi = Scores.zeros(len(TWO_X)); hi.cds[0][2] = eps
        lo = Scores.zeros(len(TWO_X)); lo.cds[0][2] = -eps
        dec = ReferenceDecoder(TABLES[1], dur)
        lh, _, _ = oracle_chain_nll(dec, TWO_X, hi, TWO_CDS, TWO_INTRON)
        ll, _, _ = oracle_chain_nll(dec, TWO_X, lo, TWO_CDS, TWO_INTRON)
        fd = (lh - ll) / (2 * eps)
        self.assertAlmostEqual(g, fd, places=5)
        self.assertGreater(g, -1.0)
        self.assertLess(g, 0.0)


@unittest.skipUnless(HAS_TORCH_LOSS, "torch chain loss not available")
class Scope(unittest.TestCase):
    def test_bad_emission_shape_rejected(self):
        with self.assertRaises(ValueError):
            torch_loss.chain_nll(TWO_X, torch.zeros((11, 3), dtype=torch.float64),
                                 TWO_CDS, TWO_INTRON)
        with self.assertRaises(ValueError):
            torch_loss.TorchScores(torch.zeros((10, 5)))


@unittest.skipUnless(HAS_TORCH_LOSS, "torch chain loss not available")
class InputValidation(unittest.TestCase):
    """Both public entry points must restore the scalar oracle's
    finite-or-``-inf`` contract: ``_partition`` drops an ``isinf`` path, so a
    ``+inf`` would silently discard a legal transition and a ``NaN`` poison
    ``logsumexp``. The scalar oracle raises ``ValueError`` on each witness below;
    the torch loss must too, while still admitting legitimate ``-inf`` support."""

    ATG = "ATGTAA"        # single-exon complete chain, one CDS interval
    ATG_CDS = [(0, 6)]
    ATG_INTRON: list = []

    def _zeros(self):
        return torch.zeros((11, len(self.ATG)), dtype=torch.float64)

    def test_positive_infinity_rejected(self):
        # start channel (row 7) at base 0; the scalar oracle rejects +inf.
        for entry in (self._partition_entry, self._chain_entry):
            e = self._zeros(); e[7, 0] = float("inf")
            with self.assertRaises(ValueError):
                entry(e)

    def test_nan_on_used_channel_rejected(self):
        for entry in (self._partition_entry, self._chain_entry):
            e = self._zeros(); e[1, 0] = float("nan")   # cds[0] at base 0
            with self.assertRaises(ValueError):
                entry(e)

    def test_nan_on_unused_channel_rejected(self):
        # donor is never consumed by this single-exon chain, yet a NaN there must
        # still be rejected (the CCCCCC/donor[0]=NaN witness).
        for entry in (self._partition_entry, self._chain_entry):
            e = self._zeros(); e[9, 0] = float("nan")   # donor at base 0
            with self.assertRaises(ValueError):
                entry(e)

    def test_negative_infinity_preserved(self):
        # A hard mask (-inf) is legitimate support and must not be rejected.
        e = self._zeros(); e[0, 0] = float("-inf")      # forbid U at base 0
        torch_loss.partition(self.ATG, e)               # does not raise
        torch_loss.chain_nll(self.ATG, e, self.ATG_CDS, self.ATG_INTRON)

    def test_short_emissions_rejected(self):
        for entry in (self._partition_entry, self._chain_entry):
            with self.assertRaises(ValueError):
                entry(torch.zeros((11, 3), dtype=torch.float64))

    def test_long_emissions_rejected(self):
        for entry in (self._partition_entry, self._chain_entry):
            with self.assertRaises(ValueError):
                entry(torch.zeros((11, 7), dtype=torch.float64))

    def _partition_entry(self, e):
        return torch_loss.partition(self.ATG, e)

    def _chain_entry(self, e):
        return torch_loss.chain_nll(self.ATG, e, self.ATG_CDS, self.ATG_INTRON)


if __name__ == "__main__":
    unittest.main()
