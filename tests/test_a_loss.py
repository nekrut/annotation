"""Tests for candidate A's chain-loss oracle (proposal sections 3.1-3.2, 3.6).

All checks here use the standard-library ``ReferenceDecoder`` and run anywhere:
they pin the numerator support mask, the ``loss = log Z - log Z_num >= 0``
contract, that the masked numerator decodes back to exactly the gold chain, and
the gradient *signs* a training loss must have. The PyTorch loss (a later
T-human-014 increment) must reproduce these same numbers on the same fixtures;
the torch-gated test below is skipped until it exists.
"""
import importlib.util
import unittest
from math import isfinite, log

from model.grammar import DurationMixture, EdgePrior, ReferenceDecoder, Scores, TABLES
from model.a.loss import chain_nll, numerator_scores

# Gate the torch-parity test on both torch *and* the optional ``model.a.torch_loss``
# submodule existing. ``import torch`` raises ModuleNotFoundError when torch is
# absent; ``from model.a import torch_loss`` raises a plain ImportError (not a
# ModuleNotFoundError) once torch is present but the submodule has not been
# written, which the old ``except ModuleNotFoundError`` let abort discovery on a
# torch host. Probing with ``find_spec`` skips cleanly when the submodule is
# missing while still surfacing errors from a genuinely broken implementation.
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

# A short-intron duration so a legal intron fits a tiny hand-checkable window;
# the default m=20 would need a >=20 bp intron. Component weights are the
# section-3.4 placeholders.
SHORT = DurationMixture(m=2, pi=((0.5, 0.3, 0.2),) * 3, q=((0.9, 0.99, 0.999),) * 3)

# Single-exon gene: TT | ATG AAA AAA TAA | TT  -> M K K *  (no intron).
SINGLE_X = "TT" + "ATGAAAAAA" + "TAA" + "TT"
SINGLE_CDS = [(2, 14)]
SINGLE_INTRON: list = []

# Two-exon gene: TT | ATG AA | GTAC (intron) | A TAA | TT, spliced ATG AAA TAA.
TWO_X = "TT" + "ATGAA" + "GTAC" + "ATAA" + "TT"
TWO_CDS = [(2, 7), (11, 15)]
TWO_INTRON = [(7, 11)]


class SupportMask(unittest.TestCase):
    def test_mask_allows_only_gold_emissions(self):
        n = len(SINGLE_X)
        sc = numerator_scores(None, n, SINGLE_CDS, SINGLE_INTRON)
        NEG = float("-inf")
        # Outside the CDS only U may emit.
        for t in list(range(0, 2)) + list(range(14, n)):
            self.assertEqual(sc.u[t], 0.0)
            self.assertEqual(sc.cds[0][t], NEG)
            self.assertEqual(sc.start[t], NEG)
        # Inside the CDS U is forbidden and every phase channel is allowed.
        for t in range(2, 14):
            self.assertEqual(sc.u[t], NEG)
            for p in range(3):
                self.assertEqual(sc.cds[p][t], 0.0)
        # start only on the first CDS base, stop only on the last.
        self.assertEqual(sc.start[2], 0.0)
        self.assertEqual(sc.start[3], NEG)
        self.assertEqual(sc.stop[13], 0.0)
        self.assertEqual(sc.stop[12], NEG)
        # No intron here, so no donor/acceptor and no intron channel anywhere.
        for t in range(n):
            self.assertEqual(sc.donor[t], NEG)
            self.assertEqual(sc.acceptor[t], NEG)
            for p in range(3):
                self.assertEqual(sc.intron[p][t], NEG)

    def test_mask_places_intron_junctions(self):
        sc = numerator_scores(None, len(TWO_X), TWO_CDS, TWO_INTRON)
        NEG = float("-inf")
        # donor on the intron's first base, acceptor on the first CDS base after.
        self.assertEqual(sc.donor[7], 0.0)
        self.assertEqual(sc.acceptor[11], 0.0)
        for t in range(len(TWO_X)):
            if t != 7:
                self.assertEqual(sc.donor[t], NEG)
            if t != 11:
                self.assertEqual(sc.acceptor[t], NEG)
        # intron channels only on the four intronic bases.
        for t in range(7, 11):
            self.assertEqual(sc.intron[0][t], 0.0)
            self.assertEqual(sc.u[t], NEG)

    def test_base_scores_are_preserved_where_allowed_and_copied(self):
        base = Scores.zeros(len(SINGLE_X))
        base.cds[0][2] = 2.5      # gold: kept
        base.u[0] = 1.5           # gold: kept
        base.u[5] = 9.0           # off-gold (U inside CDS): masked to -inf
        sc = numerator_scores(base, len(SINGLE_X), SINGLE_CDS, SINGLE_INTRON)
        self.assertEqual(sc.cds[0][2], 2.5)
        self.assertEqual(sc.u[0], 1.5)
        self.assertEqual(sc.u[5], float("-inf"))
        # the mask must not mutate the caller's scores
        self.assertEqual(base.u[5], 9.0)


class Loss(unittest.TestCase):
    def setUp(self):
        self.dec = ReferenceDecoder(TABLES[1], SHORT)
        self.default = ReferenceDecoder(TABLES[1], DurationMixture())

    def test_single_exon_loss_is_log_two(self):
        # With zero emissions the only two legal paths are the gene and all-U,
        # each weight 0, so log Z = log 2 and the numerator (gene only) is 0.
        loss, z, zn = chain_nll(self.default, SINGLE_X, None, SINGLE_CDS, SINGLE_INTRON)
        self.assertAlmostEqual(zn, 0.0, places=9)
        self.assertAlmostEqual(z, log(2), places=9)
        self.assertAlmostEqual(loss, log(2), places=9)

    def test_loss_nonnegative_and_numerator_finite(self):
        for dec, x, cds, itr in [
            (self.default, SINGLE_X, SINGLE_CDS, SINGLE_INTRON),
            (self.dec, TWO_X, TWO_CDS, TWO_INTRON),
        ]:
            loss, z, zn = chain_nll(dec, x, None, cds, itr)
            self.assertTrue(isfinite(zn), "gold support must admit a legal path")
            self.assertGreaterEqual(loss, -1e-9)
            self.assertLessEqual(zn, z + 1e-9)

    def test_numerator_decodes_back_to_gold_chain(self):
        for dec, x, cds, itr in [
            (self.default, SINGLE_X, SINGLE_CDS, SINGLE_INTRON),
            (self.dec, TWO_X, TWO_CDS, TWO_INTRON),
        ]:
            sc = numerator_scores(None, len(x), cds, itr)
            v, chains = dec.viterbi(x, sc)
            self.assertTrue(isfinite(v))
            self.assertEqual(len(chains), 1)
            c = chains[0]
            self.assertFalse(c.partial_5 or c.partial_3)
            self.assertEqual([(s.start, s.end) for s in c.cds()], cds)
            self.assertEqual([(s.start, s.end) for s in c.introns()], itr)

    def test_gradient_signs(self):
        # Rewarding a gold-path emission lowers the loss; rewarding an emission
        # only a competing path can use (U inside the CDS) raises it. This is the
        # sign structure autograd must yield for the torch loss.
        base_loss, _, _ = chain_nll(self.dec, TWO_X, None, TWO_CDS, TWO_INTRON)
        reward_gold = Scores.zeros(len(TWO_X))
        reward_gold.cds[0][2] = 1.0                 # the gold start base
        lower, _, _ = chain_nll(self.dec, TWO_X, reward_gold, TWO_CDS, TWO_INTRON)
        self.assertLess(lower, base_loss)

        reward_off = Scores.zeros(len(TWO_X))
        reward_off.u[5] = 1.0                       # U on a CDS-interior base
        higher, _, _ = chain_nll(self.dec, TWO_X, reward_off, TWO_CDS, TWO_INTRON)
        self.assertGreater(higher, base_loss)

    def test_edge_enabled_decoder_is_rejected(self):
        # The support mask constrains emissions, not boundary states, so an
        # enabled EdgePrior would let the numerator claim extra entry/exit
        # hypotheses and quietly supervise the wrong path set. The complete-target
        # oracle rejects it; edge-partial targets are a section-3.6 increment.
        edged = ReferenceDecoder(TABLES[1], SHORT, edges=EdgePrior())
        with self.assertRaises(ValueError):
            chain_nll(edged, TWO_X, None, TWO_CDS, TWO_INTRON)

    def test_finite_difference_gradient_matches_marginal_difference(self):
        # d loss / d e_c[t] = P_free(emit c at t) - P_num(emit c at t). Rewarding
        # the gold start base sits only on the numerator support, so the marginal
        # difference is bounded in (-1, 0]; a central difference must land there.
        eps = 1e-4
        hi = Scores.zeros(len(TWO_X)); hi.cds[0][2] = eps
        lo = Scores.zeros(len(TWO_X)); lo.cds[0][2] = -eps
        lh, _, _ = chain_nll(self.dec, TWO_X, hi, TWO_CDS, TWO_INTRON)
        ll, _, _ = chain_nll(self.dec, TWO_X, lo, TWO_CDS, TWO_INTRON)
        grad = (lh - ll) / (2 * eps)
        self.assertGreater(grad, -1.0)
        self.assertLess(grad, 0.0)


@unittest.skipUnless(HAS_TORCH_LOSS, "torch chain loss not yet implemented")
class TorchParity(unittest.TestCase):
    def test_torch_matches_oracle(self):  # pragma: no cover - needs torch loss
        loss, _, _ = chain_nll(
            ReferenceDecoder(TABLES[1], SHORT), TWO_X, None, TWO_CDS, TWO_INTRON)
        t = torch_loss.chain_nll(TWO_X, None, TWO_CDS, TWO_INTRON, duration=SHORT)
        self.assertAlmostEqual(float(t), loss, places=5)


if __name__ == "__main__":
    unittest.main()
