"""Tests for the learned pooled decoder (``model.a.pooled``): the fast loss
and tensor Viterbi run under ``DecoderParams``' duration tables and
dinucleotide bias agree with the reference kernels run under the equivalent
concrete ``DurationMixture`` and biased emissions, and every one of the 50
consumed decoder scalars gets a finite gradient. Skips without torch, like the
other candidate-A tensor suites.
"""
import importlib.util
import random
import unittest
from math import isinf

from model.grammar import TABLES
from model.grammar.delayed import DelayedEntryDecoder

try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ModuleNotFoundError:
    _HAS_TORCH = False

if _HAS_TORCH and importlib.util.find_spec("model.a.pooled") is not None:
    from model.a import fast_loss, fast_viterbi, pooled, torch_loss
    from model.a.encoder import DecoderParams
    from tests.test_a_fast_viterbi import _random_case, _scores
    HAS_POOLED = True
else:
    HAS_POOLED = False

TWO_X = "TT" + "ATGAA" + "GTAC" + "ATAA" + "TT"
TWO_CDS = [(2, 7), (11, 15)]
TWO_INTRON = [(7, 11)]


def _random_decoder(seed):
    torch.manual_seed(seed)
    dec = DecoderParams()
    with torch.no_grad():
        for p in dec.parameters():
            p.normal_(0.0, 1.0)
    return dec


@unittest.skipUnless(HAS_POOLED, "pooled decoder not available")
class Tables(unittest.TestCase):
    def test_tables_match_mixture(self):
        dec = _random_decoder(1)
        t = pooled.duration_tables(dec)
        mix = pooled.as_mixture(dec, m=3)
        self.assertEqual(mix.m, 3)
        self.assertEqual(mix.R, 3)
        for p in range(3):
            for r in range(3):
                self.assertAlmostEqual(float(t.log_pi[p, r]), mix.log_pi(p, r), places=12)
                self.assertAlmostEqual(float(t.log_q[p, r]), mix.log_q(p, r), places=12)
                self.assertAlmostEqual(float(t.log_1mq[p, r]), mix.log_1mq(p, r), places=12)
        self.assertEqual(pooled.structure(dec, m=3).R, 3)
        self.assertEqual(pooled.structure(dec, m=3).m, 3)
        # The structure key is value-stable across parameter changes.
        self.assertEqual(pooled.structure(dec), pooled.structure(_random_decoder(2)))

    def test_fresh_decoder_is_uniform_mixture_with_distinct_hazards(self):
        mix = pooled.as_mixture(DecoderParams())
        self.assertEqual(mix.m, pooled.DEFAULT_M)
        for p in range(3):
            self.assertAlmostEqual(sum(mix.pi[p]), 1.0, places=12)
            self.assertTrue(all(abs(v - 1 / 3) < 1e-12 for v in mix.pi[p]))
            self.assertEqual(len(set(mix.q[p])), 3)

    def test_motif_bias_indexing(self):
        dec = _random_decoder(3)
        x = "GTAGNNAG"  # n = 8
        b = pooled.motif_bias(x, dec)
        self.assertEqual(tuple(b.shape), (11, 8))
        self.assertEqual(float(b[:9].abs().sum()), 0.0)  # only donor/acceptor rows
        d, a = dec.donor_dinuc.detach(), dec.acceptor_dinuc.detach()
        self.assertEqual(float(b[9, 0]), float(d[pooled.DINUC_INDEX["GT"]]))    # x[0:2]
        self.assertEqual(float(b[9, 3]), 0.0)                                   # "GN"
        self.assertEqual(float(b[9, 7]), 0.0)                                   # past the end
        self.assertEqual(float(b[10, 2]), float(a[pooled.DINUC_INDEX["GT"]]))   # x[0:2]
        self.assertEqual(float(b[9, 6]), float(d[pooled.DINUC_INDEX["AG"]]))    # x[6:8]
        self.assertEqual(float(b[10, 7]), 0.0)                                  # "NA"
        self.assertEqual(float(b[10, 4]), float(a[pooled.DINUC_INDEX["AG"]]))   # x[2:4]
        self.assertEqual(float(b[10, 0]), 0.0)                                  # before the start
        self.assertEqual(float(b[10, 1]), 0.0)
        # Batched form pads with zeros and matches the single form.
        bb = pooled.motif_bias_batch([x, "AG"], dec)
        self.assertEqual(tuple(bb.shape), (2, 11, 8))
        self.assertTrue(torch.equal(bb[0], b))
        self.assertEqual(float(bb[1, :, 2:].abs().sum()), 0.0)
        self.assertEqual(float(bb[1, 9, 0]), float(d[pooled.DINUC_INDEX["AG"]]))


@unittest.skipUnless(HAS_POOLED, "pooled decoder not available")
class LossParity(unittest.TestCase):
    def test_fixture_matches_reference_and_all_scalars_get_gradient(self):
        dec = _random_decoder(4)
        torch.manual_seed(4)
        e = torch.randn(11, len(TWO_X), dtype=torch.float64, requires_grad=True)
        biased = e + pooled.motif_bias(TWO_X, dec)
        fast = fast_loss.chain_nll(TWO_X, biased, TWO_CDS, TWO_INTRON,
                                   duration=pooled.structure(dec, m=2),
                                   tables=pooled.duration_tables(dec))
        ref = torch_loss.chain_nll(TWO_X, biased.detach(), TWO_CDS, TWO_INTRON,
                                   duration=pooled.as_mixture(dec, m=2))
        self.assertAlmostEqual(float(fast), float(ref), places=10)
        fast.backward()
        for name, p in dec.named_parameters():
            if name == "partial_families":
                self.assertIsNone(p.grad)  # not consumed until section 3.6
                continue
            self.assertIsNotNone(p.grad, name)
            self.assertTrue(torch.isfinite(p.grad).all(), name)
        # Mixture and hazard gradients are nonzero for the intron's phase.
        self.assertGreater(float(dec.mixture_logits.grad.abs().sum()), 0.0)
        self.assertGreater(float(dec.hazard_logits.grad.abs().sum()), 0.0)
        self.assertGreater(float(dec.donor_dinuc.grad.abs().sum()), 0.0)
        self.assertGreater(float(dec.acceptor_dinuc.grad.abs().sum()), 0.0)

    def test_gradient_matches_finite_differences(self):
        dec = _random_decoder(5).double()  # float32 parameters would quantise the step
        torch.manual_seed(5)
        e = torch.randn(11, len(TWO_X), dtype=torch.float64)

        def loss():
            biased = e + pooled.motif_bias(TWO_X, dec)
            return fast_loss.chain_nll(TWO_X, biased, TWO_CDS, TWO_INTRON,
                                       duration=pooled.structure(dec, m=2),
                                       tables=pooled.duration_tables(dec))

        loss().backward()
        h = 1e-6
        for name, p in dec.named_parameters():
            if name == "partial_families":
                continue
            flat = p.data.view(-1)
            for k in (0, flat.numel() - 1):
                old = float(flat[k])
                flat[k] = old + h
                up = float(loss())
                flat[k] = old - h
                down = float(loss())
                flat[k] = old
                self.assertAlmostEqual(float(p.grad.view(-1)[k]), (up - down) / (2 * h),
                                       places=5, msg=(name, k))

    def test_random_lattices_match_reference(self):
        rng = random.Random(6)
        dec = _random_decoder(6)
        finite = 0
        for _ in range(30):
            m = rng.choice([1, 2, 3])
            x, code, e = _random_case(rng)
            biased = e + pooled.motif_bias(x, dec)
            fast = fast_loss.partition(x, biased, code=code, duration=pooled.structure(dec, m),
                                       tables=pooled.duration_tables(dec))
            ref = torch_loss.partition(x, biased, code=code, duration=pooled.as_mixture(dec, m))
            if isinf(float(ref)):
                self.assertTrue(isinf(float(fast)))
                continue
            finite += 1
            self.assertAlmostEqual(float(fast), float(ref), places=9, msg=(x, m))
        self.assertGreater(finite, 10)


@unittest.skipUnless(HAS_POOLED, "pooled decoder not available")
class ViterbiParity(unittest.TestCase):
    def test_random_lattices_match_python_decoder(self):
        rng = random.Random(8)
        dec = _random_decoder(8)
        finite = 0
        cases = []
        for _ in range(40):
            m = rng.choice([1, 2, 3])
            x, code, e = _random_case(rng)
            biased = (e + pooled.motif_bias(x, dec)).detach()
            a, ca = DelayedEntryDecoder(code, pooled.as_mixture(dec, m)).viterbi(x, _scores(biased))
            b, cb = fast_viterbi.viterbi(x, biased, code=code, duration=pooled.structure(dec, m),
                                         tables=pooled.duration_tables(dec))
            if isinf(a):
                self.assertTrue(isinf(b))
            else:
                finite += 1
                self.assertAlmostEqual(a, b, places=9, msg=(x, m, code))
            self.assertEqual(ca, cb)
            if m == 2:
                cases.append((x, code, biased, b, cb))
        self.assertGreater(finite, 10)
        # Batched decoding under the learned tables equals the single calls.
        ws, cs, es, scores, chains = zip(*cases)
        batched = fast_viterbi.viterbi_windows(ws, es, codes=cs, duration=pooled.structure(dec, 2),
                                               tables=pooled.duration_tables(dec), batch_size=5)
        for (b, cb), s, ch in zip(batched, scores, chains):
            if isinf(s):
                self.assertTrue(isinf(b))
            else:
                self.assertAlmostEqual(b, s, places=9)
            self.assertEqual(cb, ch)

    def test_wrong_table_shape_rejected(self):
        dec = _random_decoder(9)
        e = torch.zeros(11, len(TWO_X), dtype=torch.float64)
        with self.assertRaises(ValueError):
            fast_viterbi.viterbi(TWO_X, e, tables=pooled.duration_tables(dec))  # R = 1 grammar


if __name__ == "__main__":
    unittest.main()
