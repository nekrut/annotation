"""Tests for candidate A's encoder (proposal section 3.5).

The arithmetic and featurizer tests use the standard library only and run
anywhere. The torch tests instantiate the module and check that the *built*
parameter total, the emission shape and the dependency radius agree with the
spec; they are skipped where torch is not installed.
"""
import unittest

from model.a import SECTION_35_PARAM_COUNT, section_35_param_breakdown
from model.a.features import GC_WINDOW, encode_sequence, gc_track

try:
    import torch
    from model.a.encoder import CandidateA, DNAEncoder, DEPENDENCY_RADIUS
    HAS_TORCH = True
except ModuleNotFoundError:
    HAS_TORCH = False


class Inventory(unittest.TestCase):
    def test_total_matches_spec(self):
        # Proposal section 3.5 fixes the inventory at exactly 455,841 scalars.
        self.assertEqual(SECTION_35_PARAM_COUNT, 455_841)

    def test_breakdown_rows(self):
        b = section_35_param_breakdown()
        self.assertEqual(b["stem_convolution"], 1_168)
        self.assertEqual(b["residual_blocks"], 1_392)
        self.assertEqual(b["pooled_context_projection"], 1_632)
        self.assertEqual(b["attention_blocks"], 447_616)
        self.assertEqual(b["fine_context_fusion"], 3_616)
        self.assertEqual(b["emission_projection"], 363)
        self.assertEqual(b["pooled_decoder"], 54)
        self.assertEqual(b["total"], 455_841)


class Featurizer(unittest.TestCase):
    def test_gc_track_ignores_ambiguous_and_centres_window(self):
        self.assertEqual(GC_WINDOW, 129)
        self.assertEqual(gc_track("GGGGG"), [1.0] * 5)          # all G/C
        self.assertEqual(gc_track("NNNNN"), [0.0] * 5)          # no unambiguous
        self.assertEqual(gc_track("GNC", window=3), [1.0] * 3)  # N excluded
        self.assertEqual(gc_track("ATGC", window=129), [0.5] * 4)

    def test_gc_track_excludes_unavailable_positions(self):
        # An unavailable neighbour must not count toward a real base's window,
        # whatever letter the padding invented (engels-0059 P2). The real A has
        # no unambiguous available neighbour, so its GC is the fixed 0.0.
        avail = [True, False]
        for padded in ("AC", "AG", "AT", "AN"):
            self.assertEqual(gc_track(padded, available=avail)[0], 0.0)
        # With the padded base marked available it does count (GC of A,C = 0.5).
        self.assertEqual(gc_track("AC", available=[True, True])[0], 0.5)


@unittest.skipUnless(HAS_TORCH, "torch not installed")
class BuiltModule(unittest.TestCase):
    def test_exactly_455841_parameters(self):
        self.assertEqual(CandidateA().num_parameters(), 455_841)

    def test_under_the_5m_cap(self):
        self.assertLess(CandidateA().num_parameters(), 5_000_000)

    def test_forward_eleven_channels_at_core_resolution(self):
        model = DNAEncoder().eval()
        L = 12 * 32  # multiple of the pooling stride; 32 context tokens
        with torch.no_grad():
            y = model(torch.zeros(2, 8, L))
        self.assertEqual(tuple(y.shape), (2, 11, L))

    def test_chunk_length_must_be_multiple_of_pool_stride(self):
        with self.assertRaises(ValueError):
            DNAEncoder()(torch.zeros(1, 8, 12 * 4 + 1))

    def test_featurizer_channels_and_padding(self):
        feats = encode_sequence("ACgtN")
        self.assertEqual(tuple(feats.shape), (8, 5))
        self.assertEqual(feats[0, 0].item(), 1.0)  # A
        self.assertEqual(feats[1, 1].item(), 1.0)  # C
        self.assertEqual(feats[2, 2].item(), 1.0)  # g -> G indicator
        self.assertEqual(feats[5, 2].item(), 1.0)  # g soft-masked
        self.assertEqual(feats[3, 3].item(), 1.0)  # t -> T indicator
        self.assertEqual(feats[4, 4].item(), 1.0)  # N ambiguity
        self.assertEqual(feats[0:4, 4].sum().item(), 0.0)  # N sets no ACGT
        self.assertEqual([v.item() for v in feats[7]], [1.0] * 5)  # all real
        padded = encode_sequence("AC", available=[True, False])
        self.assertEqual(padded[:, 1].abs().sum().item(), 0.0)

    def test_padding_letter_does_not_change_real_gc(self):
        # The GC channel at the real base is invariant to the invented padding
        # letter, and matches the unavailable-excluded 0.0 (engels-0059 P2).
        vals = {
            seq: encode_sequence(seq, available=[True, False])[6, 0].item()
            for seq in ("AC", "AT", "AN")
        }
        self.assertEqual(set(vals.values()), {0.0})

    def test_vectorised_featurizer_equals_reference(self):
        # The torch featurizer must be bit-identical to the per-base Python
        # reference (a-pilot 3.2 item 6 revision step 1) over every character
        # class it can meet: both cases of ACGT, N, IUPAC codes, non-letters,
        # with and without an availability mask, including the empty window.
        import random
        from model.a.features import encode_sequence_reference
        rng = random.Random(1406)
        alphabet = "ACGTacgtNnRYKMSWrykmsw-*"
        cases = [("", None), ("A", None), ("n", [False]), ("AC", [True, False])]
        for trial in range(40):
            n = rng.randint(1, 3 * GC_WINDOW)
            seq = "".join(rng.choices(alphabet, k=n))
            avail = None if trial % 2 else [rng.random() < 0.85 for _ in range(n)]
            cases.append((seq, avail))
        for seq, avail in cases:
            got, want = encode_sequence(seq, avail), encode_sequence_reference(seq, avail)
            self.assertEqual(got.dtype, torch.float32)
            self.assertTrue(torch.equal(got, want), (seq, avail))

    def test_featurizer_rejects_non_ascii_and_bad_mask(self):
        with self.assertRaises(ValueError):
            encode_sequence("AC\u00e9")
        with self.assertRaises(ValueError):
            encode_sequence("ACG", available=[True, True])

    def test_local_attention_matches_dense_oracle(self):
        # The windowed forward must equal the masked dense product it replaces,
        # in both value and gradient (engels-0059 P2).
        from model.a.encoder import LocalAttention
        from model.a.inventory import ATTN_OFFSETS, CONTEXT_WIDTH, N_HEADS
        torch.manual_seed(0)
        attn = LocalAttention(CONTEXT_WIDTH, N_HEADS, ATTN_OFFSETS)
        x = torch.randn(2, 40, CONTEXT_WIDTH, requires_grad=True)
        y_win = attn(x)
        y_win.sum().backward()
        g_win = x.grad.clone()
        x2 = x.detach().clone().requires_grad_(True)
        y_dense = attn._dense_forward(x2)
        y_dense.sum().backward()
        self.assertTrue(torch.allclose(y_win, y_dense, atol=1e-6))
        self.assertTrue(torch.allclose(g_win, x2.grad, atol=1e-6))

    def test_duration_components_start_distinct(self):
        # Symmetric hazards collapse the mixture to one geometric: mixture-logit
        # gradients are zero and hazard-logit gradients are equal across
        # components, so the components cannot differentiate under symmetric
        # updates (stalin-0062, engels-0060); the seed must break that.
        from model.a.encoder import DecoderParams
        h = DecoderParams().hazard_logits
        self.assertEqual(tuple(h.shape), (3, 3))
        for phase in range(3):
            row = h[phase]
            self.assertEqual(len(set(row.tolist())), 3)  # distinct components

    def test_dependency_radius_matches_spec(self):
        self.assertEqual(DEPENDENCY_RADIUS, 32 + 64 + 11 + 384)
        self.assertEqual(DEPENDENCY_RADIUS, 491)
        self.assertLess(DEPENDENCY_RADIUS, 516)  # below the proposed halo


if __name__ == "__main__":
    unittest.main()
