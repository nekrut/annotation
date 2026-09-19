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

    def test_dependency_radius_matches_spec(self):
        self.assertEqual(DEPENDENCY_RADIUS, 32 + 64 + 11 + 384)
        self.assertEqual(DEPENDENCY_RADIUS, 491)
        self.assertLess(DEPENDENCY_RADIUS, 516)  # below the proposed halo


if __name__ == "__main__":
    unittest.main()
