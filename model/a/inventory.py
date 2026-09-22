"""Candidate A parameter inventory (proposal section 3.5), pure arithmetic.

No torch import, so the spec formula and its 455,841-scalar total are usable
and testable without a GPU. `model.a.encoder` reuses these constants to build
the module, and a test checks that the built module matches this total.
"""
from __future__ import annotations

STEM_WIDTH = 16
CONTEXT_WIDTH = 96
POOL_STRIDE = 12
N_RESIDUAL_BLOCKS = 3
RESIDUAL_DILATIONS = (1, 2, 4)
N_ATTN_BLOCKS = 4
N_HEADS = 4
ATTN_OFFSETS = tuple(range(-8, 8))  # -8 .. +7, i.e. 16 relative positions
MLP_EXPANSION = 4
EMISSION_CHANNELS = 11
INPUT_CHANNELS = 8
DEPENDENCY_RADIUS = 491  # 4 + 4*(1+2+4) stem, +64 GC, +11 pool/repeat, +384 ctx


def section_35_param_breakdown() -> dict:
    """Trainable-scalar counts per proposal section 3.5, from the formulas."""
    d = CONTEXT_WIDTH
    stem_conv = INPUT_CHANNELS * STEM_WIDTH * 9 + STEM_WIDTH
    residual = N_RESIDUAL_BLOCKS * (
        STEM_WIDTH * 9 + STEM_WIDTH             # depthwise conv + bias
        + STEM_WIDTH * STEM_WIDTH + STEM_WIDTH  # pointwise conv + bias
        + 2 * STEM_WIDTH                        # one layer norm (weight + bias)
    )
    pooled_proj = STEM_WIDTH * d + d
    attn_blocks = N_ATTN_BLOCKS * (
        12 * d * d                     # 4*d^2 attention (qkv+out) + 8*d^2 MLP
        + 9 * d                        # attention + MLP biases
        + 4 * d                        # two layer norms (weight + bias each)
        + N_HEADS * len(ATTN_OFFSETS)  # learned relative-offset bias per head
    )
    fusion = (STEM_WIDTH + CONTEXT_WIDTH) * 32 + 32
    emission = 32 * EMISSION_CHANNELS + EMISSION_CHANNELS
    pooled_decoder = 3 * 3 * 2 + 2 * STEM_WIDTH + 4  # mixtures/hazards, dinuc, partials
    breakdown = {
        "stem_convolution": stem_conv,
        "residual_blocks": residual,
        "pooled_context_projection": pooled_proj,
        "attention_blocks": attn_blocks,
        "fine_context_fusion": fusion,
        "emission_projection": emission,
        "pooled_decoder": pooled_decoder,
    }
    breakdown["total"] = sum(breakdown.values())
    return breakdown


SECTION_35_PARAM_COUNT = section_35_param_breakdown()["total"]  # == 455_841
