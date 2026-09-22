"""Candidate A encoder: compact DNA stem + pooled context (proposal section 3.5).

Maps an eight-channel DNA chunk (see `features`) to the eleven per-base
emission channels the `model.grammar` decoders consume. The exact layer
inventory and its 455,841-scalar total are fixed by proposal section 3.5;
`section_35_param_breakdown` reproduces that arithmetic from the formulas so
a test can check both the spec and the instantiated module without a GPU.

Shapes (B = batch, L = chunk length, a multiple of the pooling stride 12):
  input   x            [B, 8, L]
  stem                 [B, 16, L]
  pooled context       [B, 96, L/12]  -> repeated to [B, 96, L]
  fused                [B, 112, L]     (stem 16 ++ context 96)
  emissions            [B, 11, L]

The last heads are meaningful only on the retained central ("core")
positions; halo cropping for chunk-seam determinism is a loader concern
handled with the structured adapter, not here. The dependency radius is
491 bases (`DEPENDENCY_RADIUS`), below the proposed 516-base halo.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from .inventory import (
    ATTN_OFFSETS,
    CONTEXT_WIDTH,
    DEPENDENCY_RADIUS,
    EMISSION_CHANNELS,
    INPUT_CHANNELS,
    MLP_EXPANSION,
    N_ATTN_BLOCKS,
    N_HEADS,
    POOL_STRIDE,
    RESIDUAL_DILATIONS,
    SECTION_35_PARAM_COUNT,
    STEM_WIDTH,
    section_35_param_breakdown,
)

__all__ = [
    "CandidateA",
    "DNAEncoder",
    "DecoderParams",
    "SECTION_35_PARAM_COUNT",
    "section_35_param_breakdown",
    "DEPENDENCY_RADIUS",
]


class ChannelLayerNorm(nn.Module):
    """LayerNorm over the channel dimension of an [B, C, L] tensor.

    Section 3.5 normalizes over channels at each position, never over spatial
    positions, so the receptive-field argument holds.
    """

    def __init__(self, channels: int, eps: float = 1e-5):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(channels))
        self.bias = nn.Parameter(torch.zeros(channels))
        self.eps = eps

    def forward(self, x):  # x: [B, C, L]
        mean = x.mean(dim=1, keepdim=True)
        var = x.var(dim=1, keepdim=True, unbiased=False)
        x = (x - mean) / torch.sqrt(var + self.eps)
        return x * self.weight[None, :, None] + self.bias[None, :, None]


class ResidualBlock(nn.Module):
    """LN -> dilated depthwise conv (k=9) -> GELU -> pointwise conv, residual."""

    def __init__(self, width: int, dilation: int):
        super().__init__()
        self.norm = ChannelLayerNorm(width)
        pad = 4 * dilation  # keep length; kernel 9 has radius 4
        self.depthwise = nn.Conv1d(
            width, width, 9, padding=pad, dilation=dilation, groups=width
        )
        self.pointwise = nn.Conv1d(width, width, 1)

    def forward(self, x):
        h = self.norm(x)
        h = self.depthwise(h)
        h = F.gelu(h)
        h = self.pointwise(h)
        return x + h


class LocalAttention(nn.Module):
    """Multi-head attention restricted to relative offsets -8..+7.

    A single qkv projection and an output projection (4*d^2 + 4*d), plus a
    per-head learned bias over the 16 permitted offsets. The forward pass only
    forms scores and value products over the ``W = len(offsets)`` permitted
    offsets per query (via a gathered local window), so its attention work is
    O(T*W*hd) rather than the O(T^2*hd) of a masked dense product (engels-0059
    P2). ``_dense_forward`` keeps the equivalent dense implementation as an
    oracle the tests check against; both share the same weights and biases.
    """

    def __init__(self, dim: int, n_heads: int, offsets):
        super().__init__()
        assert dim % n_heads == 0
        self.dim = dim
        self.n_heads = n_heads
        self.head_dim = dim // n_heads
        self.offsets = tuple(offsets)
        self.lo, self.hi = self.offsets[0], self.offsets[-1]
        # Contiguous offsets let the window index equal (offset - lo), so a
        # single unfold gathers the per-query neighbourhood without a scatter.
        assert self.offsets == tuple(range(self.lo, self.hi + 1))
        self.qkv = nn.Linear(dim, 3 * dim)
        self.out = nn.Linear(dim, dim)
        # rel_bias[h, k] is the bias for offset offsets[k], per head.
        self.rel_bias = nn.Parameter(torch.zeros(n_heads, len(self.offsets)))

    def _window_mask(self, T, device):
        # valid[i, w] is True when key i + offsets[w] is a real position.
        i = torch.arange(T, device=device)[:, None]
        w = torch.arange(len(self.offsets), device=device)[None, :]
        key = i + w + self.lo
        return (key >= 0) & (key < T)  # [T, W]

    def forward(self, x):  # x: [B, T, C]
        B, T, C = x.shape
        W = len(self.offsets)
        qkv = self.qkv(x).reshape(B, T, 3, self.n_heads, self.head_dim)
        q, k, v = qkv.permute(2, 0, 3, 1, 4)  # each [B, H, T, hd]
        left, right = max(0, -self.lo), max(0, self.hi)
        # Gather, for every query, only the W keys/values at permitted offsets.
        # Pad along T so window w of query i reads padded index i+w == i+offset.
        kpad = F.pad(k, (0, 0, left, right))  # [B, H, T+W-1, hd]
        vpad = F.pad(v, (0, 0, left, right))
        kw = kpad.unfold(2, W, 1).permute(0, 1, 2, 4, 3)  # [B, H, T, W, hd]
        vw = vpad.unfold(2, W, 1).permute(0, 1, 2, 4, 3)
        scores = (q.unsqueeze(3) * kw).sum(-1) / (self.head_dim ** 0.5)  # [B,H,T,W]
        scores = scores + self.rel_bias[None, :, None, :]
        valid = self._window_mask(T, x.device)          # [T, W]
        scores = scores.masked_fill(~valid[None, None], float("-inf"))
        attn = torch.softmax(scores, dim=-1)             # [B, H, T, W]
        ctx = (attn.unsqueeze(-1) * vw).sum(3)           # [B, H, T, hd]
        ctx = ctx.transpose(1, 2).reshape(B, T, C)
        return self.out(ctx)

    def _dense_forward(self, x):  # oracle: identical output, O(T^2) work
        B, T, C = x.shape
        qkv = self.qkv(x).reshape(B, T, 3, self.n_heads, self.head_dim)
        q, k, v = qkv.permute(2, 0, 3, 1, 4)
        scores = torch.matmul(q, k.transpose(-1, -2)) / (self.head_dim ** 0.5)
        i = torch.arange(T, device=x.device)[:, None]
        j = torch.arange(T, device=x.device)[None, :]
        off = j - i
        in_window = (off >= self.lo) & (off <= self.hi)
        idx = (off - self.lo).clamp(0, len(self.offsets) - 1)
        scores = scores + self.rel_bias[:, idx][None]
        scores = scores.masked_fill(~in_window[None, None], float("-inf"))
        attn = torch.softmax(scores, dim=-1)
        ctx = torch.matmul(attn, v)
        ctx = ctx.transpose(1, 2).reshape(B, T, C)
        return self.out(ctx)


class ContextBlock(nn.Module):
    """Pre-norm local-attention + pre-norm expansion-4 MLP (section 3.5)."""

    def __init__(self, dim: int, n_heads: int, offsets, expansion: int):
        super().__init__()
        self.norm1 = nn.LayerNorm(dim)
        self.attn = LocalAttention(dim, n_heads, offsets)
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, expansion * dim),
            nn.GELU(),
            nn.Linear(expansion * dim, dim),
        )

    def forward(self, x):  # x: [B, T, C]
        x = x + self.attn(self.norm1(x))
        x = x + self.mlp(self.norm2(x))
        return x


class DNAEncoder(nn.Module):
    """Section 3.5 encoder: 8-channel input -> 11 emission channels per base."""

    def __init__(self):
        super().__init__()
        self.stem_conv = nn.Conv1d(INPUT_CHANNELS, STEM_WIDTH, 9, padding=4)
        self.residual = nn.ModuleList(
            ResidualBlock(STEM_WIDTH, d) for d in RESIDUAL_DILATIONS
        )
        self.pool_proj = nn.Conv1d(STEM_WIDTH, CONTEXT_WIDTH, 1)
        self.context = nn.ModuleList(
            ContextBlock(CONTEXT_WIDTH, N_HEADS, ATTN_OFFSETS, MLP_EXPANSION)
            for _ in range(N_ATTN_BLOCKS)
        )
        self.fuse = nn.Conv1d(STEM_WIDTH + CONTEXT_WIDTH, 32, 1)
        self.emission = nn.Conv1d(32, EMISSION_CHANNELS, 1)

    def forward(self, x):  # x: [B, 8, L], L % POOL_STRIDE == 0
        B, _, L = x.shape
        if L % POOL_STRIDE != 0:
            raise ValueError("chunk length must be a multiple of %d" % POOL_STRIDE)
        stem = self.stem_conv(x)
        for block in self.residual:
            stem = block(stem)
        # Pooled context path, downsampled by 12.
        pooled = F.avg_pool1d(stem, POOL_STRIDE)      # [B, 16, L/12]
        ctx = self.pool_proj(pooled)                  # [B, 96, L/12]
        ctx = ctx.transpose(1, 2)                     # [B, T, 96]
        for block in self.context:
            ctx = block(ctx)
        ctx = ctx.transpose(1, 2)                     # [B, 96, T]
        ctx = ctx.repeat_interleave(POOL_STRIDE, dim=2)  # [B, 96, L]
        fused = torch.cat([stem, ctx], dim=1)         # [B, 112, L]
        fused = F.gelu(self.fuse(fused))              # [B, 32, L]
        return self.emission(fused)                   # [B, 11, L]


class DecoderParams(nn.Module):
    """The 54 pooled-decoder scalars of section 3.5.

    Phase/component mixture and hazard logits, donor/acceptor dinucleotide
    score tables, and four partial-entry/exit family scores. Start/stop
    compatibility and prefix transitions are fixed by the genetic code, not
    learned. These feed the `model.grammar` decoders; the recurrence itself
    lives there.
    """

    # Distinct per-component hazard logits recorded for the run manifest.
    # Equal mixture weights are fine, but identical hazards make the three
    # geometric components coincide: under the softmax/sigmoid duration law the
    # mixture-logit gradients are then zero and the hazard-logit gradients are
    # equal across components (they can be nonzero but move the components
    # together), so the components cannot differentiate under symmetric updates
    # (the mixture stays a single geometric regardless of component count; see
    # stalin-0062, engels-0060). Seeding distinct hazards per component,
    # broadcast across phases, breaks that symmetry deterministically without
    # adding parameters.
    HAZARD_LOGIT_INIT = (-1.0, 0.0, 1.0)  # per component, same for all phases

    def __init__(self):
        super().__init__()
        self.mixture_logits = nn.Parameter(torch.zeros(3, 3))  # phase x component
        hazard_init = torch.tensor(self.HAZARD_LOGIT_INIT, dtype=torch.float32)
        self.hazard_logits = nn.Parameter(
            hazard_init.unsqueeze(0).repeat(3, 1)  # [phase, component], distinct
        )
        self.donor_dinuc = nn.Parameter(torch.zeros(16))
        self.acceptor_dinuc = nn.Parameter(torch.zeros(16))
        self.partial_families = nn.Parameter(torch.zeros(4))   # coding/intron x entry/exit


class CandidateA(nn.Module):
    """Full candidate A parameter set: encoder + pooled decoder = 455,841."""

    def __init__(self):
        super().__init__()
        self.encoder = DNAEncoder()
        self.decoder = DecoderParams()

    def forward(self, x):
        return self.encoder(x)

    def num_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters())
