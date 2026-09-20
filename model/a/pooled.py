"""The learned pooled-decoder scalars of proposal 3.5 as decoder tables.

:class:`model.a.encoder.DecoderParams` holds 54 scalars: phase x component
mixture logits and hazard logits (the duration law of proposal 3.2), a
16-entry donor and a 16-entry acceptor dinucleotide score table, and four
partial-entry/exit family scores (section 3.6 boundary support, not consumed
yet). This module turns them into what the kernels consume, keeping the
gradient path:

* :func:`duration_tables` gives ``(log pi, log q, log(1 - q))`` as ``(3, R)``
  tensors -- ``pi[p] = softmax(mixture_logits[p])``, ``q = sigmoid(hazard)``
  -- which :func:`model.a.fast_loss.log_partition_batch` and
  :func:`model.a.fast_viterbi.viterbi_batch` take in place of the fixed
  :class:`model.grammar.DurationMixture` values. The grammar's *shape*
  (``m``, ``R``) stays a :class:`DurationMixture` (:func:`structure`), so the
  cached :class:`model.a.fast_loss.Grammar` tables are built once per code
  and not once per gradient step.
* :func:`motif_bias` adds the dinucleotide scores to the ``donor`` and
  ``acceptor`` emission channels: a donor at boundary ``t`` pays
  ``donor_dinuc[x[t] x[t+1]]`` (the first two intron bases), an acceptor at
  boundary ``t`` pays ``acceptor_dinuc[x[t-2] x[t-1]]`` (the last two).
  Dinucleotides that are not four concrete bases, or reach past the window,
  score 0. Adding them to the emissions before the scan keeps both kernels
  and the reference decoders unchanged.
* :func:`as_mixture` gives the same duration law as a concrete
  :class:`DurationMixture` for the Python reference decoders (parity runs,
  ``measure --decoder python``).
"""
from __future__ import annotations

from typing import NamedTuple, Optional, Sequence

import torch
import torch.nn.functional as F

from model.grammar import DurationMixture

from .torch_loss import CHANNEL_ORDER, EMISSION_CHANNELS

# Minimum intron length of the pilot grammar (proposal 3.2: m = 20).
DEFAULT_M = 20

DINUC_INDEX = {a + b: 4 * i + j for i, a in enumerate("ACGT") for j, b in enumerate("ACGT")}
_NO_SCORE = len(DINUC_INDEX)  # padded table row that scores 0
_DONOR, _ACCEPTOR = CHANNEL_ORDER.index("donor"), CHANNEL_ORDER.index("acceptor")


class DurationTables(NamedTuple):
    """``(3, R)`` log tables by phase: ``log pi``, ``log q``, ``log (1 - q)``."""
    log_pi: torch.Tensor
    log_q: torch.Tensor
    log_1mq: torch.Tensor

    def to(self, dtype=None, device=None) -> "DurationTables":
        return DurationTables(*(t.to(dtype=dtype, device=device) for t in self))


def duration_tables(decoder, *, dtype=torch.float64) -> DurationTables:
    """Differentiable duration tables from ``DecoderParams``."""
    mix = decoder.mixture_logits.to(dtype)
    haz = decoder.hazard_logits.to(dtype)
    return DurationTables(F.log_softmax(mix, dim=1), F.logsigmoid(haz), F.logsigmoid(-haz))


def structure(decoder, m: int = DEFAULT_M) -> DurationMixture:
    """The grammar shape ``(m, R)`` of ``decoder`` as a fixed
    :class:`DurationMixture`; only its ``m`` and ``R`` are read when the
    learned tables are supplied, so it is a stable cache key."""
    R = int(decoder.mixture_logits.shape[1])
    return DurationMixture(m=m, pi=((1.0 / R,) * R,) * 3, q=((0.5,) * R,) * 3)


def as_mixture(decoder, m: int = DEFAULT_M) -> DurationMixture:
    """The decoder's current duration law as a concrete
    :class:`DurationMixture` for the Python reference decoders."""
    t = duration_tables(decoder)
    pi = t.log_pi.exp().detach().cpu().tolist()
    q = t.log_q.exp().detach().cpu().tolist()
    pi = tuple(tuple(v / sum(row) for v in row) for row in pi)  # exact normalisation
    return DurationMixture(m=m, pi=pi, q=tuple(tuple(row) for row in q))


def _dinuc_indices(x: str, n: int) -> "tuple[list, list]":
    up = x.upper()
    donor = [DINUC_INDEX.get(up[t:t + 2], _NO_SCORE) if t + 2 <= n else _NO_SCORE
             for t in range(n)]
    acceptor = [DINUC_INDEX.get(up[t - 2:t], _NO_SCORE) if t >= 2 else _NO_SCORE
                for t in range(n)]
    return donor, acceptor


def motif_bias(x: str, decoder, *, dtype=torch.float64, device=None) -> torch.Tensor:
    """``(11, n)`` additive emission bias: the dinucleotide scores on the
    ``donor`` and ``acceptor`` rows, zero elsewhere."""
    return motif_bias_batch([x], decoder, dtype=dtype, device=device)[0, :, :len(x)]


def motif_bias_batch(windows: Sequence[str], decoder, *, L: Optional[int] = None,
                     dtype=torch.float64, device=None) -> torch.Tensor:
    """``(B, 11, L)`` :func:`motif_bias` for padded windows (zero past each
    window's end)."""
    B = len(windows)
    L = max((len(w) for w in windows), default=0) if L is None else L
    donor_idx = torch.full((B, L), _NO_SCORE, dtype=torch.long)
    acceptor_idx = torch.full((B, L), _NO_SCORE, dtype=torch.long)
    for b, w in enumerate(windows):
        n = len(w)
        if n > L:
            raise ValueError(f"window {b} is longer than L={L}")
        if n:
            d, a = _dinuc_indices(w, n)
            donor_idx[b, :n] = torch.tensor(d)
            acceptor_idx[b, :n] = torch.tensor(a)
    zero = torch.zeros(1, dtype=dtype, device=device)
    donor_tab = torch.cat((decoder.donor_dinuc.to(dtype=dtype, device=device), zero))
    acceptor_tab = torch.cat((decoder.acceptor_dinuc.to(dtype=dtype, device=device), zero))
    bias = torch.zeros((B, EMISSION_CHANNELS, L), dtype=dtype, device=device)
    bias[:, _DONOR, :] = donor_tab[donor_idx.to(device)]
    bias[:, _ACCEPTOR, :] = acceptor_tab[acceptor_idx.to(device)]
    return bias
