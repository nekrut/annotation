"""The learned pooled-decoder scalars of proposal 3.5 as decoder tables.

:class:`model.a.encoder.DecoderParams` holds 54 scalars: phase x component
mixture logits and hazard logits (the duration law of proposal 3.2), a
16-entry donor and a 16-entry acceptor dinucleotide score table, and four
partial-entry/exit family scores (section 3.6 boundary support: consumed by
the edge-enabled decoders, not yet by the loss). This module turns them into
what the kernels consume, keeping the gradient path:

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
  ``measure --decoder python``). It detaches: values only, no gradient.
* :func:`edge_prior` reads the four partial-family scalars as the
  :class:`model.grammar.EdgePrior` of proposal 3.1 (coding entry, coding
  exit, intron entry, intron exit, in the parameter's order) that the
  reference decoders and :func:`model.a.fast_viterbi.viterbi_batch_edges`
  take at an actual sequence edge. Values only, no gradient: the edge
  prior enters the training loss with the section-3.6 edge-partial
  numerators, not yet.
* :func:`as_torch_duration` gives the same law as a :class:`TorchDuration`,
  a duck-typed :class:`DurationMixture` whose ``log_pi``/``log_q``/``log_1mq``
  return 0-d tensors from :func:`duration_tables`. The reference recurrence
  (:mod:`model.a.torch_loss`) only ever *adds* those to emission tensors, so
  under it the reference chain loss carries the duration gradient too
  (engels-0083 P2: ``loss_kernel: reference`` used to fit with the 18
  duration scalars silently frozen).
"""
from __future__ import annotations

from typing import NamedTuple, Optional, Sequence

import torch
import torch.nn.functional as F

from model.grammar import DurationMixture, EdgePrior

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


class TorchDuration:
    """Differentiable stand-in for :class:`DurationMixture` (same ``m``, ``R``,
    ``log_pi``, ``log_q``, ``log_1mq``, ``log_prob`` surface) whose entries are
    0-d tensors sharing the graph of :func:`duration_tables`. For the Python
    reference decoders only; the fast kernels take the tables directly."""

    def __init__(self, tables: DurationTables, m: int = DEFAULT_M):
        if m < 1:
            raise ValueError("m must be at least 1")
        if tables.log_pi.shape != tables.log_q.shape or tables.log_pi.shape[0] != 3:
            raise ValueError("tables must be (3, R) with one row per phase")
        self.m = int(m)
        self.tables = tables

    @property
    def R(self) -> int:
        return int(self.tables.log_pi.shape[1])

    def log_pi(self, p, r):
        return self.tables.log_pi[p, r]

    def log_q(self, p, r):
        return self.tables.log_q[p, r]

    def log_1mq(self, p, r):
        return self.tables.log_1mq[p, r]

    def log_prob(self, p, length):
        if length < self.m:
            return torch.full((), float("-inf"), dtype=self.tables.log_pi.dtype,
                              device=self.tables.log_pi.device)
        t = self.tables
        return torch.logsumexp(t.log_pi[p] + t.log_1mq[p] + (length - self.m) * t.log_q[p], dim=0)


# ``DecoderParams.partial_families`` order: coding entry, coding exit, intron
# entry, intron exit (E0 / S, E / J / I, T of proposal 3.1).
PARTIAL_FAMILIES = ("entry", "exit", "intron_entry", "intron_exit")


def edge_prior(decoder) -> EdgePrior:
    """The decoder's four partial-family scalars as the sequence-edge prior
    of the reference and tensor decoders (values only, no gradient)."""
    fam = decoder.partial_families.detach().cpu()
    if tuple(fam.shape) != (len(PARTIAL_FAMILIES),):
        raise ValueError(f"partial_families must be ({len(PARTIAL_FAMILIES)},), got {tuple(fam.shape)}")
    return EdgePrior(**dict(zip(PARTIAL_FAMILIES, fam.tolist())))


def as_torch_duration(decoder, m: int = DEFAULT_M, *, dtype=torch.float64) -> TorchDuration:
    """The decoder's duration law with its gradient path, for the reference
    (torch) chain loss."""
    return TorchDuration(duration_tables(decoder, dtype=dtype), m)


def _dinuc_indices(x: str, n: int) -> "tuple[list, list]":
    """Per-position donor (``x[t:t+2]``) and acceptor (``x[t-2:t]``) row
    indices, ``_NO_SCORE`` where the dinucleotide is off the end or contains a
    non-ACGT character; pure-Python reference for :func:`_dinuc_index_tensors`."""
    up = x.upper()
    donor = [DINUC_INDEX.get(up[t:t + 2], _NO_SCORE) if t + 2 <= n else _NO_SCORE
             for t in range(n)]
    acceptor = [DINUC_INDEX.get(up[t - 2:t], _NO_SCORE) if t >= 2 else _NO_SCORE
                for t in range(n)]
    return donor, acceptor


def _dinuc_index_tensors(x: str) -> "tuple[torch.Tensor, torch.Tensor]":
    """Vectorised :func:`_dinuc_indices` over the base codes of ``x``
    (:func:`model.a.features.base_codes`); equal to the reference lists."""
    from .features import base_codes

    n = len(x)
    index, _ = base_codes(x)
    donor = torch.full((n,), _NO_SCORE, dtype=torch.long)
    acceptor = torch.full((n,), _NO_SCORE, dtype=torch.long)
    if n >= 2:
        first, second = index[:-1], index[1:]
        pair = torch.where((first < 4) & (second < 4), 4 * first + second, torch.full_like(first, _NO_SCORE))
        donor[:-1] = pair       # donor at t reads x[t:t+2]
        acceptor[2:] = pair[:-1] if n > 2 else pair[:0]  # acceptor at t reads x[t-2:t]
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
            donor_idx[b, :n], acceptor_idx[b, :n] = _dinuc_index_tensors(w)
    zero = torch.zeros(1, dtype=dtype, device=device)
    donor_tab = torch.cat((decoder.donor_dinuc.to(dtype=dtype, device=device), zero))
    acceptor_tab = torch.cat((decoder.acceptor_dinuc.to(dtype=dtype, device=device), zero))
    bias = torch.zeros((B, EMISSION_CHANNELS, L), dtype=dtype, device=device)
    bias[:, _DONOR, :] = donor_tab[donor_idx.to(device)]
    bias[:, _ACCEPTOR, :] = acceptor_tab[acceptor_idx.to(device)]
    return bias
