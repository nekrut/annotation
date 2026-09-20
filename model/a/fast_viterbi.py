"""Tensor Viterbi for candidate A's grammar: the max-product twin of
:func:`model.a.fast_loss.log_partition_batch` with back-pointers and a
traceback that returns the same :class:`model.grammar.reference.Chain`
objects as :meth:`model.grammar.delayed.DelayedEntryDecoder.viterbi`.

The pure-Python delayed-entry decoder costs ~66 CPU-s/Mb on yeast chr I
(a-pilot.md section 3.2), the failing stage of the CPU regime. This scan
replaces it: the same ``n``-step loop over ``(B, K)`` and ``(B, K, R)``
tensors as the training kernel, ``max``/``argmax`` in place of
``logsumexp``, sharing the grammar tables of :mod:`model.a.fast_loss` but
building the coding transition per step from its sparse parts
(:func:`_operands`) instead of the dense ``(B, L, K, K)`` stack.

Back-pointers per boundary ``t + 1`` (``t`` is the consumed base):

* ``prev[t + 1, j]`` (int8): the coding state ``i`` at boundary ``t`` whose
  transition reached coding state ``j``;
* ``exit_r[t, i]`` (int8): ``-1`` if coding state ``i`` at boundary ``t``
  scored from the coding layer itself, else the component ``r`` of the tail
  ``T(i, r)`` whose exit at ``t`` beat it;
* ``entered[t + 1, c, r]`` (bool): tail ``T(c, r)`` at ``t + 1`` was entered
  by the donor parked at boundary ``t + 1 - m`` (its ``m - 1`` mandatory
  intronic positions are re-expanded to ``I(c, k)`` states on traceback,
  exactly as the Python decoder does), else it continued.

Storage is ``L * K * (2 + R)`` bytes per window, so a whole chromosome would
want the checkpoint-and-replay traceback of proposal 3.3 (not offered here);
``measure`` decodes annotation-selected windows and does not need it.

Complete-target only, like the training kernel: the path must end in ``U``
at boundary ``n``. Ties are broken by the first maximal index, which can
differ from the Python decoder's insertion order; the score is tie-free.
"""
from __future__ import annotations

from itertools import groupby
from typing import List, Optional, Sequence, Tuple

import torch

from model.grammar import DurationMixture, GeneticCode, TABLES
from model.grammar.reference import Chain, ReferenceDecoder

from .fast_loss import _CH, FLOOR, SYMBOL_SETS, Grammar, symbol_indices
from .torch_loss import EMISSION_CHANNELS, _check_input


def _operands(emissions: torch.Tensor, symbols: torch.Tensor, lengths: torch.Tensor,
              g: Grammar):
    """Per-step operands of the max-product scan, without materializing the
    dense ``(B, L, K, K)`` transition stack of the training kernel (at 16
    windows of 12 kb that stack is 0.9 GB and ~80% of the decode time).

    The coding transition ``L[t][i, j]`` decomposes as ``base[t][i]`` (the
    ``cds[p(i)]`` emission of the predecessor; 0 for ``U``) plus the
    symbol-conditional prior ``prior_max[x[t]][i, j]`` plus three special
    entries: ``u[t]`` on ``U -> U``, ``stop[t]`` on ``coding -> U`` (both
    folded into ``ucol[t][i]``, the column-``U`` extra by predecessor) and
    ``start[t] + cds[0][t]`` on ``U -> S(b)`` (``post[t][j]``, added after the
    max since only ``U`` reaches a length-1 initiator prefix)."""
    B, C, L = emissions.shape
    K, m = g.K, g.dur.m
    dtype, device = emissions.dtype, emissions.device
    e = emissions.clamp(min=FLOOR)
    valid = (torch.arange(L, device=device)[None, :] < lengths[:, None])      # (B, L)
    zero = torch.zeros((B, L), dtype=dtype, device=device)
    neg = torch.full((B, L), FLOOR, dtype=dtype, device=device)

    def masked(ch, fill):
        return torch.where(valid, e[:, _CH[ch], :], fill)

    sym_pad = torch.where(valid, symbols, torch.full_like(symbols, len(SYMBOL_SETS)))
    cds = torch.stack([masked(f"cds[{p}]", zero) for p in range(3)], 1)      # (B, 3, L)
    base = cds[:, g.phase, :]                                                 # (B, K, L)
    base[:, g.u_index, :] = 0.0
    ucol = masked("stop", zero)[:, None, :].expand(B, K, L).clone()          # (B, K, L)
    ucol[:, g.u_index, :] = masked("u", zero)
    post = torch.zeros((B, K, L), dtype=dtype, device=device)
    s1 = torch.tensor([s[0] == "S" and len(s[1]) == 1 for s in g.states], device=device)
    post[:, s1, :] = (masked("start", zero) + masked("cds[0]", zero))[:, None, :]
    donor = masked("donor", neg)[:, None, :].expand(B, K, L).clone()
    donor[:, ~g.coding, :] = FLOOR
    acceptor = masked("acceptor", neg)
    intron = torch.where(valid[:, None, :], e[:, 4:7, :], neg[:, None, :])    # (B, 3, L)
    if L >= m:
        window_sum = intron.unfold(2, m, 1).sum(-1)                           # (B, 3, L-m+1)
    else:
        window_sum = intron.new_full((B, 3, 0), FLOOR)
    intron_q = intron[:, g.phase, :, None] + g.log_q[None, :, None, :]        # (B, K, L, R)
    window_pi = window_sum[:, g.phase, :, None] + g.log_pi[None, :, None, :]  # (B, K, L-m+1, R)
    return (sym_pad.unbind(1), base.unbind(2), ucol.unbind(2), post.unbind(2),
            donor.unbind(2), acceptor.unbind(1), intron_q.unbind(2), window_pi.unbind(2))


def viterbi_batch(emissions: torch.Tensor, symbols: torch.Tensor,
                  lengths: torch.Tensor, grammar: Grammar):
    """Max-product scan. Returns ``(score, prev, exit_r, entered)``: ``score``
    is ``(B,)`` with ``-inf`` where no complete path exists; the back-pointer
    tensors are ``(B, L + 1, K)`` int8, ``(B, L + 1, K)`` int8 and
    ``(B, L + 1, K, R)`` bool (row 0 unused). ``exit_r[t, i]`` is indexed by
    the *predecessor* ``i`` at boundary ``t``: the component of the tail
    ``T(i, r)`` whose exit gave ``i`` its merged score at ``t``, or ``-1`` if
    the coding layer did."""
    B, C, L = emissions.shape
    if C != EMISSION_CHANNELS:
        raise ValueError(f"emissions must be (B, {EMISSION_CHANNELS}, L)")
    g = grammar
    K, R, m, u = g.K, g.R, g.dur.m, g.u_index
    dtype, device = emissions.dtype, emissions.device
    sym_t, base_t, ucol_t, post_t, donor_t, acceptor_t, intron_q, window_pi = _operands(
        emissions, symbols, lengths, g)
    prior = g.prior_max
    log_1mq = g.log_1mq[None]
    neg_one = torch.full((B, K), -1, dtype=torch.int8, device=device)

    prev = torch.zeros((B, L + 1, K), dtype=torch.int8, device=device)
    exit_r = torch.full((B, L + 1, K), -1, dtype=torch.int8, device=device)
    entered = torch.zeros((B, L + 1, K, R), dtype=torch.bool, device=device)

    alpha = torch.full((B, K), FLOOR, dtype=dtype, device=device)
    alpha[:, u] = 0.0
    tau = torch.full((B, K, R), FLOOR, dtype=dtype, device=device)
    pending: List[torch.Tensor] = []
    with torch.no_grad():
        for t in range(L):
            pending.append(alpha + donor_t[t])
            exits, best_r = (tau + log_1mq).max(dim=-1)                       # (B, K)
            exits += acceptor_t[t][:, None]
            from_tail = exits > alpha
            exit_r[:, t] = torch.where(from_tail, best_r.to(torch.int8), neg_one)
            merged = torch.where(from_tail, exits, alpha) + base_t[t]
            trans = merged[:, :, None] + prior[sym_t[t]]                      # (B, K, K)
            trans[:, :, u] += ucol_t[t]
            alpha, best_i = trans.max(dim=1)
            alpha += post_t[t]
            alpha.clamp_(min=FLOOR)
            prev[:, t + 1] = best_i.to(torch.int8)
            tau += intron_q[t]
            if len(pending) == m:
                entry = pending.pop(0)[:, :, None] + window_pi[t - m + 1]
                enter = entry > tau
                entered[:, t + 1] = enter
                tau = torch.where(enter, entry, tau)
            tau.clamp_(min=FLOOR)
    z = alpha[:, u]
    score = torch.where(z > FLOOR / 2, z, torch.full_like(z, float("-inf")))
    return score, prev, exit_r, entered


def traceback(n: int, grammar: Grammar, prev: torch.Tensor, exit_r: torch.Tensor,
              entered: torch.Tensor) -> List[tuple]:
    """State tuples at boundaries ``0..n`` of one window from its back-pointers
    (each ``(L + 1, ...)``), in the reference decoder's vocabulary:
    ``("U",)``, ``("S", q)``, ``("E", q)``, ``("I", c, k)``, ``("T", c, r)``."""
    g = grammar
    m = g.dur.m
    # numpy scalar reads are ~10x cheaper than tensor indexing or ``tolist``.
    prev, exit_r, entered = (a.cpu().numpy() for a in (prev, exit_r, entered))
    states: List[Optional[tuple]] = [None] * (n + 1)
    t, j, r = n, g.u_index, -1                        # r >= 0: in tail T(j, r)
    while t > 0:
        if r < 0:
            states[t] = g.states[j]
            i = int(prev[t, j])
            t -= 1
            j, r = i, int(exit_r[t, i])
        else:
            states[t] = ("T", g.states[j], r)
            if entered[t, j, r]:
                s = t - m
                for k in range(s + 1, t):
                    states[k] = ("I", g.states[j], k - s)
                t, r = s, -1
            else:
                t -= 1
    states[0] = g.states[j] if r < 0 else ("T", g.states[j], r)
    return states


def viterbi(x: str, emissions: torch.Tensor, *, code: GeneticCode = TABLES[1],
            duration: DurationMixture = DurationMixture()) -> Tuple[float, List[Chain]]:
    """Best complete-path score and its chains for one window; the tensor twin
    of ``DelayedEntryDecoder(code, duration).viterbi(x, scores)`` without edge
    priors. Returns ``(-inf, [])`` when the grammar admits no path."""
    _check_input(x, emissions)
    g = Grammar.get(code, duration, dtype=emissions.dtype, device=emissions.device)
    sym = torch.tensor([symbol_indices(x)], dtype=torch.long, device=emissions.device)
    lengths = torch.tensor([len(x)], device=emissions.device)
    score, prev, exit_r, entered = viterbi_batch(emissions[None], sym, lengths, g)
    best = float(score[0])
    if best == float("-inf"):
        return best, []
    states = traceback(len(x), g, prev[0], exit_r[0], entered[0])
    return best, ReferenceDecoder._chains(states, [None] * len(x), x)


def viterbi_windows(windows: Sequence[str], emissions: Sequence[torch.Tensor], *,
                    codes: Optional[Sequence[GeneticCode]] = None,
                    duration: DurationMixture = DurationMixture(),
                    batch_size: int = 16) -> List[Tuple[float, List[Chain]]]:
    """:func:`viterbi` over many windows, batched to amortize the per-step
    launch cost of the scan (on one CPU core the scan over a ``(B, K)`` layer
    costs about the same as over ``(1, K)``, so decoding ``B`` windows
    together is nearly ``B`` times cheaper per base). Windows are grouped by
    genetic code and sorted by length so padding stays small; each batch is
    padded to its longest window with ``-inf`` emissions past the end, which
    :func:`model.a.fast_loss.prepare_batch` turns into identity transitions.
    ``emissions[i]`` is ``(11, len(windows[i]))``. Results are returned in the
    input order."""
    if len(windows) != len(emissions):
        raise ValueError("windows and emissions must agree in length")
    if codes is None:
        codes = [TABLES[1]] * len(windows)
    if len(codes) != len(windows):
        raise ValueError("codes must have one entry per window")
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    out: List[Optional[Tuple[float, List[Chain]]]] = [None] * len(windows)
    order = sorted(range(len(windows)), key=lambda i: (codes[i].table, len(windows[i])))
    for code, group in groupby(order, key=lambda i: codes[i]):
        idx = list(group)
        for start in range(0, len(idx), batch_size):
            ids = idx[start:start + batch_size]
            dtype, device = emissions[ids[0]].dtype, emissions[ids[0]].device
            g = Grammar.get(code, duration, dtype=dtype, device=device)
            lengths = [len(windows[i]) for i in ids]
            L = max(lengths)
            e = torch.full((len(ids), emissions[ids[0]].shape[0], L), float("-inf"),
                           dtype=dtype, device=device)
            sym = torch.zeros((len(ids), L), dtype=torch.long, device=device)
            for b, i in enumerate(ids):
                _check_input(windows[i], emissions[i])
                e[b, :, :lengths[b]] = emissions[i]
                if lengths[b]:
                    sym[b, :lengths[b]] = torch.tensor(symbol_indices(windows[i]), device=device)
            score, prev, exit_r, entered = viterbi_batch(
                e, sym, torch.tensor(lengths, device=device), g)
            for b, i in enumerate(ids):
                best = float(score[b])
                if best == float("-inf"):
                    out[i] = (best, [])
                    continue
                states = traceback(lengths[b], g, prev[b], exit_r[b], entered[b])
                out[i] = (best, ReferenceDecoder._chains(states, [None] * lengths[b], windows[i]))
    return out  # type: ignore[return-value]
