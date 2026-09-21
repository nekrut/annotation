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
(:func:`_operands`) instead of the dense ``(B, L, K, K)`` stack, and taking
the per-step maximum over each state's *predecessors* only
(:func:`_sparse`): every coding state has exactly one except ``U`` (itself
and the stop-completing two-base prefixes) and the empty prefix ``E("")``
(every two-base prefix and the initiator-completing ``S`` prefixes), so the
step gathers ``K - 2 + 2 P`` candidates per window instead of the dense
``K x K`` block. :func:`viterbi_batch_reference` keeps the dense block; the
two agree bit for bit on scores and on every back-pointer a traceback can
read (``tests/test_a_fast_viterbi.py``).

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
Back-pointers of states no finite path reaches, and of boundaries past a
window's end inside a padded batch, are unspecified (the dense reference
fills them with the identity; the sparse scan does not).
"""
from __future__ import annotations

from itertools import groupby
from typing import List, NamedTuple, Optional, Sequence, Tuple

import torch

from model.grammar import DurationMixture, GeneticCode, TABLES
from model.grammar.reference import Chain, ReferenceDecoder

from .fast_loss import (_CH, FLOOR, SYMBOL_SETS, Grammar, duration_by_phase, duration_by_state,
                        symbol_index_tensor, symbol_indices)
from .torch_loss import EMISSION_CHANNELS, _check_input


def _operands(emissions: torch.Tensor, symbols: torch.Tensor, lengths: torch.Tensor,
              g: Grammar, tables=None):
    """Per-step operands of the max-product scan, without materializing the
    dense ``(B, L, K, K)`` transition stack of the training kernel (at 16
    windows of 12 kb that stack is 0.9 GB and ~80% of the decode time).

    The coding transition ``L[t][i, j]`` decomposes as ``base[t][i]`` (the
    ``cds[p(i)]`` emission of the predecessor; 0 for ``U``) plus the
    symbol-conditional prior ``prior_max[x[t]][i, j]`` plus three special
    entries: ``u[t]`` on ``U -> U``, ``stop[t]`` on ``coding -> U`` (both
    folded into ``ucol[t][i]``, the column-``U`` extra by predecessor) and
    ``start[t] + cds[0][t]`` on ``U -> S(b)`` (``post[t][j]``, added after the
    max since only ``U`` reaches a length-1 initiator prefix). ``ucol`` is
    returned at ``U``'s predecessor candidates only, ``(B, P)`` per step, the
    one column the sparse scan adds it to. Each operand is a tuple of ``L``
    contiguous per-step tensors."""
    B, C, L = emissions.shape
    K, m = g.K, g.dur.m
    sp = _sparse(g)
    log_pi, log_q, _ = duration_by_phase(g, tables)
    dtype, device = emissions.dtype, emissions.device
    e = emissions.clamp(min=FLOOR)
    valid = (torch.arange(L, device=device)[None, :] < lengths[:, None])      # (B, L)
    zero = torch.zeros((B, L), dtype=dtype, device=device)
    neg = torch.full((B, L), FLOOR, dtype=dtype, device=device)

    # Every operand is built step-major, ``(L, B, ...)``, so each step's slice
    # is contiguous: a slice of ``(B, K, L)`` along ``L`` is strided, and
    # ``index_select`` on such a slice (or with such an index) clones it
    # every step; transposing the ``K``-wide stacks afterwards costs more
    # than the scan saves, so only the ``(B, L)`` channels are transposed.
    def masked_t(ch, fill):
        return torch.where(valid, e[:, _CH[ch], :], fill).T                   # (L, B)

    sym_pad = torch.where(valid, symbols, torch.full_like(symbols, len(SYMBOL_SETS))).T
    cds = torch.stack([masked_t(f"cds[{p}]", zero) for p in range(3)], 2)    # (L, B, 3)
    base = cds.index_select(2, g.phase)                                       # (L, B, K)
    base[:, :, g.u_index] = 0.0
    # Column-``U`` extras at ``U``'s predecessor candidates: ``u`` on ``U``
    # (also the padding slots), ``stop`` on the stop-completing prefixes.
    ucol = torch.where((sp.cand[0] == g.u_index)[None, None, :],
                       masked_t("u", zero)[:, :, None], masked_t("stop", zero)[:, :, None])
    post = torch.zeros((L, B, K), dtype=dtype, device=device)
    s1 = torch.tensor([s[0] == "S" and len(s[1]) == 1 for s in g.states], device=device)
    post[:, :, s1] = (masked_t("start", zero) + masked_t("cds[0]", zero))[:, :, None]
    donor = masked_t("donor", neg)[:, :, None].expand(L, B, K).clone()
    donor[:, :, ~g.coding] = FLOOR
    acceptor = masked_t("acceptor", neg)
    intron = torch.where(valid[:, None, :], e[:, 4:7, :], neg[:, None, :])    # (B, 3, L)
    if L >= m:
        window_sum = intron.unfold(2, m, 1).sum(-1)                           # (B, 3, L-m+1)
    else:
        window_sum = intron.new_full((B, 3, 0), FLOOR)
    # Per phase, not per state: the scan expands ``(B, 3, R)`` to ``(B, K, R)``
    # with one ``index_select`` per step, which costs less than building and
    # holding the ``K / 3``-fold larger ``(B, K, L, R)`` stacks (2 x 113 MB at
    # 16 windows of 12 kb, 2 x 450 MB at 64).
    intron_q = intron.permute(2, 0, 1).contiguous()[..., None] + log_q            # (L, B, 3, R)
    window_pi = window_sum.permute(2, 0, 1).contiguous()[..., None] + log_pi      # (L-m+1, B, 3, R)
    return tuple(a.contiguous().unbind(0)
                 for a in (sym_pad, base, ucol, post, donor, acceptor, intron_q, window_pi))



class _Sparse(NamedTuple):
    gather: torch.Tensor    # (K - 2 + 2 P,) predecessor index per candidate slot
    prior: torch.Tensor     # (S + 1, K - 2 + 2 P) prior per symbol row and slot
    cand: torch.Tensor      # (2, P) candidate states of U (row 0) and E("") (row 1)
    e0: int                 # index of E(""); U is 0, single-predecessor states fill the rest


def _sparse(g: Grammar) -> _Sparse:
    """Predecessor tables of the sparse scan for one grammar, cached on it.
    Slot order is the single-predecessor states in index order (``1 .. e0 - 1``,
    ``e0 + 1 .. K - 1``), then ``U``'s candidates, then ``E("")``'s, each
    padded to ``P`` with ``U`` under a ``FLOOR`` prior; candidates ascend in
    state index so the first maximum agrees with the dense scan's."""
    tab = getattr(g, "_viterbi_sparse", None)
    if tab is not None:
        return tab
    K, u, S = g.K, g.u_index, len(SYMBOL_SETS)
    e0 = g.index[("E", "")]
    dtype, device = g.prior_max.dtype, g.prior_max.device
    permitted = (g.prior_max[:S] > FLOOR / 2).any(0)                           # (K, K) i -> j
    preds = [permitted[:, j].nonzero().flatten() for j in range(K)]
    single = list(range(1, e0)) + list(range(e0 + 1, K))
    if u != 0 or any(len(preds[j]) != 1 for j in single):
        raise AssertionError("grammar states are not U, single-predecessor prefixes, E('')")
    one_pred = torch.tensor([int(preds[j][0]) for j in single], device=device)
    P = max(len(preds[u]), len(preds[e0]))
    cand = torch.full((2, P), u, dtype=torch.long, device=device)
    prior = torch.full((S + 1, 2, P), FLOOR, dtype=dtype, device=device)
    for row, j in enumerate((u, e0)):
        cand[row, :len(preds[j])] = preds[j]
        prior[:, row, :len(preds[j])] = g.prior_max[:, preds[j], j]
    tab = _Sparse(gather=torch.cat([one_pred, cand.flatten()]),
                  prior=torch.cat([g.prior_max[:, one_pred, torch.tensor(single, device=device)],
                                   prior.view(S + 1, 2 * P)], 1).contiguous(),
                  cand=cand, e0=e0)
    g._viterbi_sparse = tab
    return tab


def viterbi_batch(emissions: torch.Tensor, symbols: torch.Tensor,
                  lengths: torch.Tensor, grammar: Grammar, tables=None):
    """Max-product scan. Returns ``(score, prev, exit_r, entered)``: ``score``
    is ``(B,)`` with ``-inf`` where no complete path exists; the back-pointer
    tensors are ``(B, L + 1, K)`` int8, ``(B, L + 1, K)`` int8 and
    ``(B, L + 1, K, R)`` bool (row 0 unused). ``exit_r[t, i]`` is indexed by
    the *predecessor* ``i`` at boundary ``t``: the component of the tail
    ``T(i, r)`` whose exit gave ``i`` its merged score at ``t``, or ``-1`` if
    the coding layer did. ``tables`` (learned
    :class:`model.a.pooled.DurationTables`) replaces the grammar's fixed
    duration values.

    The coding transition is the sparse-predecessor form (:func:`_sparse`):
    one ``index_select`` of the merged layer at every slot, the symbol's
    prior row added, ``U``'s column extras on its candidates, one ``max``
    over the two ``(B, P)`` candidate rows, and the ``(B, K)`` layer
    reassembled by ``cat``. Past a window's end the pad symbol row keeps
    ``U -> U`` at 0 and everything else at ``FLOOR``, so ``alpha[U]`` is
    preserved to boundary ``L`` while the other states are unspecified."""
    B, C, L = emissions.shape
    if C != EMISSION_CHANNELS:
        raise ValueError(f"emissions must be (B, {EMISSION_CHANNELS}, L)")
    g = grammar
    K, R, m, u = g.K, g.R, g.dur.m, g.u_index
    dtype, device = emissions.dtype, emissions.device
    sym_t, base_t, ucol_t, post_t, donor_t, acceptor_t, intron_q, window_pi = _operands(
        emissions, symbols, lengths, g, tables)
    sp = _sparse(g)
    gather, prior, e0 = sp.gather, sp.prior, sp.e0
    K1, P = K - 2, sp.cand.shape[1]
    log_1mq = duration_by_state(g, tables)[2][None]
    phase = g.phase
    neg_one = torch.full((B, K), -1, dtype=torch.long, device=device)

    prev2 = torch.zeros((B, L + 1, 2), dtype=torch.int8, device=device)
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
            # ``copy_`` into the int8 stores converts; no separate ``to``.
            exit_r[:, t] = torch.where(from_tail, best_r, neg_one)
            merged = torch.where(from_tail, exits, alpha) + base_t[t]
            v = merged.index_select(1, gather) + prior.index_select(0, sym_t[t])  # (B, K1 + 2P)
            cand = v[:, K1:].view(B, 2, P)
            cand[:, 0] += ucol_t[t]
            best, best_p = cand.max(dim=-1)                                   # (B, 2)
            alpha = torch.cat([best[:, :1], v[:, :e0 - 1], best[:, 1:], v[:, e0 - 1:K1]], 1)
            alpha += post_t[t]
            alpha.clamp_(min=FLOOR)
            prev2[:, t + 1] = best_p
            tau += intron_q[t].index_select(1, phase)                          # (B, K, R)
            if len(pending) == m:
                entry = pending.pop(0)[:, :, None] + window_pi[t - m + 1].index_select(1, phase)
                enter = entry > tau
                entered[:, t + 1] = enter
                tau = torch.where(enter, entry, tau)
            tau.clamp_(min=FLOOR)
    z = alpha[:, u]
    score = torch.where(z > FLOOR / 2, z, torch.full_like(z, float("-inf")))
    # Expand the two argmax slots to the dense ``prev`` layout once, after
    # the loop: single-predecessor states point at their one predecessor.
    prev = torch.empty((B, L + 1, K), dtype=torch.int8, device=device)
    prev[:, :, 1:e0] = gather[:e0 - 1].to(torch.int8)
    prev[:, :, e0 + 1:] = gather[e0 - 1:K1].to(torch.int8)
    chosen = sp.cand.gather(1, prev2.long().view(-1, 2).T).T.view(B, L + 1, 2)
    prev[:, :, u] = chosen[..., 0]
    prev[:, :, e0] = chosen[..., 1]
    return score, prev, exit_r, entered


def viterbi_batch_reference(emissions: torch.Tensor, symbols: torch.Tensor,
                            lengths: torch.Tensor, grammar: Grammar, tables=None):
    """The dense-transition scan :func:`viterbi_batch` replaced: the same
    contract, with the per-step ``(B, K, K)`` block ``merged + prior_max[x[t]]``
    maximised over all ``K`` predecessors. Kept as the test oracle of the
    sparse scan and for profiling; not used by the pipeline."""
    B, C, L = emissions.shape
    if C != EMISSION_CHANNELS:
        raise ValueError(f"emissions must be (B, {EMISSION_CHANNELS}, L)")
    g = grammar
    K, R, m, u = g.K, g.R, g.dur.m, g.u_index
    dtype, device = emissions.dtype, emissions.device
    sym_t, base_t, ucol_t, post_t, donor_t, acceptor_t, intron_q, window_pi = _operands(
        emissions, symbols, lengths, g, tables)
    ucol_t = tuple(torch.zeros((B, K), dtype=dtype, device=device).index_copy_(1, _sparse(g).cand[0], c)
                   for c in ucol_t)
    prior = g.prior_max
    log_1mq = duration_by_state(g, tables)[2][None]
    phase = g.phase
    neg_one = torch.full((B, K), -1, dtype=torch.long, device=device)

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
            # ``copy_`` into the int8 stores converts; no separate ``to``.
            exit_r[:, t] = torch.where(from_tail, best_r, neg_one)
            merged = torch.where(from_tail, exits, alpha) + base_t[t]
            trans = merged[:, :, None] + prior[sym_t[t]]                      # (B, K, K)
            trans[:, :, u] += ucol_t[t]
            alpha, best_i = trans.max(dim=1)
            alpha += post_t[t]
            alpha.clamp_(min=FLOOR)
            prev[:, t + 1] = best_i
            tau += intron_q[t].index_select(1, phase)                          # (B, K, R)
            if len(pending) == m:
                entry = pending.pop(0)[:, :, None] + window_pi[t - m + 1].index_select(1, phase)
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
            duration: DurationMixture = DurationMixture(),
            tables=None) -> Tuple[float, List[Chain]]:
    """Best complete-path score and its chains for one window; the tensor twin
    of ``DelayedEntryDecoder(code, duration).viterbi(x, scores)`` without edge
    priors. Returns ``(-inf, [])`` when the grammar admits no path. ``tables``
    overrides ``duration``'s values (``duration`` then fixes ``m`` and ``R``)."""
    _check_input(x, emissions)
    g = Grammar.get(code, duration, dtype=emissions.dtype, device=emissions.device)
    sym = symbol_index_tensor(x)[None].to(emissions.device)
    lengths = torch.tensor([len(x)], device=emissions.device)
    score, prev, exit_r, entered = viterbi_batch(emissions[None], sym, lengths, g, tables)
    best = float(score[0])
    if best == float("-inf"):
        return best, []
    states = traceback(len(x), g, prev[0], exit_r[0], entered[0])
    return best, ReferenceDecoder._chains(states, [None] * len(x), x)


def viterbi_windows(windows: Sequence[str], emissions: Sequence[torch.Tensor], *,
                    codes: Optional[Sequence[GeneticCode]] = None,
                    duration: DurationMixture = DurationMixture(), tables=None,
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
                    sym[b, :lengths[b]] = symbol_index_tensor(windows[i]).to(device)
            score, prev, exit_r, entered = viterbi_batch(
                e, sym, torch.tensor(lengths, device=device), g, tables)
            for b, i in enumerate(ids):
                best = float(score[b])
                if best == float("-inf"):
                    out[i] = (best, [])
                    continue
                states = traceback(lengths[b], g, prev[b], exit_r[b], entered[b])
                out[i] = (best, ReferenceDecoder._chains(states, [None] * lengths[b], windows[i]))
    return out  # type: ignore[return-value]
