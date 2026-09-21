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
  ``T(i, r)`` whose exit at ``t`` beat it (``R + r`` when the exiting tail
  is the residual intron ``J(i, r)`` of the edge mode);
* ``entered[t + 1, c, r]`` (bool): tail ``T(c, r)`` at ``t + 1`` was entered
  by the donor parked at boundary ``t + 1 - m`` (its ``m - 1`` mandatory
  intronic positions are re-expanded to ``I(c, k)`` states on traceback,
  exactly as the Python decoder does), else it continued.

Storage is ``L * K * (2 + R)`` bytes per window in this dense layout (120 at
``K = 24``, ``R = 3``); the tiled scan over segments (:func:`scan_segments`)
holds them packed instead (:class:`PackedBackPointers`, 23 bytes per base)
and expands one row at a time for its traceback.

Complete-target only by default, like the training kernel: the path must
end in ``U`` at boundary ``n``. :func:`viterbi_batch_edges` adds the
sequence-edge partials of proposal 3.1 under an :class:`EdgePrior` (the four
partial-family scalars, :func:`model.a.pooled.edge_prior`): at boundary 0
every ``E(q)`` may be entered as ``E0(q)`` (coding entry plus the
normalized phase/prefix prior; its first base must be CDS, so no donor at
``t = 0``) and every tail ``T(E(q), r)`` as the residual intron ``J``
(intron entry plus the prefix prior and ``log pi``; it must consume one
intronic base before it may close, so no acceptor at ``t = 0``); at a
window's own end ``n`` any ``S``/``E`` state may exit with the coding exit,
any tail ``T`` that was really entered by a donor with the intron exit and
no ``(1 - q)`` factor, and any donor parked at ``s > n - m`` as the
censored ``I(c, n - s)`` with its intronic bases summed and the intron
exit. ``J`` is never terminal, and it lives in its own ``(B, K, R)`` layer:
it advances with the same intronic emissions as the donor-entered tails and
competes with them at every acceptor, but never receives donor entries and
is not offered at ``n``. Folding it into the tail layer (the first edge
implementation) lost the best donor-entered tail whenever ``J`` outscored
the entry on its row, and with it the window's only valid terminal path
(engels-0088, stalin-0090; ``tests/test_a_fast_viterbi.py::EdgeParity``
``test_terminal_tail_survives_residual_intron``). The extra layer costs one
add and one max per step in edge mode only; the per-window final choice is
made once, at the step that reaches its length. Ties are broken by the
first maximal index, which can
differ from the Python decoder's insertion order; the score is tie-free.
Back-pointers of states no finite path reaches, and of boundaries past a
window's end inside a padded batch, are unspecified (the dense reference
fills them with the identity; the sparse scan does not).
"""
from __future__ import annotations

from itertools import groupby
from typing import List, NamedTuple, Optional, Sequence, Tuple

import torch

from model.grammar import DurationMixture, EdgePrior, GeneticCode, TABLES
from model.grammar.reference import Chain, ReferenceDecoder, prefix_prior

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
    """Max-product scan over complete paths (``U`` to ``U``). Returns
    ``(score, prev, exit_r, entered)``: ``score``
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
    return _scan(emissions, symbols, lengths, grammar, tables, None)[:4]


# Final-state kinds of :func:`viterbi_batch_edges`: ``final[b] = (kind, j, x)``.
FINAL_COMPLETE, FINAL_CODING, FINAL_TAIL, FINAL_PENDING = 0, 1, 2, 3


def viterbi_batch_edges(emissions: torch.Tensor, symbols: torch.Tensor,
                        lengths: torch.Tensor, grammar: Grammar, edges: EdgePrior,
                        tables=None):
    """:func:`viterbi_batch` with the sequence-edge partials of proposal 3.1
    (module docstring) under ``edges``. Returns ``(score, prev, exit_r,
    entered, final)``: ``score`` is the best complete-or-partial path (an
    empty window scores 0 with no chain), and ``final`` is ``(B, 3)`` int64,
    the state the best path holds at the window's own boundary ``n``:
    ``(FINAL_COMPLETE, U, -1)``, ``(FINAL_CODING, j, -1)`` for an ``S``/``E``
    exit, ``(FINAL_TAIL, c, r)`` for a ``T(c, r)`` exit, or
    ``(FINAL_PENDING, c, k)`` for the censored ``I(c, k)`` of the donor parked
    at ``n - k``. :func:`traceback` takes it as ``final=``; the entry
    states at boundary 0 come out as ``("E0", q)`` and ``("J", c, r)``
    (``exit_r`` holds ``R + r`` where ``J`` closed)."""
    return _scan(emissions, symbols, lengths, grammar, tables, edges)


class Carry(NamedTuple):
    """Scan state at a tile seam, carried into the next tile of the same
    segments (proposal 3.3): the coding layer ``alpha`` ``(B, K)``, the
    donor-entered tails ``tau`` ``(B, K, R)``, the residual-intron layer
    ``tau_j`` (edge mode, else ``None``), the ``pending`` donors parked at the
    last ``len(pending) <= m - 1`` boundaries (each ``(B, K)``; their ``m``
    mandatory intronic bases straddle the seam, so the next tile re-reads the
    matching ``context`` emissions ``(B, C, len(pending))`` and starts its
    loop at ``t0 = len(pending)``), the running edge-mode ``score`` /
    ``final``, the tails a donor has entered so far (``ever`` ``(B, K, R)``;
    ``_finish`` offers only those) and the segment ``done`` mask (rows whose
    own end has been reached; their later tiles have length 0)."""
    alpha: torch.Tensor
    tau: torch.Tensor
    tau_j: Optional[torch.Tensor]
    pending: List[torch.Tensor]
    context: torch.Tensor
    score: Optional[torch.Tensor]
    final: Optional[torch.Tensor]
    ever: Optional[torch.Tensor]
    done: torch.Tensor


def _scan(emissions: torch.Tensor, symbols: torch.Tensor, lengths: torch.Tensor,
          grammar: Grammar, tables, edges: Optional[EdgePrior],
          carry: Optional[Carry] = None, keep: bool = False):
    """One scan over a batch of windows (``carry`` is ``None``) or over the
    next tile of a batch of segments whose state at the seam is ``carry``
    (:func:`scan_segments`). Returns ``(score, prev, exit_r, entered, final)``
    and, with ``keep``, the :class:`Carry` at the tile's end as a sixth
    element. Under a ``carry`` the back-pointer tensors cover the tile's own
    boundaries ``0..L`` (boundary 0 is the seam) and ``exit_r[:, L]`` is
    unwritten (the seam exit is written by the next tile at its boundary 0);
    the boundary-0 entries and the ``t = 0`` floors of the edge mode apply to
    the first tile only, and ``_finish`` runs for a row at the tile in which
    its segment ends."""
    B, C, L = emissions.shape
    if C != EMISSION_CHANNELS:
        raise ValueError(f"emissions must be (B, {EMISSION_CHANNELS}, L)")
    g = grammar
    K, R, m, u = g.K, g.R, g.dur.m, g.u_index
    dtype, device = emissions.dtype, emissions.device
    t0 = 0
    if carry is not None:
        # Prepend the seam context so the pending donors' mandatory windows
        # (``window_pi[t - m + 1]`` spans ``[t - m + 1, t]``) and the censored
        # suffix of ``_finish`` index the same operands as an unbroken scan.
        t0 = carry.context.shape[-1]
        if len(carry.pending) != t0:
            raise ValueError("carry context must match its pending donors")
        emissions = torch.cat([carry.context, emissions], 2)
        symbols = torch.cat([torch.zeros((B, t0), dtype=symbols.dtype, device=device), symbols], 1)
        lengths = torch.where(lengths > 0, lengths + t0, lengths)
        L += t0
    sym_t, base_t, ucol_t, post_t, donor_t, acceptor_t, intron_q, window_pi = _operands(
        emissions, symbols, lengths, g, tables)
    sp = _sparse(g)
    gather, prior, e0 = sp.gather, sp.prior, sp.e0
    K1, P = K - 2, sp.cand.shape[1]
    log_pi, _, log_1mq = duration_by_state(g, tables)
    log_1mq = log_1mq[None]
    phase = g.phase
    neg_one = torch.full((B, K), -1, dtype=torch.long, device=device)

    prev2 = torch.zeros((B, L + 1, 2), dtype=torch.int8, device=device)
    exit_r = torch.full((B, L + 1, K), -1, dtype=torch.int8, device=device)
    entered = torch.zeros((B, L + 1, K, R), dtype=torch.bool, device=device)

    if carry is None:
        alpha = torch.full((B, K), FLOOR, dtype=dtype, device=device)
        alpha[:, u] = 0.0
        tau = torch.full((B, K, R), FLOOR, dtype=dtype, device=device)
        tau_j = None
        pending: List[torch.Tensor] = []
        done = torch.zeros((B,), dtype=torch.bool, device=device)
        ever = None
    else:
        alpha, tau, tau_j = carry.alpha.clone(), carry.tau.clone(), carry.tau_j
        tau_j = None if tau_j is None else tau_j.clone()
        pending = list(carry.pending)
        done, ever = carry.done, carry.ever
    ends = None
    final = None
    score = None
    if edges is not None:
        if carry is None:
            # Entry at boundary 0: E0(q) in the coding layer, J(E(q), r) in its
            # own tail layer, kept apart from the donor-entered ``tau`` because
            # the two differ at the window's end (``T`` may exit there, ``J`` may
            # not) and a ``J`` that outscored a donor entry would otherwise erase
            # the only terminal candidate on that row (engels-0088, stalin-0090).
            is_e = torch.tensor([st[0] == "E" for st in g.states], device=device)
            pp = torch.tensor([prefix_prior(st[1]) if st[0] == "E" else 0.0 for st in g.states],
                              dtype=dtype, device=device)
            alpha = torch.where(is_e, pp + edges.entry, alpha)                                # (B, K)
            tau_j = torch.where(is_e[:, None], pp[:, None] + log_pi.to(dtype) + edges.intron_entry,
                                tau)                                                          # (B, K, R)
            if L:
                donor_t[0].fill_(FLOOR)            # E0's first base is CDS: no donor at t = 0
                acceptor_t[0].fill_(FLOOR)         # J consumes one intronic base before closing
            score = torch.full((B,), FLOOR, dtype=dtype, device=device)
            final = torch.zeros((B, 3), dtype=torch.long, device=device)
            final[:, 1], final[:, 2] = u, -1
        else:
            score, final = carry.score.clone(), carry.final.clone()
        valid = torch.arange(L, device=device)[None, :] < lengths[:, None]
        neg = torch.full((B, L), FLOOR, dtype=dtype, device=device)
        intron = torch.where(valid[:, None, :], emissions[:, 4:7, :], neg[:, None, :])  # (B, 3, L)
        exit_c = torch.full((K,), edges.exit, dtype=dtype, device=device)
        exit_c[u] = FLOOR
        exit_i = torch.tensor(edges.intron_exit, dtype=dtype, device=device)
        ends = {}
        for b, n in enumerate(lengths.tolist()):
            if n > 0 or carry is None:
                ends.setdefault(n, []).append(b)
        ends = {n: torch.tensor(rows, device=device) for n, rows in ends.items()}
        if 0 in ends:
            score[ends.pop(0)] = 0.0           # only U -> U: partition 1, no chain
    with torch.no_grad():
        for t in range(t0, L):
            pending.append(alpha + donor_t[t])
            exits, best_r = (tau + log_1mq).max(dim=-1)                       # (B, K)
            if tau_j is not None:
                # The residual intron closes like a tail; ``exit_r`` records
                # it as ``R + r`` so the traceback knows it reaches boundary 0.
                exits_j, best_rj = (tau_j + log_1mq).max(dim=-1)
                from_j = exits_j > exits
                exits = torch.where(from_j, exits_j, exits)
                best_r = torch.where(from_j, best_rj + R, best_r)
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
            step = intron_q[t].index_select(1, phase)                          # (B, K, R)
            tau += step
            if tau_j is not None:
                tau_j += step
                tau_j.clamp_(min=FLOOR)
            if len(pending) == m:
                entry = pending.pop(0)[:, :, None] + window_pi[t - m + 1].index_select(1, phase)
                enter = entry > tau
                entered[:, t + 1] = enter
                tau = torch.where(enter, entry, tau)
            tau.clamp_(min=FLOOR)
            if ends is not None and t + 1 in ends:
                rows = ends[t + 1]
                _finish(rows, t + 1, alpha, tau, pending, entered, intron, phase, m, u,
                        exit_c, exit_i, score, final, ever)
    if keep:
        ctx = len(pending)
        if edges is not None:
            seen = entered[:, t0 + 1:].any(dim=1)
            ever = seen if ever is None else ever | seen
        carry_out = Carry(alpha=alpha, tau=tau, tau_j=tau_j, pending=pending,
                          context=emissions[:, :, L - ctx:L].clone() if ctx else emissions[:, :, :0].clone(),
                          score=score, final=final, ever=ever, done=done | (lengths < L))
    if ends is None:
        z = alpha[:, u]
        score = torch.where(z > FLOOR / 2, z, torch.full_like(z, float("-inf")))
    else:
        score = torch.where(score > FLOOR / 2, score, torch.full_like(score, float("-inf")))
    prev2, exit_r, entered = prev2[:, t0:], exit_r[:, t0:], entered[:, t0:]
    if keep:
        # The tiled scan packs the raw two-slot ``prev2`` and expands one row
        # at a time on traceback (:class:`PackedBackPointers`).
        return score, prev2, exit_r, entered, final, carry_out
    return score, _expand_prev(sp, prev2), exit_r, entered, final


def _expand_prev(sp: "_Sparse", prev2: torch.Tensor) -> torch.Tensor:
    """The dense ``prev`` ``(..., K)`` from the two argmax slots ``prev2``
    ``(..., 2)`` of ``U`` and ``E("")``: single-predecessor states point at
    their one predecessor."""
    K, e0 = sp.gather.shape[0] - 2 * sp.cand.shape[1] + 2, sp.e0             # gather is (K - 2 + 2 P,)
    K1 = K - 2
    lead = prev2.shape[:-1]
    prev = torch.empty(lead + (K,), dtype=torch.int8, device=prev2.device)
    prev[..., 1:e0] = sp.gather[:e0 - 1].to(torch.int8)
    prev[..., e0 + 1:] = sp.gather[e0 - 1:K1].to(torch.int8)
    chosen = sp.cand.gather(1, prev2.reshape(-1, 2).long().T).T.view(lead + (2,))
    prev[..., 0] = chosen[..., 0]
    prev[..., e0] = chosen[..., 1]
    return prev


def _finish(rows, n, alpha, tau, pending, entered, intron, phase, m, u, exit_c, exit_i,
            score, final, ever=None):
    """Best terminal choice at boundary ``n`` for the windows ``rows`` of
    that length (edge mode): complete, coding exit, tail exit (``tau`` holds
    donor-entered tails only; rows no donor ever entered are floored, and the
    residual intron ``J`` lives in its own layer and is never terminal), or a
    censored pending donor. ``ever`` is the entered mask carried from earlier
    tiles of the same scan (:func:`scan_segments`)."""
    a = alpha[rows]                                                       # (b, K)
    cands = [a[:, u], (a + exit_c).max(dim=1)]
    from_tail = ~entered[rows, 1:n + 1].any(dim=1)                        # (b, K, R): never entered
    if ever is not None:
        from_tail &= ~ever[rows]
    tail = torch.where(from_tail, torch.full_like(tau[rows], FLOOR), tau[rows] + exit_i)
    cands.append(tail.flatten(1).max(dim=1))
    # Donors parked at s in (n - m, n): I(c, n - s) exits with the intronic
    # bases s .. n - 1 summed (a -inf inside stays -inf under the reversed cumsum).
    first = n - len(pending)
    suffix = intron[rows, :, first:n].flip(-1).cumsum(-1).flip(-1)         # (b, 3, len(pending))
    pend = [(p[rows] + suffix[:, phase, i] + exit_i).max(dim=1) for i, p in enumerate(pending)]
    if pend:
        pv = torch.stack([v for v, _ in pend], 1)                         # (b, len(pending))
        pj = torch.stack([j for _, j in pend], 1)
        best_pv, best_i = pv.max(dim=1)
        cands.append((best_pv, best_i))
    values = torch.stack([cands[0], cands[1][0], cands[2][0]] + ([cands[3][0]] if pend else []), 1)
    best, kind = values.max(dim=1)
    score[rows] = best
    R = tau.shape[-1]
    j = torch.where(kind == FINAL_CODING, cands[1][1],
                    torch.where(kind == FINAL_TAIL, cands[2][1] // R, torch.full_like(kind, u)))
    x = torch.where(kind == FINAL_TAIL, cands[2][1] % R, torch.full_like(kind, -1))
    if pend:
        pj_best = pj.gather(1, best_i[:, None])[:, 0]
        j = torch.where(kind == FINAL_PENDING, pj_best, j)
        x = torch.where(kind == FINAL_PENDING, n - (first + best_i), x)
    final[rows, 0], final[rows, 1], final[rows, 2] = kind, j, x


class PackedBackPointers:
    """The back-pointers of a batch of segments, packed per tile as the scan
    produces them and expanded one row at a time for :func:`traceback`:
    ``prev`` as the two argmax slots of ``U`` and ``E("")`` (2 bytes per
    boundary; the other states have one predecessor), ``exit_r`` as nibbles
    (``K / 2`` bytes when ``2 R < 16``, else the int8 column), ``entered`` as
    ``K R`` bits (``ceil(K R / 8)`` bytes). At ``K = 24``, ``R = 3`` this is
    23 bytes per base instead of the 120 of the dense tensors. ``dense(b)``
    returns ``(prev, exit_r, entered)`` of row ``b`` in the dense layout the
    traceback reads (``(L + 1, K)`` int8, ``(L + 1, K)`` int8, ``(L + 1, K,
    R)`` bool), so one row of dense pointers is alive during its traceback,
    never the batch."""

    def __init__(self, sp: "_Sparse", K: int, R: int, device):
        self.sp, self.K, self.R = sp, K, R
        self.nibble = 2 * R < 16
        self.prev2: List[torch.Tensor] = []
        self.exit_r: List[torch.Tensor] = []
        self.entered: List[torch.Tensor] = []
        self.bits = torch.tensor([1, 2, 4, 8, 16, 32, 64, 128], dtype=torch.uint8, device=device)

    def append(self, prev2: torch.Tensor, exit_r: torch.Tensor, entered: torch.Tensor,
               first: bool) -> None:
        """Pack one tile's ``(B, T + 1, ...)`` pointers. A tile after the
        first shares its boundary 0 with the previous tile's last boundary:
        ``prev`` / ``entered`` at the seam are unread (the traceback reads
        ``prev[t]`` for ``t >= 1``) and kept once; ``exit_r`` at the seam is
        written by the later tile, so the earlier tile's last column is
        dropped when the next arrives."""
        if not first:
            prev2, entered = prev2[:, 1:], entered[:, 1:]
            self.exit_r[-1] = self.exit_r[-1][:, :-1]
        self.prev2.append(prev2.contiguous())
        self.exit_r.append(self._pack_exit(exit_r))
        self.entered.append(self._pack_bits(entered))

    def _pack_exit(self, exit_r: torch.Tensor) -> torch.Tensor:
        if not self.nibble:
            return exit_r.contiguous()
        v = (exit_r + 1).to(torch.uint8)                                        # -1 .. 2R-1 -> 0 .. 2R
        if v.shape[-1] % 2:
            v = torch.cat([v, torch.zeros_like(v[..., :1])], -1)
        return v[..., 0::2] | (v[..., 1::2] << 4)

    def _unpack_exit(self, packed: torch.Tensor) -> torch.Tensor:
        if not self.nibble:
            return packed
        v = torch.stack([packed & 15, packed >> 4], -1).flatten(-2)[..., :self.K]
        return v.to(torch.int8) - 1

    def _pack_bits(self, entered: torch.Tensor) -> torch.Tensor:
        flat = entered.flatten(-2).to(torch.uint8)                             # (B, T, K R)
        nbits = flat.shape[-1]
        pad = (-nbits) % 8
        if pad:
            flat = torch.cat([flat, torch.zeros(flat.shape[:-1] + (pad,), dtype=torch.uint8,
                                                device=flat.device)], -1)
        return (flat.view(flat.shape[:-1] + (-1, 8)) * self.bits).sum(-1, dtype=torch.uint8)

    def _unpack_bits(self, packed: torch.Tensor) -> torch.Tensor:
        flat = (packed[..., None] & self.bits) != 0
        return flat.flatten(-2)[..., :self.K * self.R].view(packed.shape[:-1] + (self.K, self.R))

    def dense(self, b: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        prev2 = torch.cat([p[b] for p in self.prev2], 0) if len(self.prev2) > 1 else self.prev2[0][b]
        exit_r = torch.cat([self._unpack_exit(p[b]) for p in self.exit_r], 0)
        entered = torch.cat([self._unpack_bits(p[b]) for p in self.entered], 0)
        return _expand_prev(self.sp, prev2), exit_r, entered

    def nbytes(self) -> int:
        """Bytes held by the packed store."""
        return sum(t.numel() * t.element_size()
                   for ts in (self.prev2, self.exit_r, self.entered) for t in ts)


def scan_segments(emissions: torch.Tensor, symbols: torch.Tensor, lengths: torch.Tensor,
                  grammar: Grammar, tile: int, edges: Optional[EdgePrior] = None, tables=None):
    """:func:`viterbi_batch` / :func:`viterbi_batch_edges` over a batch of
    *segments* scanned ``tile`` bases at a time with the state carried across
    the seams (proposal 3.3): the scores and back-pointers are those of one
    unbroken scan over each segment, but at most one tile of operands is
    alive at a time, so host memory for the operands is bounded by ``B *
    tile`` rather than the segment length, and the back-pointers of the
    whole segment are held packed (:class:`PackedBackPointers`, 23 bytes
    per base at ``K = 24``, ``R = 3``) rather than dense (120). ``emissions``
    is ``(B, C, L)`` with each row's segment in ``[0, lengths[b])``; interior
    seams are not sequence edges, so a gene crossing one is decoded whole,
    and ``edges`` applies to the two real ends of each segment only. Returns
    ``(score, packed, final)``; ``packed.dense(b)`` gives row ``b``'s
    ``(prev, exit_r, entered)`` for :func:`traceback`."""
    if tile < 1:
        raise ValueError("tile must be positive")
    B, C, L = emissions.shape
    device = emissions.device
    carry = None
    packed = PackedBackPointers(_sparse(grammar), grammar.K, grammar.R, device)
    score = final = None
    for start in range(0, max(L, 1), tile):
        end = min(start + tile, L)
        chunk = (lengths - start).clamp(min=0, max=end - start)
        out = _scan(emissions[:, :, start:end], symbols[:, start:end], chunk, grammar, tables,
                    edges, carry, keep=True)
        score, prev2, exit_r, entered, final, carry = out
        packed.append(prev2, exit_r, entered, first=start == 0)
    # A row that ended before the last tile kept ``alpha[U]`` through the pad
    # steps of its later tiles (``U -> U`` at 0), so the free-grammar score
    # read at the end is its own, as in a padded batch.
    return score, packed, final


def viterbi_segments(segments: Sequence[str], emissions: Sequence[torch.Tensor], *,
                     tile: int, code: GeneticCode = TABLES[1],
                     duration: DurationMixture = DurationMixture(), tables=None,
                     edges: Optional[EdgePrior] = None) -> List[Tuple[float, List[Chain]]]:
    """:func:`viterbi` over a batch of segments of one genetic code, scanned
    ``tile`` bases at a time with the state carried across seams
    (:func:`scan_segments`); ``emissions[i]`` is ``(11, len(segments[i]))``.
    Equal to :func:`viterbi` on each whole segment."""
    if len(segments) != len(emissions):
        raise ValueError("segments and emissions must agree in length")
    if not segments:
        return []
    dtype, device = emissions[0].dtype, emissions[0].device
    g = Grammar.get(code, duration, dtype=dtype, device=device)
    lengths = [len(x) for x in segments]
    L = max(lengths)
    e = torch.full((len(segments), emissions[0].shape[0], L), float("-inf"), dtype=dtype, device=device)
    sym = torch.zeros((len(segments), L), dtype=torch.long, device=device)
    for b, (x, em) in enumerate(zip(segments, emissions)):
        _check_input(x, em)
        e[b, :, :lengths[b]] = em
        if lengths[b]:
            sym[b, :lengths[b]] = symbol_index_tensor(x).to(device)
    score, packed, final = scan_segments(
        e, sym, torch.tensor(lengths, device=device), g, tile, edges, tables)
    out = []
    for b, x in enumerate(segments):
        best = float(score[b])
        if best == float("-inf"):
            out.append((best, []))
            continue
        states = traceback(lengths[b], g, *packed.dense(b),
                           None if final is None else final[b])
        out.append((best, ReferenceDecoder._chains(states, [None] * lengths[b], x)))
    return out


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
              entered: torch.Tensor, final=None) -> List[tuple]:
    """State tuples at boundaries ``0..n`` of one window from its back-pointers
    (each ``(L + 1, ...)``), in the reference decoder's vocabulary:
    ``("U",)``, ``("S", q)``, ``("E", q)``, ``("I", c, k)``, ``("T", c, r)``,
    and with ``final`` (one row of :func:`viterbi_batch_edges`) also the
    edge states ``("E0", q)`` and ``("J", c, r)`` at the start."""
    g = grammar
    m = g.dur.m
    # numpy scalar reads are ~10x cheaper than tensor indexing or ``tolist``.
    prev, exit_r, entered = (a.cpu().numpy() for a in (prev, exit_r, entered))
    states: List[Optional[tuple]] = [None] * (n + 1)
    t, j, r = n, g.u_index, -1                        # r >= 0: in tail T(j, r)
    if final is not None:
        kind, fj, fx = (int(v) for v in final)
        if kind == FINAL_CODING:
            j = fj
        elif kind == FINAL_TAIL:
            j, r = fj, fx
        elif kind == FINAL_PENDING:                   # censored I(c, 1..k) from the donor at n - k
            for k in range(1, fx + 1):
                states[n - fx + k] = ("I", g.states[fj], k)
            t, j = n - fx, fj
    R = entered.shape[-1]
    while t > 0:
        if r < 0:
            states[t] = g.states[j]
            i = int(prev[t, j])
            t -= 1
            j, r = i, int(exit_r[t, i])
        elif r >= R:                                  # the residual intron J: boundary 0 to here
            for k in range(0, t + 1):
                states[k] = ("J", g.states[j], r - R)
            return states
        else:
            states[t] = ("T", g.states[j], r)
            if entered[t, j, r]:
                s = t - m
                for k in range(s + 1, t):
                    states[k] = ("I", g.states[j], k - s)
                t, r = s, -1
            else:
                t -= 1
    if r >= 0:
        raise ValueError("traceback reached boundary 0 inside a tail no donor entered")
    states[0] = g.states[j] if j == g.u_index else ("E0", g.states[j][1])
    return states


def viterbi(x: str, emissions: torch.Tensor, *, code: GeneticCode = TABLES[1],
            duration: DurationMixture = DurationMixture(),
            tables=None, edges: Optional[EdgePrior] = None) -> Tuple[float, List[Chain]]:
    """Best path score and its chains for one window; the tensor twin of
    ``DelayedEntryDecoder(code, duration, edges).viterbi(x, scores)``:
    complete paths only without ``edges``, sequence-edge partials with them.
    Returns ``(-inf, [])`` when the grammar admits no path. ``tables``
    overrides ``duration``'s values (``duration`` then fixes ``m`` and ``R``)."""
    _check_input(x, emissions)
    g = Grammar.get(code, duration, dtype=emissions.dtype, device=emissions.device)
    sym = symbol_index_tensor(x)[None].to(emissions.device)
    lengths = torch.tensor([len(x)], device=emissions.device)
    score, prev, exit_r, entered, final = _scan(emissions[None], sym, lengths, g, tables, edges)
    best = float(score[0])
    if best == float("-inf"):
        return best, []
    states = traceback(len(x), g, prev[0], exit_r[0], entered[0],
                       None if final is None else final[0])
    return best, ReferenceDecoder._chains(states, [None] * len(x), x)


def viterbi_windows(windows: Sequence[str], emissions: Sequence[torch.Tensor], *,
                    codes: Optional[Sequence[GeneticCode]] = None,
                    duration: DurationMixture = DurationMixture(), tables=None,
                    batch_size: int = 16,
                    edges: Optional[EdgePrior] = None) -> List[Tuple[float, List[Chain]]]:
    """:func:`viterbi` over many windows, batched to amortize the per-step
    launch cost of the scan (on one CPU core the scan over a ``(B, K)`` layer
    costs about the same as over ``(1, K)``, so decoding ``B`` windows
    together is nearly ``B`` times cheaper per base). Windows are grouped by
    genetic code and sorted by length so padding stays small; each batch is
    padded to its longest window with ``-inf`` emissions past the end, which
    :func:`model.a.fast_loss.prepare_batch` turns into identity transitions.
    ``emissions[i]`` is ``(11, len(windows[i]))``. Results are returned in the
    input order. ``edges`` enables the sequence-edge partials for every
    window (each window is then a whole sequence with two real edges)."""
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
            score, prev, exit_r, entered, final = _scan(
                e, sym, torch.tensor(lengths, device=device), g, tables, edges)
            for b, i in enumerate(ids):
                best = float(score[b])
                if best == float("-inf"):
                    out[i] = (best, [])
                    continue
                states = traceback(lengths[b], g, prev[b], exit_r[b], entered[b],
                                   None if final is None else final[b])
                out[i] = (best, ReferenceDecoder._chains(states, [None] * lengths[b], windows[i]))
    return out  # type: ignore[return-value]
