"""Vectorized delayed-entry forward for candidate A's chain loss (proposal 3.2,
3.3; the training kernel the smoke profile of a-pilot.md section 3.1 showed is
mandatory).

The reference torch forward in :mod:`model.a.torch_loss` walks the expanded
grammar (every mandatory intronic position is a state, every boundary a Python
dict) and costs 4-7 s/kb and ~0.8 MB/base on a training window. This module is
the same partition function, ``log Z``, computed with the *delayed-entry*
recurrence of :class:`model.grammar.DelayedEntryDecoder` on dense tensors:

* the coding layer ``alpha[t]`` is one vector over the ``K`` coding states
  (``U``, the initiator prefixes ``S(q)``, the 21 ordinary prefixes ``E(q)``),
* the intron tails ``tau[t]`` are one ``(K, R)`` matrix (tail of coding state
  ``c`` in duration component ``r``; the ``U`` row is never live),
* a donor parked at boundary ``s`` is the vector ``alpha[s] + donor[s]``; it
  enters every tail exactly ``m`` boundaries later with the rolling
  ``m``-window sum of the phase's intron emissions and ``log pi[p, r]``
  (the ``m-1`` mandatory positions are never visited).

The per-base coding transition is a ``(K, K)`` log-matrix ``L[t]`` built once
per window from a symbol-indexed prior table (uniform over the permitted bases
of an IUPAC code, ``-log k``, summed over bases that reach the same successor)
plus a fixed ``(K, K, 11)`` channel-coefficient tensor contracted with the
emissions -- so the whole ``(n, K, K)`` transition stack is one ``einsum``.
The scan itself is a Python loop of ``n`` steps over ``(B, K)`` and
``(B, K, R)`` tensors and runs batched over windows: padded positions get the
identity transition and floored intron/donor/acceptor emissions, so a shorter
window's ``alpha[U]`` is frozen at its own end.

Numerics: ``-inf`` never enters the scan. Emissions are clamped at
:data:`FLOOR`, so every forbidden path has weight exactly ``exp(FLOOR - z) = 0``
against any finite competitor and autograd assigns it zero posterior, while a
window whose support admits *no* legal path yields ``log Z <= FLOOR / 2`` and is
reported as ``-inf`` (``torch.logsumexp`` over an all ``-inf`` layer would
backpropagate NaN). Parity with :func:`model.a.torch_loss.chain_nll` and its
gradient is pinned by ``tests/test_a_fast_loss.py`` on the section-3.4
fixtures, random emissions, ambiguous bases, ``m = 1`` and multi-component
mixtures, and batched-vs-single execution.

Complete-target only, like the reference torch loss: edge partials (``E0``,
``J``, censored exits) are the section-3.6 boundary-support increment and are
not offered here.
"""
from __future__ import annotations

from math import log
from typing import Dict, List, Optional, Sequence, Tuple

import torch

from model.grammar import DurationMixture, GeneticCode, TABLES
from model.grammar.codes import IUPAC, permitted_bases

from .torch_loss import CHANNEL_ORDER, EMISSION_CHANNELS, _check_input, support_mask

Range = Tuple[int, int]

# Finite stand-in for -inf inside the scan (see the module docstring).
FLOOR = -1e30

# Distinct permitted-base sets an observed symbol can denote; index 0 is the
# fully ambiguous "ACGT" so unknown symbols map there (``permitted_bases``).
SYMBOL_SETS: List[str] = ["ACGT"] + sorted({v for v in IUPAC.values()} - {"ACGT"})
_SYMBOL_INDEX: Dict[str, int] = {s: i for i, s in enumerate(SYMBOL_SETS)}

_CH = {name: i for i, name in enumerate(CHANNEL_ORDER)}


def symbol_indices(x: str) -> List[int]:
    """Map an oriented sequence to indices into :data:`SYMBOL_SETS`."""
    return [_SYMBOL_INDEX[permitted_bases(ch)] for ch in x]


def _build_symbol_table():
    """256-entry byte -> :data:`SYMBOL_SETS` index lookup (unknown bytes map to
    the fully ambiguous set, as :func:`permitted_bases` does)."""
    table = torch.zeros((256,), dtype=torch.long)
    for byte in range(256):
        table[byte] = _SYMBOL_INDEX[permitted_bases(chr(byte))]
    return table


_SYMBOL_TABLE = _build_symbol_table()


def symbol_index_tensor(x: str) -> torch.Tensor:
    """:func:`symbol_indices` as an int64 tensor, through a byte lookup
    (equal to the list reference; non-ASCII input raises)."""
    try:
        raw = x.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError("sequence must be ASCII") from exc
    if not raw:
        return torch.zeros((0,), dtype=torch.long)
    return _SYMBOL_TABLE[torch.frombuffer(bytearray(raw), dtype=torch.uint8).long()]


class Grammar:
    """The dense tables of one (genetic code, duration mixture) pair.

    ``states`` lists the coding states in index order: ``U`` first, then the
    initiator prefixes, then the 21 ordinary prefixes. ``prior[sym, i, j]`` is
    the symbol-conditional log transition prior, ``coef[i, j, c]`` the 0/1
    channel coefficients, so ``L[t] = prior[x[t]] + coef . e[:, t]``. The
    intron phase of state ``i`` is ``phase[i]``; ``coding[i]`` marks the
    ``S``/``E`` states that may open an intron.
    """

    _cache: Dict[Tuple[int, DurationMixture, torch.dtype, torch.device], "Grammar"] = {}

    def __init__(self, code: GeneticCode, duration: DurationMixture, *,
                 dtype=torch.float64, device=None):
        self.code, self.dur = code, duration
        prefixes = sorted(code.initiator_prefixes(), key=lambda q: (len(q), q))
        states: List[tuple] = [("U",)] + [("S", q) for q in prefixes]
        states.append(("E", ""))
        for a in "ACGT":
            states.append(("E", a))
        for a in "ACGT":
            for b in "ACGT":
                states.append(("E", a + b))
        self.states = states
        self.index = {s: i for i, s in enumerate(states)}
        K = len(states)
        self.K = K
        self.u_index = 0
        self.phase = torch.tensor([len(s[1]) if s[0] != "U" else 0 for s in states],
                                  device=device)
        self.coding = torch.tensor([s[0] in ("S", "E") for s in states], device=device)

        # Symbol-conditional prior (counts of bases reaching each successor).
        counts = torch.zeros((len(SYMBOL_SETS), K, K), dtype=dtype, device=device)
        coef = torch.zeros((K, K, EMISSION_CHANNELS), dtype=dtype, device=device)
        prefix_set = code.initiator_prefixes()
        for si, bases in enumerate(SYMBOL_SETS):
            for b in bases:
                for i, s in enumerate(states):
                    for j in self._successors(s, b, prefix_set):
                        counts[si, i, j] += 1.0
        prior = torch.where(counts > 0, torch.log(counts.clamp(min=1.0)),
                            torch.full_like(counts, FLOOR))
        for si, bases in enumerate(SYMBOL_SETS):
            prior[si] = torch.where(counts[si] > 0, prior[si] - log(len(bases)), prior[si])
        self.prior = prior
        # Max-product prior (:mod:`model.a.fast_viterbi`): an ambiguous symbol
        # whose ``k`` permitted bases reach the same successor contributes
        # ``-log k`` once (the reference decoder's per-base uniform prior), not
        # ``log(count) - log k``; ``U -> U`` does not branch over bases and
        # carries no prior. Row ``len(SYMBOL_SETS)`` is the identity, used
        # past a window's end.
        prior_max = torch.full((len(SYMBOL_SETS) + 1, K, K), FLOOR, dtype=dtype, device=device)
        for si, bases in enumerate(SYMBOL_SETS):
            prior_max[si] = torch.where(counts[si] > 0, torch.full_like(counts[si], -log(len(bases))),
                                        prior_max[si])
        prior_max[:, 0, 0] = 0.0
        prior_max[len(SYMBOL_SETS)].fill_diagonal_(0.0)
        self.prior_max = prior_max
        # Channel coefficients depend on (i, j) only.
        for i, s in enumerate(states):
            for j, s2 in enumerate(states):
                if s[0] == "U":
                    if s2[0] == "U":
                        coef[i, j, _CH["u"]] = 1.0
                    elif s2[0] == "S":
                        coef[i, j, _CH["start"]] = 1.0
                        coef[i, j, _CH["cds[0]"]] = 1.0
                else:
                    p = len(s[1])
                    coef[i, j, _CH[f"cds[{p}]"]] = 1.0
                    if s2[0] == "U":
                        coef[i, j, _CH["stop"]] = 1.0
        self.coef = coef

        R = duration.R
        self.R = R
        log_pi = torch.tensor([[duration.log_pi(p, r) for r in range(R)] for p in range(3)],
                              dtype=dtype, device=device)
        log_q = torch.tensor([[duration.log_q(p, r) for r in range(R)] for p in range(3)],
                             dtype=dtype, device=device)
        log_1mq = torch.tensor([[duration.log_1mq(p, r) for r in range(R)] for p in range(3)],
                               dtype=dtype, device=device)
        # The (3, R) per-phase tables and their per-state (K, R) copies,
        # indexed by the state's phase.
        self.phase_log_pi, self.phase_log_q, self.phase_log_1mq = log_pi, log_q, log_1mq
        self.log_pi = log_pi[self.phase]
        self.log_q = log_q[self.phase]
        self.log_1mq = log_1mq[self.phase]
        self.identity = torch.full((K, K), FLOOR, dtype=dtype, device=device)
        self.identity.fill_diagonal_(0.0)

    def _successors(self, s, b, prefix_set):
        """Coding successors of state ``s`` on concrete base ``b``
        (:meth:`ReferenceDecoder.transitions`, coding moves only)."""
        kind = s[0]
        if kind == "U":
            yield self.index[("U",)]
            if b in prefix_set:
                yield self.index[("S", b)]
            return
        q = s[1]
        codon = q + b
        if kind == "S":
            if len(codon) < 3:
                if codon in prefix_set:
                    yield self.index[("S", codon)]
            elif self.code.is_initiator(codon):
                yield self.index[("E", "")]
        else:
            if len(codon) < 3:
                yield self.index[("E", codon)]
            elif self.code.is_stop(codon):
                yield self.index[("U",)]
            else:
                yield self.index[("E", "")]

    @classmethod
    def get(cls, code: GeneticCode, duration: DurationMixture, *, dtype, device) -> "Grammar":
        # Keyed by the frozen ``GeneticCode`` *value*, not its table number:
        # ``TABLES[1]`` and ``TABLES[1].with_alternative_initiators()`` share a
        # table but differ in initiators, and must not share a grammar
        # (engels-0080, stalin-0081).
        key = (code, duration, dtype, torch.device(device) if device is not None else None)
        g = cls._cache.get(key)
        if g is None:
            g = cls(code, duration, dtype=dtype, device=device)
            cls._cache[key] = g
        return g


def _floor(x: torch.Tensor) -> torch.Tensor:
    return x.clamp(min=FLOOR)


def duration_by_state(grammar: Grammar, tables=None):
    """Per-state ``(K, R)`` ``log pi``, ``log q``, ``log (1 - q)``: the
    grammar's fixed :class:`DurationMixture` values, or the learned ``(3, R)``
    :class:`model.a.pooled.DurationTables` indexed by each state's phase (the
    gradient path to the pooled duration scalars)."""
    g = grammar
    if tables is None:
        return g.log_pi, g.log_q, g.log_1mq
    if tuple(tables.log_pi.shape) != (3, g.R):
        raise ValueError(f"duration tables must be (3, {g.R}), got {tuple(tables.log_pi.shape)}")
    return tuple(t.to(dtype=g.log_pi.dtype, device=g.log_pi.device)[g.phase] for t in tables)


def duration_by_phase(grammar: Grammar, tables=None):
    """The ``(3, R)`` per-phase ``log pi``, ``log q``, ``log (1 - q)`` behind
    :func:`duration_by_state` (the grammar's fixed law or the learned
    ``tables``), for kernels that expand to states per step."""
    g = grammar
    if tables is None:
        return g.phase_log_pi, g.phase_log_q, g.phase_log_1mq
    if tuple(tables.log_pi.shape) != (3, g.R):
        raise ValueError(f"duration tables must be (3, {g.R}), got {tuple(tables.log_pi.shape)}")
    return tuple(t.to(dtype=g.log_pi.dtype, device=g.log_pi.device) for t in tables)


def log_partition_batch(emissions: torch.Tensor, symbols: torch.Tensor,
                        lengths: torch.Tensor, grammar: Grammar, tables=None) -> torch.Tensor:
    """``log Z`` of the free grammar for a batch of windows.

    ``emissions`` is ``(B, 11, L)`` (finite or ``-inf``; ``-inf`` is floored),
    ``symbols`` ``(B, L)`` long indices into :data:`SYMBOL_SETS`, ``lengths``
    ``(B,)`` the real window lengths. Returns ``(B,)`` with ``-inf`` where the
    grammar admits no path. ``tables`` optionally replaces the grammar's fixed
    duration law with learned :class:`model.a.pooled.DurationTables`.
    """
    B, C, L = emissions.shape
    if C != EMISSION_CHANNELS:
        raise ValueError(f"emissions must be (B, {EMISSION_CHANNELS}, L)")
    g = grammar
    K, R, m = g.K, g.R, g.dur.m
    log_pi, log_q, log_1mq = duration_by_state(g, tables)
    dtype, device = emissions.dtype, emissions.device
    e = _floor(emissions)
    valid = (torch.arange(L, device=device)[None, :] < lengths[:, None])      # (B, L)

    # (B, L, K, K) coding transitions; identity after each window's end.
    trans = torch.einsum("ijc,bct->btij", g.coef, e) + g.prior[symbols]
    trans = torch.where(valid[:, :, None, None], trans, g.identity)
    trans = _floor(trans)

    neg = torch.full((B, L), FLOOR, dtype=dtype, device=device)
    intron = torch.where(valid[:, None, :], e[:, 4:7, :], neg[:, None, :])    # (B, 3, L)
    donor = torch.where(valid, e[:, _CH["donor"], :], neg)                    # (B, L)
    acceptor = torch.where(valid, e[:, _CH["acceptor"], :], neg)              # (B, L)
    # Rolling m-window sums of the intron emissions, ending at t (inclusive):
    # W[b, p, t] = sum_{j = t-m+1}^{t} intron[b, p, j], defined for t >= m-1.
    if L >= m:
        window_sum = intron.unfold(2, m, 1).sum(-1)                           # (B, 3, L-m+1)
    else:
        window_sum = intron.new_full((B, 3, 0), FLOOR)
    intron_state = intron[:, g.phase, :]                                      # (B, K, L)
    window_state = window_sum[:, g.phase, :]                                  # (B, K, L-m+1)
    # Per-step slices via ``unbind``: its backward is one ``stack``, whereas
    # indexing ``trans[:, t]`` inside the loop would scatter into a full-size
    # zero gradient at every step and make the backward pass quadratic in L.
    trans_t = trans.unbind(1)
    donor_t, acceptor_t = donor.unbind(1), acceptor.unbind(1)
    intron_t = intron_state.unbind(2)
    window_t = window_state.unbind(2)

    alpha = torch.full((B, K), FLOOR, dtype=dtype, device=device)
    alpha[:, g.u_index] = 0.0
    tau = torch.full((B, K, R), FLOOR, dtype=dtype, device=device)
    pending: List[torch.Tensor] = []
    coding_mask = g.coding[None, :]
    for t in range(L):
        # 1. park donors from every S/E state at boundary t
        parked = torch.where(coding_mask, alpha + donor_t[t][:, None], torch.full_like(alpha, FLOOR))
        pending.append(parked)
        # 2. merge tail exits into the coding layer, then the coding transition
        exits = torch.logsumexp(tau + log_1mq[None], dim=-1) + acceptor_t[t][:, None]  # (B, K)
        merged = torch.logsumexp(torch.stack((alpha, exits), dim=-1), dim=-1)
        alpha = _floor(torch.logsumexp(merged[:, :, None] + trans_t[t], dim=1))
        # 3. tails continue; the donor parked m boundaries ago enters
        tau = tau + intron_t[t][:, :, None] + log_q[None]
        if len(pending) == m:
            entry = pending.pop(0)[:, :, None] + window_t[t - m + 1][:, :, None] + log_pi[None]
            tau = torch.logsumexp(torch.stack((tau, entry), dim=-1), dim=-1)
        tau = _floor(tau)
    z = alpha[:, g.u_index]
    return torch.where(z > FLOOR / 2, z, torch.full_like(z, float("-inf")))


def partition(x: str, emissions: torch.Tensor, *, code: GeneticCode = TABLES[1],
              duration: DurationMixture = DurationMixture(), tables=None) -> torch.Tensor:
    """Differentiable ``log Z`` of one window; the fast twin of
    :func:`model.a.torch_loss.partition`. ``tables`` (learned
    :class:`model.a.pooled.DurationTables`) overrides ``duration``'s values;
    ``duration`` then only fixes ``m`` and ``R``."""
    _check_input(x, emissions)
    g = Grammar.get(code, duration, dtype=emissions.dtype, device=emissions.device)
    # ``dtype=torch.long`` so a zero-length window (``log Z = 0``, as in the
    # reference) does not infer a float index tensor.
    sym = torch.tensor([symbol_indices(x)], dtype=torch.long, device=emissions.device)
    lengths = torch.tensor([len(x)], device=emissions.device)
    return log_partition_batch(emissions[None], sym, lengths, g, tables)[0]


def chain_nll(x: str, emissions: Optional[torch.Tensor],
              cds_ranges: Sequence[Range], intron_ranges: Sequence[Range], *,
              code: GeneticCode = TABLES[1], duration: DurationMixture = DurationMixture(),
              device=None, dtype=torch.float64, tables=None) -> torch.Tensor:
    """``log Z - log Z_num`` for one complete admitted chain; same contract as
    :func:`model.a.torch_loss.chain_nll`, delayed-entry recurrence. ``tables``
    as in :func:`partition`."""
    n = len(x)
    if emissions is None:
        emissions = torch.zeros((EMISSION_CHANNELS, n), device=device, dtype=dtype)
    else:
        _check_input(x, emissions)
    g = Grammar.get(code, duration, dtype=emissions.dtype, device=emissions.device)
    mask = support_mask(n, cds_ranges, intron_ranges, device=emissions.device, dtype=emissions.dtype)
    both = torch.stack((emissions, emissions + mask))                       # (2, 11, n)
    sym = torch.tensor([symbol_indices(x)] * 2, dtype=torch.long, device=emissions.device)
    lengths = torch.tensor([n, n], device=emissions.device)
    z = log_partition_batch(both, sym, lengths, g, tables)
    return z[0] - z[1]


def batch_chain_nll(windows: Sequence[str], emissions: torch.Tensor,
                    chains: Sequence[Tuple[Sequence[Range], Sequence[Range]]], *,
                    code: GeneticCode = TABLES[1],
                    duration: DurationMixture = DurationMixture(), tables=None) -> torch.Tensor:
    """Per-window ``log Z - log Z_num`` for a padded batch. ``tables`` as in
    :func:`partition`.

    ``emissions`` is ``(B, 11, L)`` with ``L >= max(len(w))``; columns past a
    window's length are ignored. ``chains[b]`` is ``(cds_ranges, intron_ranges)``
    in window coordinates. Returns ``(B,)``; a window whose support admits no
    path is ``+inf`` (drop it, section 3.6), never clamped.
    """
    B, C, L = emissions.shape
    if len(windows) != B or len(chains) != B:
        raise ValueError("windows, emissions and chains must agree on the batch size")
    device, dtype = emissions.device, emissions.dtype
    g = Grammar.get(code, duration, dtype=dtype, device=device)
    sym = torch.zeros((B, L), dtype=torch.long, device=device)
    masks = torch.full((B, C, L), float("-inf"), dtype=dtype, device=device)
    lengths = torch.tensor([len(w) for w in windows], device=device)
    for b, (w, (cds, itr)) in enumerate(zip(windows, chains)):
        n = len(w)
        if n > L:
            raise ValueError(f"window {b} is longer than the emission width {L}")
        _check_input(w, emissions[b, :, :n])
        sym[b, :n] = torch.tensor(symbol_indices(w), device=device)
        masks[b, :, :n] = support_mask(n, cds, itr, device=device, dtype=dtype)
    both = torch.cat((emissions, emissions + masks))                        # (2B, 11, L)
    z = log_partition_batch(both, torch.cat((sym, sym)), torch.cat((lengths, lengths)), g,
                            tables)
    return z[:B] - z[B:]
