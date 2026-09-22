"""Differentiable PyTorch chain loss for candidate A (proposal sections 3.1-3.2,
admission per 3.6).

This is the torch counterpart of the standard-library oracle in
:mod:`model.a.loss`. The loss of one admitted *complete* gold chain is

    L(x, e) = log Z(x, e)  -  log Z_num(x, e; chain)

with ``e`` the encoder's eleven per-base emission channels, ``Z`` the free
partition over every legal path of A's grammar and ``Z_num`` the partition
restricted to the gold chain's hard emission support (:func:`support_mask`).
Because ``Z`` and ``Z_num`` are computed with :func:`torch.logsumexp`, autograd
yields exactly the CRF marginal difference ``dL/de = P_free(e) - P_num(e)`` --
the gradient a training loop applies to the emission head.

**Parity by construction.** Rather than re-derive the recurrence, the forward
here reuses the *reviewed* grammar state machine unchanged:
:meth:`ReferenceDecoder.initial`, :meth:`ReferenceDecoder.transitions` and
:meth:`ReferenceDecoder.terminal`. Those methods only ever *read* emission
channels as ``sc.u[t]``, ``sc.cds[p][t]`` and so on and add Python-float grammar
priors (uniform base priors, duration ``log_pi``/``log_q``/``log_1mq``); feeding
them a :class:`TorchScores` view over a ``(11, n)`` torch tensor makes every
yielded transition score a torch tensor, so the only substitution against the
standard-library oracle is Python ``float`` arithmetic + ``math`` ``logsumexp``
for torch arithmetic + :func:`torch.logsumexp`. The value matches
``model.a.loss.chain_nll`` to floating-point tolerance on the section-3.4
fixtures (``tests/test_a_torch_loss.py``), exactly as ``DelayedEntryDecoder`` is
checked against ``ReferenceDecoder``.

**Scope and cost.** This forward has the *reference* recurrence's cost -- an
intron carries its ``m-1`` mandatory positions as explicit states and each
boundary is a Python dict, so with the default ``m=20`` it is far too slow for a
training window. It is the differentiable *reference* the fast vectorized
kernel (a delayed-entry forward, the next T-human-014 increment) is checked
against, and it drives the small-window parity and gradient tests today. Like
``model.a.loss.chain_nll`` it covers only *complete* targets: an edge-enabled
decoder is rejected, because the emission mask constrains channels, not the
boundary states an :class:`EdgePrior` would open. Edge-partial numerators are
the same section-3.6 boundary-support increment named there.
"""
from __future__ import annotations

from typing import Optional, Sequence, Tuple

import torch

from model.grammar import DurationMixture, GeneticCode, ReferenceDecoder, Scores, TABLES

Range = Tuple[int, int]

# Emission channel order, identical to ``model.grammar.Scores.channels()`` and to
# the encoder's ``[B, 11, L]`` head (``model.a.encoder``). This is the contract
# that binds an emission tensor row to a grammar channel.
CHANNEL_ORDER = ("u", "cds[0]", "cds[1]", "cds[2]",
                 "intron[0]", "intron[1]", "intron[2]",
                 "start", "stop", "donor", "acceptor")
EMISSION_CHANNELS = len(CHANNEL_ORDER)


class TorchScores:
    """A ``Scores``-shaped read-only view over a ``(11, n)`` torch tensor.

    Exposes exactly the attributes the grammar transition generator reads --
    ``u``, ``cds``, ``intron`` (lists of three phase rows), ``start``, ``stop``,
    ``donor`` and ``acceptor`` -- as tensor rows, so ``sc.cds[p][t]`` is a 0-dim
    tensor that shares storage and autograd history with the emission tensor.
    Rows in ``CHANNEL_ORDER``.
    """

    def __init__(self, emissions: torch.Tensor):
        if emissions.dim() != 2 or emissions.shape[0] != EMISSION_CHANNELS:
            raise ValueError(
                f"emissions must be ({EMISSION_CHANNELS}, n); got {tuple(emissions.shape)}")
        self.n = emissions.shape[1]
        self._e = emissions
        self.u = emissions[0]
        self.cds = [emissions[1], emissions[2], emissions[3]]
        self.intron = [emissions[4], emissions[5], emissions[6]]
        self.start = emissions[7]
        self.stop = emissions[8]
        self.donor = emissions[9]
        self.acceptor = emissions[10]


def support_mask(n: int, cds_ranges: Sequence[Range], intron_ranges: Sequence[Range],
                 *, device=None, dtype=torch.float64) -> torch.Tensor:
    """Additive ``(11, n)`` support mask for one *complete* admitted chain,
    or for a gene-free window when both range lists are empty.

    Every emission a gold path cannot use is ``-inf`` and every emission it can
    use is ``0``; adding this to an emission tensor keeps the allowed scores and
    drives the rest to ``-inf``, exactly as :func:`model.a.loss.numerator_scores`
    copies base values where allowed and masks the rest. ``cds_ranges`` and
    ``intron_ranges`` are merged, oriented, half-open intervals (the
    ``model.labels.numerator_check.oriented_chain`` convention); everything
    outside their cover is intergenic ``U``. The auxiliary channels sit at the
    base that *consumes* them: ``start`` on the first CDS base, ``stop`` on the
    last, ``donor`` on each intron's first base, ``acceptor`` on the first CDS
    base after each intron.
    """
    if not cds_ranges and intron_ranges:
        raise ValueError("introns without a CDS interval are not a chain")
    # An empty chain is the gene-free (background) window: every base is
    # intergenic ``U`` (proposal 3, "gene-free/background windows").
    cds = sorted(cds_ranges)
    introns = sorted(intron_ranges)
    cds_bases = {t for a, b in cds for t in range(a, b)}
    intron_bases = {t for a, b in introns for t in range(a, b)}
    if cds_bases & intron_bases:
        raise ValueError("CDS and intron intervals overlap")
    donors = {a for a, _ in introns}
    acceptors = {b for _, b in introns}
    start_t = cds[0][0] if cds else -1
    stop_t = cds[-1][1] - 1 if cds else -1

    neg = float("-inf")
    mask = torch.full((EMISSION_CHANNELS, n), neg, device=device, dtype=dtype)
    for t in range(n):
        if t in cds_bases:
            mask[1, t] = mask[2, t] = mask[3, t] = 0.0     # cds[0..2]
            if t == start_t:
                mask[7, t] = 0.0
            if t == stop_t:
                mask[8, t] = 0.0
            if t in acceptors:
                mask[10, t] = 0.0
        elif t in intron_bases:
            mask[4, t] = mask[5, t] = mask[6, t] = 0.0     # intron[0..2]
            if t in donors:
                mask[9, t] = 0.0
        else:
            mask[0, t] = 0.0                                # intergenic U
    return mask


def _check_input(x: str, emissions: torch.Tensor) -> None:
    """Enforce the scalar oracle's finite-or-``-inf`` emission contract.

    The torch forward reuses the grammar recurrence, which -- like
    :meth:`ReferenceDecoder._check_input` -- assumes every emission channel is
    finite or ``-inf`` (a hard-masked transition). ``_partition`` treats
    ``isinf`` as "drop this path", so a stray ``+inf`` would silently discard a
    legal transition and a ``NaN`` would poison ``logsumexp``; neither public
    entry point may accept them. Legitimate ``-inf`` support is preserved. The
    emission width must also equal ``len(x)``, since the recurrence iterates over
    the emission columns and would otherwise truncate or over-run the sequence.
    """
    if emissions.dim() != 2 or emissions.shape != (EMISSION_CHANNELS, len(x)):
        raise ValueError(
            f"emissions must be ({EMISSION_CHANNELS}, {len(x)}); "
            f"got {tuple(emissions.shape)}")
    if torch.isnan(emissions).any() or (emissions == float("inf")).any():
        raise ValueError(
            "emission channels must be finite or -inf; got NaN or +inf")


def _partition(decoder: ReferenceDecoder, x: str, sc: TorchScores,
               ref: torch.Tensor) -> torch.Tensor:
    """log Z over every legal path 0->n of the reference grammar, on torch
    emissions. ``ref`` is any tensor from the emission graph, used to source the
    device/dtype of the Python-float grammar constants."""
    def const(v: float) -> torch.Tensor:
        return ref.new_tensor(v)

    n = sc.n
    # Boundary-0 states; with edges disabled this is just {U: 0}.
    layer = {state: const(w) for state, w in decoder.initial().items()}
    for t in range(n):
        cur: dict = {}
        for state, score in layer.items():
            if torch.isinf(score):
                continue
            for nxt, s, _tag in decoder.transitions(state, x, t, sc):
                tot = score + s                       # torch + torch (or + float)
                if torch.isinf(tot):
                    continue                          # a -inf-masked emission: contributes 0
                cur.setdefault(nxt, []).append(tot)
        layer = {k: torch.logsumexp(torch.stack(v), dim=0) for k, v in cur.items()}
    finals = []
    for state, v in layer.items():
        term = decoder.terminal(state)
        if term == float("-inf"):
            continue
        finals.append(v + term)
    if not finals:
        return const(float("-inf"))
    return torch.logsumexp(torch.stack(finals), dim=0)


def partition(x: str, emissions: torch.Tensor, *, code: GeneticCode = TABLES[1],
              duration: DurationMixture = DurationMixture()) -> torch.Tensor:
    """Differentiable log Z of the free grammar over ``emissions`` ``(11, n)``."""
    _check_input(x, emissions)
    decoder = ReferenceDecoder(code, duration)
    return _partition(decoder, x, TorchScores(emissions), emissions)


def chain_nll(x: str, emissions: Optional[torch.Tensor],
              cds_ranges: Sequence[Range], intron_ranges: Sequence[Range], *,
              code: GeneticCode = TABLES[1], duration: DurationMixture = DurationMixture(),
              device=None, dtype=torch.float64) -> torch.Tensor:
    """``log Z - log Z_num`` for one *complete* admitted gold chain.

    ``emissions`` is a ``(11, n)`` torch tensor in :data:`CHANNEL_ORDER`; ``None``
    means zeros (the pure grammar-plus-duration score, matching the section-3.4
    cases). The returned scalar is differentiable w.r.t. ``emissions``: its
    gradient is the free-minus-numerator emission posterior. A non-finite
    numerator partition means the support admits no legal path -- the adapter
    must drop that crop (section 3.6), never clamp an infinite loss.

    Complete-target oracle only: an edge-enabled decoder is rejected for the
    same reason as :func:`model.a.loss.chain_nll` -- the emission mask does not
    constrain the boundary states an :class:`~model.grammar.EdgePrior` opens, so
    the numerator would claim extra entry/exit/phase hypotheses.
    """
    decoder = ReferenceDecoder(code, duration)  # edges default None
    if getattr(decoder, "edges", None) is not None:  # defensive: future default change
        raise ValueError("chain_nll is the complete-target oracle; edges must be disabled")
    n = len(x)
    if emissions is None:
        emissions = torch.zeros((EMISSION_CHANNELS, n), device=device, dtype=dtype)
    else:
        _check_input(x, emissions)

    log_z = _partition(decoder, x, TorchScores(emissions), emissions)
    mask = support_mask(n, cds_ranges, intron_ranges,
                        device=emissions.device, dtype=emissions.dtype)
    log_z_num = _partition(decoder, x, TorchScores(emissions + mask), emissions)
    return log_z - log_z_num
