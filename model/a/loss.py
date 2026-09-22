"""Chain loss for candidate A: the duration-mixture CRF negative log-likelihood
of an admitted gold chain (proposal sections 3.1-3.2, admission per 3.6).

The loss of one admitted representative is

    L(x, e) = log Z(x, e)  -  log Z_num(x, e; chain)

where ``e`` are the encoder's eleven emission channels (``grammar.Scores``),
``Z`` the free partition over every legal path of A's grammar and ``Z_num`` the
partition restricted to the gold chain's *support*. The restriction is a hard
``-inf`` mask on every emission a gold path cannot use (:func:`numerator_scores`
below): U only outside the gene, a coding channel only on a CDS base, an intron
channel only on an intron base, and ``start``/``stop``/``donor``/``acceptor``
only at their gold coordinates. The grammar's own state machine then fixes the
reading frame and the junction structure, so ``Z_num`` marginalizes exactly the
two things the chain leaves latent -- the per-base phase channel (only one phase
is reachable per position given the fixed start) and the intron duration
mixture components -- and nothing else. Hence ``L >= 0`` with equality only when
the gold chain already carries all of the partition mass.

This is the standard-library *oracle* the Phase-4 PyTorch loss must reproduce on
the section 3.4 fixtures: the torch path builds the identical support mask and
computes both partitions with a ``logsumexp`` forward pass whose autograd yields
``dL/de = posterior_free - posterior_num`` (the CRF marginal difference). Keeping
the oracle in the standard library lets it run without torch on any host, exactly
as ``tests/test_a_encoder.py`` gates its torch checks. Standard library only.
"""
from __future__ import annotations

from copy import deepcopy
from math import inf
from typing import List, Optional, Sequence, Tuple

from model.grammar import ReferenceDecoder, Scores

NEG = -inf
Range = Tuple[int, int]


def _mask_all(sc: Scores, t: int) -> None:
    """Forbid every emission at boundary ``t`` (caller re-enables the allowed
    ones from the untouched base scores)."""
    sc.u[t] = NEG
    for p in range(3):
        sc.cds[p][t] = NEG
        sc.intron[p][t] = NEG
    sc.start[t] = NEG
    sc.stop[t] = NEG
    sc.donor[t] = NEG
    sc.acceptor[t] = NEG


def numerator_scores(base: Optional[Scores], n: int, cds_ranges: Sequence[Range],
                     intron_ranges: Sequence[Range]) -> Scores:
    """Hard support mask for one *complete* admitted chain, or for a gene-free
    window when both range lists are empty (every base intergenic ``U``).

    ``cds_ranges`` and ``intron_ranges`` are the merged, oriented, half-open CDS
    and intron intervals of the chain in a window of length ``n`` (the same
    convention as ``model.labels.numerator_check.oriented_chain``); everything
    outside the CDS/intron cover is intergenic ``U``. Returns a copy of ``base``
    (zeros when ``base`` is ``None``) in which every emission a gold path cannot
    use is ``-inf`` and every emission it can use keeps its base score.

    The auxiliary channels sit at the boundary that *consumes* the base carrying
    them, matching ``ReferenceDecoder.transitions``: ``start`` on the first CDS
    base, ``stop`` on the last CDS base (the base completing the terminal codon),
    ``donor`` on each intron's first base, ``acceptor`` on the first CDS base
    after each intron.
    """
    if not cds_ranges and intron_ranges:
        raise ValueError("introns without a CDS interval are not a chain")
    # Empty ranges: the gene-free (background) window, every base intergenic U.
    sc = deepcopy(base) if base is not None else Scores.zeros(n)
    if sc.n != n:
        raise ValueError("base scores length differs from window length")
    keep = deepcopy(sc)                       # untouched copy of the allowed values

    cds = sorted(cds_ranges)
    introns = sorted(intron_ranges)
    cds_bases = {t for a, b in cds for t in range(a, b)}
    intron_bases = {t for a, b in introns for t in range(a, b)}
    if cds_bases & intron_bases:
        raise ValueError("CDS and intron intervals overlap")

    donors = {a for a, _ in introns}          # first intronic base
    acceptors = {b for _, b in introns}       # first CDS base after the intron
    start_t = cds[0][0] if cds else -1         # first CDS base
    stop_t = cds[-1][1] - 1 if cds else -1     # last CDS base

    for t in range(n):
        _mask_all(sc, t)
        if t in cds_bases:
            for p in range(3):
                sc.cds[p][t] = keep.cds[p][t]
            if t == start_t:
                sc.start[t] = keep.start[t]
            if t == stop_t:
                sc.stop[t] = keep.stop[t]
            if t in acceptors:
                sc.acceptor[t] = keep.acceptor[t]
        elif t in intron_bases:
            for p in range(3):
                sc.intron[p][t] = keep.intron[p][t]
            if t in donors:
                sc.donor[t] = keep.donor[t]
        else:                                  # intergenic
            sc.u[t] = keep.u[t]
    return sc


def chain_nll(decoder: ReferenceDecoder, x: str, base: Optional[Scores],
              cds_ranges: Sequence[Range], intron_ranges: Sequence[Range]):
    """Return ``(loss, log_Z, log_Z_num)`` for one *complete* admitted gold chain.

    ``loss = log_Z - log_Z_num`` is the per-representative chain NLL. ``decoder``
    supplies the grammar, genetic code and duration mixture. A non-finite
    ``log_Z_num`` means the support admits no legal path -- the adapter must
    record the transcript id and drop the crop (section 3.6), never clamp an
    infinite loss.

    This increment covers only complete targets: the gold path enters and leaves
    in intergenic ``U``, so ``decoder`` must have its edge grammar disabled
    (``edges is None``). The support mask restricts the *emission* channels but
    not the initial/terminal *states*, so an enabled ``EdgePrior`` would let the
    numerator claim extra boundary/phase hypotheses (for example ``log Z_num`` of
    a complete gene at a real edge becomes ``log(5/4)`` instead of ``0``), quietly
    supervising the wrong path set. Declared edge-partial targets need the
    section-3.6 boundary-support interface -- compatible entry/exit families and a
    first-row phase carried alongside ``cds_ranges`` -- which does not exist yet,
    so an edge-enabled decoder is rejected rather than silently mis-scored.
    """
    if getattr(decoder, "edges", None) is not None:
        raise ValueError(
            "chain_nll is the complete-target oracle and requires the decoder's "
            "edge grammar disabled (edges is None); edge-partial numerators are a "
            "later section-3.6 increment")
    num = numerator_scores(base, len(x), cds_ranges, intron_ranges)
    log_z = decoder.partition(x, base if base is not None else Scores.zeros(len(x)))
    log_z_num = decoder.partition(x, num)
    return log_z - log_z_num, log_z, log_z_num
