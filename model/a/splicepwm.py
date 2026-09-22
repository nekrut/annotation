"""Extended splice-site scoring for the decoder (section 3.5 increment 2).

The decoder already scores a site by its context: :mod:`model.a.encoder`
emits a position-dependent ``donor`` and ``acceptor`` channel at every base,
and chromosome decoding adds :func:`model.a.pooled.motif_bias` to those rows,
so two legal GTs never had the same score. What was missing is *extended*
context. The learned dinucleotide tables of :mod:`model.a.pooled` read
exactly two bases, so once the hard mask of section 3.4 restricts decoding
to GT/GC..AG those two bases are constant across every surviving site and
that term contributes the same score to every legal donor on a chromosome,
leaving the encoder's own rows to separate true sites from false ones: chr
V's masked decode kept 1,070 true GT-AG introns but bought them with 1,970
false ones (precision 0.352). This module adds a fixed term over the wider
junction context -- a position-weight matrix over the bases around the
boundary -- on top of those existing learned scores.

Shape. A donor at boundary ``t`` (``x[t]`` is the first intron base) is
scored over offsets ``DONOR_SPAN = (-3, +6)``, i.e. the three exon bases
``x[t-3:t]`` and the six intron bases ``x[t:t+6]``; an acceptor at boundary
``t`` (``x[t-1]`` is the last intron base) over ``ACCEPTOR_SPAN = (-20, +3)``,
the twenty intron bases ``x[t-20:t]`` -- long enough to hold a polypyrimidine
tract -- and the three exon bases ``x[t:t+3]``. Both spans are half-open and
cover the dinucleotides the mask already reads, so the PWM refines the mask's
ranking rather than repeating its decision.

Estimation (:func:`fit`) counts bases per column over the *admitted reference
introns of the train split* -- the same windows :func:`model.a.train.train`
fits on, with the development sequence ids removed -- and scores a column as
the log-odds of its base against the base composition of those same windows,
with a pseudocount. Columns are treated as independent, which is the standard
weight-matrix assumption and is false at a real branch point; the resulting
score is a ranking device, not a likelihood.

Application (:func:`bias_batch`) adds the score to the ``donor``/``acceptor``
emission rows exactly as :func:`model.a.pooled.motif_bias` does, so both
decoders and both kernels stay unchanged. Two conventions are inherited from
the mask: a column whose base is not a concrete ACGT contributes 0, and a
column that reaches past the window contributes 0. A site near a window edge
is therefore scored on the part of its context that exists and is never
suppressed for lying near a seam.

Inference only. The chain-loss numerator is not touched: a reference intron
with an atypical junction must stay reachable, and fitting the encoder
against a bias estimated from the labels would double-count them.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import torch

from .pooled import _ACCEPTOR, _DONOR
from .torch_loss import EMISSION_CHANNELS

BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}
_OTHER = 4  # padded row: a non-ACGT base scores 0

# Half-open offsets relative to the boundary, in oriented window coordinates.
DONOR_SPAN = (-3, 6)
ACCEPTOR_SPAN = (-20, 3)

DEFAULT_PSEUDOCOUNT = 1.0


@dataclass(frozen=True)
class SplicePWM:
    """Log-odds columns for the donor and acceptor junction context.

    ``donor``/``acceptor`` are ``(K, 4)`` nested tuples of natural-log odds
    against ``background``, in offset order from ``donor_span[0]``. ``sites``
    records how many junctions each matrix was counted over and ``sources``
    the species and the sequence ids that were *excluded* as development, so
    a scored run can be checked against its own leakage rules.
    """

    donor: Tuple[Tuple[float, ...], ...]
    acceptor: Tuple[Tuple[float, ...], ...]
    background: Tuple[float, ...]
    donor_span: Tuple[int, int] = DONOR_SPAN
    acceptor_span: Tuple[int, int] = ACCEPTOR_SPAN
    pseudocount: float = DEFAULT_PSEUDOCOUNT
    sites: Tuple[int, int] = (0, 0)
    sources: Tuple[Dict[str, object], ...] = ()
    weighting: str = "pooled"

    def __post_init__(self):
        for name, span, mat in (("donor", self.donor_span, self.donor),
                                ("acceptor", self.acceptor_span, self.acceptor)):
            width = span[1] - span[0]
            if width <= 0:
                raise ValueError(f"{name}_span must be a non-empty half-open range")
            if len(mat) != width:
                raise ValueError(f"{name} has {len(mat)} columns, span needs {width}")
            for col in mat:
                if len(col) != 4:
                    raise ValueError(f"{name} column must have 4 entries, got {len(col)}")
        if len(self.background) != 4:
            raise ValueError("background must have 4 entries")

    # -- serialisation ----------------------------------------------------
    def to_json(self) -> dict:
        return {
            "donor": [list(c) for c in self.donor],
            "acceptor": [list(c) for c in self.acceptor],
            "background": list(self.background),
            "donor_span": list(self.donor_span),
            "acceptor_span": list(self.acceptor_span),
            "pseudocount": self.pseudocount,
            "sites": list(self.sites),
            "sources": list(self.sources),
            "weighting": self.weighting,
            "bases": BASES,
        }

    @classmethod
    def from_json(cls, d: dict) -> "SplicePWM":
        if d.get("bases", BASES) != BASES:
            raise ValueError(f"base order must be {BASES!r}")
        return cls(
            donor=tuple(tuple(float(v) for v in c) for c in d["donor"]),
            acceptor=tuple(tuple(float(v) for v in c) for c in d["acceptor"]),
            background=tuple(float(v) for v in d["background"]),
            donor_span=tuple(int(v) for v in d.get("donor_span", DONOR_SPAN)),
            acceptor_span=tuple(int(v) for v in d.get("acceptor_span", ACCEPTOR_SPAN)),
            pseudocount=float(d.get("pseudocount", DEFAULT_PSEUDOCOUNT)),
            sites=tuple(int(v) for v in d.get("sites", (0, 0))),
            sources=tuple(d.get("sources", ())),
            weighting=str(d.get("weighting", "pooled")),
        )

    @classmethod
    def load(cls, path: str) -> "SplicePWM":
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_json(json.load(fh))

    def dump(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_json(), fh, indent=2, sort_keys=True)

    # -- tables -----------------------------------------------------------
    def tables(self, *, dtype=torch.float64, device=None
               ) -> Tuple[torch.Tensor, torch.Tensor]:
        """``(K_d, 5)`` and ``(K_a, 5)`` gather tables; row 4 scores 0."""
        out = []
        for mat in (self.donor, self.acceptor):
            t = torch.zeros((len(mat), 5), dtype=dtype, device=device)
            t[:, :4] = torch.tensor(mat, dtype=dtype, device=device)
            out.append(t)
        return out[0], out[1]


# --------------------------------------------------------------------------
# estimation
# --------------------------------------------------------------------------
def _counts(width: int) -> List[List[float]]:
    return [[0.0] * 4 for _ in range(width)]


def _tally(counts: List[List[float]], window: str, t: int, span: Tuple[int, int]) -> None:
    """Add one junction's context to ``counts``; columns off the window, or on
    a base that is not a concrete ACGT, are simply not counted."""
    lo, n = span[0], len(window)
    for j in range(span[1] - span[0]):
        p = t + lo + j
        if 0 <= p < n:
            i = BASE_INDEX.get(window[p].upper())
            if i is not None:
                counts[j][i] += 1.0


def _log_odds(counts: Sequence[Sequence[float]], background: Sequence[float],
              pseudocount: float) -> Tuple[Tuple[float, ...], ...]:
    import math

    cols = []
    for col in counts:
        total = sum(col) + 4.0 * pseudocount
        if total <= 0:
            # Nothing was counted in this column and there is no pseudocount
            # to fall back on: score it flat rather than invent a preference.
            cols.append((0.0,) * 4)
            continue
        # A base that was never seen scores ``-inf`` only when the caller
        # asked for no pseudocount; the default smooths it.
        cols.append(tuple(
            (math.log(((c + pseudocount) / total) / background[i])
             if c + pseudocount > 0 else float("-inf"))
            for i, c in enumerate(col)))
    return tuple(cols)


def _tally_examples(examples: Sequence, donor_span: Tuple[int, int],
                    acceptor_span: Tuple[int, int]):
    """Column counts, ACGT background counts and junction count over windows.

    One ``intron_ranges`` entry ``(a, b)`` contributes one donor at ``a`` and
    one acceptor at ``b``, so the two matrices are always counted over the
    same junctions and the returned count applies to both.
    """
    d_counts = _counts(donor_span[1] - donor_span[0])
    a_counts = _counts(acceptor_span[1] - acceptor_span[0])
    bg = [0.0] * 4
    n = 0
    for ex in examples:
        w = ex.window
        for ch in w:
            i = BASE_INDEX.get(ch.upper())
            if i is not None:
                bg[i] += 1.0
        for a, b in ex.intron_ranges:
            _tally(d_counts, w, a, donor_span)
            _tally(a_counts, w, b, acceptor_span)
            n += 1
    return d_counts, a_counts, bg, n


def fit(examples: Sequence, *, donor_span: Tuple[int, int] = DONOR_SPAN,
        acceptor_span: Tuple[int, int] = ACCEPTOR_SPAN,
        pseudocount: float = DEFAULT_PSEUDOCOUNT,
        sources: Sequence[Dict[str, object]] = (),
        min_sites: int = 1) -> SplicePWM:
    """Estimate a :class:`SplicePWM` from ``WindowExample``s.

    Every ``intron_ranges`` entry ``(a, b)`` of every example contributes one
    donor at boundary ``a`` and one acceptor at boundary ``b``. The background
    is the ACGT composition of the same windows, so a column that carries no
    information scores about 0 and the matrix adds nothing there.

    Pass only the *train* split: see :func:`model.a.train.split_windows`.
    ``min_sites`` refuses a matrix counted over fewer junctions than that,
    because a handful of sites gives a PWM that is mostly pseudocount.
    """
    d_counts, a_counts, bg, n_sites = _tally_examples(examples, donor_span, acceptor_span)
    n_donor = n_acceptor = n_sites
    if n_donor < min_sites or n_acceptor < min_sites:
        raise ValueError(f"only {n_donor} donor and {n_acceptor} acceptor sites, "
                         f"min_sites={min_sites}")
    total = sum(bg)
    if total <= 0:
        raise ValueError("no concrete ACGT base in the given windows")
    background = tuple(c / total for c in bg)
    if min(background) <= 0:
        raise ValueError(f"a base is absent from the background: {background}")
    return SplicePWM(
        donor=_log_odds(d_counts, background, pseudocount),
        acceptor=_log_odds(a_counts, background, pseudocount),
        background=background,
        donor_span=tuple(donor_span), acceptor_span=tuple(acceptor_span),
        pseudocount=pseudocount, sites=(n_donor, n_acceptor),
        sources=tuple(sources),
    )


def _mean_frequencies(per_group: Sequence[Sequence[Sequence[float]]],
                      pseudocount: float) -> List[List[float]]:
    """Column frequencies averaged over groups, each group weighted 1/G.

    A group's column is smoothed with ``pseudocount`` and normalised *within
    the group* before the average, which is what makes 2 yeast junctions
    count as much as 60,000 nematode ones. A group that counted nothing in a
    column and was given no pseudocount contributes zero there rather than a
    uniform guess, so an empty group cannot flatten a column the others
    resolved.
    """
    width = len(per_group[0])
    g = float(len(per_group))
    out = []
    for j in range(width):
        acc = [0.0] * 4
        for counts in per_group:
            col = counts[j]
            total = sum(col) + 4.0 * pseudocount
            if total <= 0:
                continue
            for i, c in enumerate(col):
                acc[i] += ((c + pseudocount) / total) / g
        out.append(acc)
    return out


def _log_odds_freq(freqs: Sequence[Sequence[float]],
                   background: Sequence[float]) -> Tuple[Tuple[float, ...], ...]:
    """Natural-log odds of already-normalised column frequencies."""
    import math

    cols = []
    for col in freqs:
        if sum(col) <= 0:
            cols.append((0.0,) * 4)
            continue
        cols.append(tuple((math.log(f / background[i]) if f > 0 else float("-inf"))
                          for i, f in enumerate(col)))
    return tuple(cols)


def fit_grouped(groups: Sequence[Tuple[str, Sequence]], *,
                donor_span: Tuple[int, int] = DONOR_SPAN,
                acceptor_span: Tuple[int, int] = ACCEPTOR_SPAN,
                pseudocount: float = DEFAULT_PSEUDOCOUNT,
                sources: Sequence[Dict[str, object]] = (),
                min_sites: int = 1,
                min_group_sites: int = 1) -> SplicePWM:
    """Estimate a source-balanced :class:`SplicePWM` from grouped windows.

    ``groups`` is a sequence of ``(name, examples)``, one per source, each
    holding that source's *train* windows only. Every source contributes the
    same total weight: a column's probability is the unweighted mean of the
    per-source smoothed frequencies, and the background is the unweighted
    mean of the per-source ACGT compositions, so neither the matrix nor its
    null model is decided by whichever species has the most introns.

    This is the section 3.5 answer to the pooled matrix's species-independence
    failure, and it buys that at a cost worth stating: a source with a handful
    of train junctions now carries a full share of the score, so its
    pseudocount-dominated columns are pulled into every other species' decode.
    ``min_group_sites`` is the guard; it is not a substitute for judgement
    about which sources belong in the pool.

    ``min_sites`` applies to the total junction count, as in :func:`fit`.
    Site counts land in ``sites`` and, per source, in the ``train_sites``
    key of the matching ``sources`` entry.
    """
    if not groups:
        raise ValueError("fit_grouped needs at least one group")
    per_donor, per_acceptor, per_bg, counts_by_name = [], [], [], {}
    total = 0
    for name, examples in groups:
        d_counts, a_counts, bg, n = _tally_examples(examples, donor_span, acceptor_span)
        if n < min_group_sites:
            raise ValueError(f"source {name!r} has {n} train junctions, "
                             f"min_group_sites={min_group_sites}")
        bg_total = sum(bg)
        if bg_total <= 0:
            raise ValueError(f"source {name!r} has no concrete ACGT base")
        per_donor.append(d_counts)
        per_acceptor.append(a_counts)
        per_bg.append([c / bg_total for c in bg])
        counts_by_name[name] = n
        total += n
    if total < min_sites:
        raise ValueError(f"only {total} donor and {total} acceptor sites, "
                         f"min_sites={min_sites}")
    g = float(len(groups))
    background = tuple(sum(b[i] for b in per_bg) / g for i in range(4))
    if min(background) <= 0:
        raise ValueError(f"a base is absent from the background: {background}")
    annotated = []
    for src in sources:
        entry = dict(src)
        if entry.get("name") in counts_by_name:
            entry["train_sites"] = counts_by_name[entry["name"]]
        annotated.append(entry)
    return SplicePWM(
        donor=_log_odds_freq(_mean_frequencies(per_donor, pseudocount), background),
        acceptor=_log_odds_freq(_mean_frequencies(per_acceptor, pseudocount), background),
        background=background,
        donor_span=tuple(donor_span), acceptor_span=tuple(acceptor_span),
        pseudocount=pseudocount, sites=(total, total),
        sources=tuple(annotated), weighting="balanced",
    )


# --------------------------------------------------------------------------
# application
# --------------------------------------------------------------------------
def _base_index_tensor(window: str) -> torch.Tensor:
    return torch.tensor([BASE_INDEX.get(c.upper(), _OTHER) for c in window],
                        dtype=torch.long)


def _score_row(idx: torch.Tensor, table: torch.Tensor, span: Tuple[int, int],
               n: int) -> torch.Tensor:
    """``(L,)`` PWM score per boundary, ``idx`` the ``(L,)`` base indices of a
    window of length ``n`` padded to ``L`` with ``_OTHER``."""
    L = idx.shape[0]
    out = torch.zeros(L, dtype=table.dtype, device=table.device)
    lo = span[0]
    for j in range(span[1] - span[0]):
        off = lo + j
        # boundary t reads idx[t + off]; keep only t with 0 <= t+off < n.
        t0, t1 = max(0, -off), min(n, n - off)
        if t1 <= t0:
            continue
        out[t0:t1] += table[j][idx[t0 + off:t1 + off]]
    return out


def bias_batch(windows: Sequence[str], pwm: SplicePWM, *, L: Optional[int] = None,
               dtype=torch.float64, device=None) -> torch.Tensor:
    """``(B, EMISSION_CHANNELS, L)`` additive bias: the donor and acceptor PWM
    scores on their emission rows, zero elsewhere and past each window's end.

    Same contract as :func:`model.a.pooled.motif_bias_batch`, and meant to be
    added to it: the dinucleotide table stays the learned part of the score
    and this is the fixed context term on top of it.
    """
    B = len(windows)
    L = max((len(w) for w in windows), default=0) if L is None else L
    d_tab, a_tab = pwm.tables(dtype=dtype, device=device)
    bias = torch.zeros((B, EMISSION_CHANNELS, L), dtype=dtype, device=device)
    for b, w in enumerate(windows):
        n = len(w)
        if n > L:
            raise ValueError(f"window {b} is longer than L={L}")
        if not n:
            continue
        idx = torch.full((L,), _OTHER, dtype=torch.long)
        idx[:n] = _base_index_tensor(w)
        idx = idx.to(device)
        bias[b, _DONOR, :] = _score_row(idx, d_tab, pwm.donor_span, n)
        bias[b, _ACCEPTOR, :] = _score_row(idx, a_tab, pwm.acceptor_span, n)
    return bias


def bias(x: str, pwm: SplicePWM, *, dtype=torch.float64, device=None) -> torch.Tensor:
    """``(EMISSION_CHANNELS, n)`` :func:`bias_batch` for one window."""
    return bias_batch([x], pwm, dtype=dtype, device=device)[0, :, :len(x)]
