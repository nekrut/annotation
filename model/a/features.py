"""Eight-channel per-base input featurizer for candidate A (proposal section 3.5).

Channels, in order:

  0..3  ACGT indicators (one-hot; all zero for an ambiguous or absent base)
  4     ambiguity      (1 for an IUPAC ambiguity code / N over real sequence)
  5     soft masking    (1 for a soft-masked lowercase base)
  6     local GC        (fraction of G/C among *unambiguous* bases in a
                         centred 129-base window; a fixed 0.0 if none exist)
  7     real-sequence availability (1 where padding did not invent the base)

None of these features use a reference annotation. There are no trainable
parameters here, so the featurizer does not enter the section 3.5 count.

The function is pure Python plus torch only for the returned tensor, so its
logic (especially the GC window) is checkable without a GPU.
"""
from __future__ import annotations

from typing import Optional, Sequence

GC_WINDOW = 129  # centred window for the local GC channel (section 3.5)

_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _classify(seq: str):
    """Return per-position (base_index|None, is_ambiguous, is_softmasked, is_real)."""
    out = []
    for ch in seq:
        is_soft = ch.islower()
        up = ch.upper()
        if up in _BASE_INDEX:
            out.append((_BASE_INDEX[up], False, is_soft, True))
        else:  # N or any IUPAC ambiguity code over real sequence
            out.append((None, True, is_soft, True))
    return out


def gc_track(
    seq: str,
    window: int = GC_WINDOW,
    available: Optional[Sequence[bool]] = None,
) -> list:
    """Local GC fraction over a centred window, ignoring ambiguous bases.

    A position whose window contains no unambiguous base gets 0.0, per the
    section 3.5 "fixed zero value if none exist" rule. O(len(seq)) via a
    rolling count so it is cheap enough to run on whole chromosomes.

    `available[i] == False` marks a padded position that invented no base; it
    is excluded from both GC counts so a padded neighbour never contaminates a
    real position's window. When `available` is None every position counts.
    """
    n = len(seq)
    if available is None:
        available = [True] * n
    half = window // 2
    # Prefix sums of (gc, unambiguous) so any window is an O(1) difference.
    gc_ps = [0] * (n + 1)
    un_ps = [0] * (n + 1)
    for i, ch in enumerate(seq):
        up = ch.upper()
        unamb = available[i] and up in _BASE_INDEX
        is_gc = available[i] and up in ("G", "C")
        gc_ps[i + 1] = gc_ps[i] + (1 if is_gc else 0)
        un_ps[i + 1] = un_ps[i] + (1 if unamb else 0)
    track = [0.0] * n
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        unamb = un_ps[hi] - un_ps[lo]
        if unamb:
            track[i] = (gc_ps[hi] - gc_ps[lo]) / unamb
    return track


def encode_sequence(seq: str, available: Optional[Sequence[bool]] = None):
    """Build the [8, L] feature tensor for one oriented sequence.

    `available[i] == False` marks a padded position (channel 7 = 0 and every
    other channel forced to 0); padding invents no sequence. When `available`
    is None every position is treated as real.

    Requires torch; import is deferred so the module docstring/logic and
    `gc_track` are usable without it.
    """
    import torch  # local import: keep the module importable without torch

    n = len(seq)
    if available is None:
        available = [True] * n
    elif len(available) != n:
        raise ValueError("available mask length must match sequence length")

    feats = torch.zeros(8, n, dtype=torch.float32)
    classes = _classify(seq)
    gc = gc_track(seq, available=available)
    for i, (bi, is_amb, is_soft, is_real) in enumerate(classes):
        if not available[i]:
            continue  # padded: all channels stay 0, availability stays 0
        if bi is not None:
            feats[bi, i] = 1.0
        if is_amb:
            feats[4, i] = 1.0
        if is_soft:
            feats[5, i] = 1.0
        feats[6, i] = gc[i]
        feats[7, i] = 1.0
    return feats
