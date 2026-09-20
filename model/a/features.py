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

`gc_track` and `_classify` are the pure-Python reference of the channel
logic (checkable without torch); `encode_sequence` is the vectorised torch
implementation of the same rules over the byte codes of the sequence
(:func:`base_codes`), tested for exact equality against the reference. On one
core the reference costs ~3.7 CPU-s per oriented Mb and the vectorised path
well under 0.1 (a-pilot 3.2 item 6, revision step 1).
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


def base_codes(seq: str):
    """Byte-level classification of ``seq`` for the vectorised paths.

    Returns ``(index, lower)``: ``index`` is an int64 tensor of length
    ``len(seq)`` with 0..3 for A/C/G/T in either case and 4 for any other
    character (N, IUPAC ambiguity codes, anything else, as in `_classify`);
    ``lower`` is a bool tensor marking ASCII lower-case letters (soft masking,
    the same predicate as ``str.islower`` on one ASCII character).
    """
    import torch

    try:
        raw = seq.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError("sequence must be ASCII") from exc
    codes = torch.frombuffer(bytearray(raw), dtype=torch.uint8) if raw else torch.zeros(0, dtype=torch.uint8)
    lower = (codes >= 97) & (codes <= 122)
    return _BASE_TABLE[codes.long()], lower


def _build_base_table():
    """256-entry lookup: byte -> base index (0..3) or 4."""
    import torch

    table = torch.full((256,), 4, dtype=torch.long)
    for ch, idx in _BASE_INDEX.items():
        table[ord(ch)] = idx
        table[ord(ch.lower())] = idx
    return table


try:
    _BASE_TABLE = _build_base_table()
except ImportError:  # torch absent: only the pure-Python reference is usable
    _BASE_TABLE = None


def _bool_tensor(mask) -> "torch.Tensor":
    """A bool tensor from a sequence of bools (or a tensor), via the buffer
    protocol: ~10x cheaper than ``torch.tensor(list)`` for a 12 kb window."""
    import torch

    if isinstance(mask, torch.Tensor):
        return mask.to(torch.bool)
    if not len(mask):
        return torch.zeros(0, dtype=torch.bool)
    return torch.frombuffer(bytearray(bytes(bool(v) for v in mask)), dtype=torch.bool)


def encode_sequence(seq: str, available: Optional[Sequence[bool]] = None):
    """Build the [8, L] feature tensor for one oriented sequence.

    `available[i] == False` marks a padded position (channel 7 = 0 and every
    other channel forced to 0); padding invents no sequence. When `available`
    is None every position is treated as real.

    Vectorised torch implementation of `_classify` + `gc_track`; the result is
    bit-identical to `encode_sequence_reference` (the GC fraction is formed in
    float64 and rounded to float32 once, as the reference does).

    Requires torch; import is deferred so the module docstring/logic and
    `gc_track` are usable without it.
    """
    import torch  # local import: keep the module importable without torch

    n = len(seq)
    if available is None:
        avail = torch.ones(n, dtype=torch.bool)
    else:
        if len(available) != n:
            raise ValueError("available mask length must match sequence length")
        avail = _bool_tensor(available)
    index, lower = base_codes(seq)
    unamb = (index < 4) & avail
    feats = torch.zeros(8, n, dtype=torch.float32)
    if n == 0:
        return feats
    pos = torch.arange(n)
    feats[index[unamb], pos[unamb]] = 1.0
    feats[4] = ((index == 4) & avail).float()
    feats[5] = (lower & avail).float()
    # Local GC over a centred window, ambiguous and unavailable bases excluded
    # (prefix sums, as in gc_track).
    is_gc = ((index == 1) | (index == 2)) & avail
    gc_ps = torch.zeros(n + 1, dtype=torch.int64)
    un_ps = torch.zeros(n + 1, dtype=torch.int64)
    torch.cumsum(is_gc.long(), 0, out=gc_ps[1:])
    torch.cumsum(unamb.long(), 0, out=un_ps[1:])
    half = GC_WINDOW // 2
    lo = (pos - half).clamp_(min=0)
    hi = (pos + half + 1).clamp_(max=n)
    un_w = un_ps[hi] - un_ps[lo]
    gc_w = gc_ps[hi] - gc_ps[lo]
    track = torch.where(un_w > 0, gc_w.double() / un_w.clamp(min=1).double(), torch.zeros((), dtype=torch.float64))
    feats[6] = track.float() * avail
    feats[7] = avail.float()
    return feats


def encode_sequence_reference(seq: str, available: Optional[Sequence[bool]] = None):
    """Per-base Python reference for `encode_sequence` (same contract); kept
    for the parity test and for checking the channel rules without reading
    the vectorised code."""
    import torch

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
