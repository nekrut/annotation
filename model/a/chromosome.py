"""End-to-end chromosome inference for candidate A: both strands, overlapping
full-sequence windows, and GFF3 output (the section-5 row of
``docs/design/a-pilot.md``, as opposed to the annotation-selected window
profile of ``measure --profile windows``).

Every base of the sequence is processed on both strands. Each oriented strand
(the plus strand as read, the minus strand as its reverse complement) is cut
into windows of ``window`` bases that overlap by ``overlap`` bases
(:func:`tiles`); each window is featurized, run through the encoder, given
the pooled decoder's dinucleotide bias, and decoded with the batched tensor
Viterbi (:mod:`model.a.fast_viterbi`). The free grammar has no edge prior, so
a window only yields *complete* chains; a gene no longer than half
``overlap`` is therefore contained in the window that owns it, and a longer
one may be lost or truncated at this stage (the section-3.6 boundary-support
increment is where partial chains and seam-crossing genes belong). Chains
that two overlapping windows both decode are reported once: a chain belongs
to the window whose *core* (:func:`tiles`: the window minus half the overlap
on each interior side) contains its oriented 5' start. Containment is the
guarantee; the *context* around an owned chain is at least half the overlap
before its start and, in the interior, ``overlap - overlap // 2`` after its
start -- so a chain of the maximal length reaches the window edge (see
:func:`tiles`), and a chain near a true sequence end has only the sequence
that exists.

The encoder's pooling grid is anchored to the oriented chromosome origin
(proposal section 3.5), not to each window: a window's encoder input starts
at the stride multiple at or before ``tile.start`` and its emissions are
cropped back to the tile, so an interior base is pooled with the same
neighbours whichever window encodes it and ``window``/``overlap`` need not
be stride multiples. The encoder input is right-padded with unavailable
bases to a stride multiple, as in training.

Windows are encoded and decoded in groups of ``decode_batch``: each group is
flushed through the scan and core claiming before the next group is
featurized, so the number of live emission tensors never exceeds
``decode_batch`` and host memory is bounded by the window size, not the
chromosome length. Only the GFF3 rows (proportional to the genes found, not
the bases) are retained.

Chains are mapped to genomic coordinates with :func:`model.grammar.gff3_rows`
over the whole oriented chromosome (shift by the window's oriented offset,
then ``n`` = chromosome length), so the plus/minus map is the reviewed one.

Timing is the caller's business (:func:`model.a.train.measure` wraps each
stage with a :class:`StageClock`): this module only makes the stages
explicit -- ``io`` (FASTA read and the minus-strand reverse complement),
``preprocess`` (featurizer), ``encoder``, ``decode`` (bias, scan, traceback)
and ``output`` (genomic map and GFF3 serialisation).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, NamedTuple, Optional, Sequence, Tuple

from model.grammar import Chain, Segment, gff3_rows, reverse_complement

# Overlap between consecutive windows (oriented bases): a complete chain up
# to half this long is inside the window that owns it (see :func:`tiles`),
# and its start has at least half this much context on each side (the
# chain's *end* may sit at the window edge). With the 12,288-base training
# windows the step is 8,192, so each strand is decoded 1.5 times over.
DEFAULT_OVERLAP = 4096


class Tile(NamedTuple):
    """One window of an oriented sequence: ``[start, end)`` is what is decoded,
    ``[core_start, core_end)`` is the part whose chains this window reports."""
    start: int
    end: int
    core_start: int
    core_end: int


def tiles(n: int, window: int, overlap: int) -> List[Tile]:
    """Cut ``[0, n)`` into windows of ``window`` bases stepping by
    ``window - overlap``; the last window is shortened to end at ``n``. Cores
    partition ``[0, n)`` exactly: the first core starts at 0, the last ends at
    ``n``, interior core boundaries sit ``overlap // 2`` bases into the next
    window. A complete chain of at most ``overlap - overlap // 2`` bases lies
    inside the window whose core holds its start. The context guarantee is
    about the chain's *start*: at least ``overlap // 2`` bases before it (in
    the interior; the first core begins at 0) and ``overlap - overlap // 2``
    after it, which a chain of the maximal length uses up entirely, so it
    can end on the window's last base with no right context. True sequence
    edges have only the sequence that exists."""
    if window < 1:
        raise ValueError("window must be positive")
    if not 0 <= overlap < window:
        raise ValueError("overlap must be in [0, window)")
    if n <= 0:
        return []
    step = window - overlap
    out: List[Tile] = []
    start = 0
    while True:
        end = min(start + window, n)
        out.append(Tile(start, end, 0, 0))  # cores filled below
        if end >= n:
            break
        start += step
    # Interior core boundaries sit ``half`` bases into the next window, so an
    # owned chain's start has at least ``half`` bases of sequence before it
    # and ``overlap - half`` after it; a chain no longer than ``overlap -
    # half`` therefore ends inside the window (possibly on its last base).
    half = overlap // 2
    result = []
    for i, t in enumerate(out):
        cs = 0 if i == 0 else t.start + half
        ce = n if i == len(out) - 1 else out[i + 1].start + half
        result.append(Tile(t.start, t.end, cs, ce))
    return result


@dataclass
class StageClock:
    """Accumulates process CPU and wall seconds per named stage. ``sync`` is
    called before reading the clocks (CUDA synchronize on a device)."""
    cpu: Dict[str, float] = field(default_factory=dict)
    wall: Dict[str, float] = field(default_factory=dict)
    sync: Optional[callable] = None

    def __call__(self, stage: str) -> "_Timer":
        return _Timer(self, stage)

    def total_cpu(self) -> float:
        return sum(self.cpu.values())

    def total_wall(self) -> float:
        return sum(self.wall.values())


class _Timer:
    def __init__(self, clock: StageClock, stage: str):
        self.clock, self.stage = clock, stage

    def __enter__(self):
        self.c0, self.w0 = time.process_time(), time.perf_counter()
        return self

    def __exit__(self, *exc):
        if self.clock.sync is not None:
            self.clock.sync()
        self.clock.cpu[self.stage] = self.clock.cpu.get(self.stage, 0.0) + time.process_time() - self.c0
        self.clock.wall[self.stage] = self.clock.wall.get(self.stage, 0.0) + time.perf_counter() - self.w0
        return False


def shifted(chain: Chain, offset: int) -> Chain:
    """``chain`` with every segment moved by ``offset`` oriented bases."""
    return Chain(segments=[Segment(s.kind, s.start + offset, s.end + offset, s.prefix_len, s.component)
                           for s in chain.segments],
                 partial_5=chain.partial_5, partial_3=chain.partial_3, uncertain=chain.uncertain)


def claimed(chains: Sequence[Chain], tile: Tile) -> List[Chain]:
    """The chains of ``tile`` (oriented to the window) whose 5' start lies in
    the tile's core, shifted to oriented-chromosome coordinates."""
    out = []
    for c in chains:
        if not c.segments:
            continue
        start = c.segments[0].start + tile.start
        if tile.core_start <= start < tile.core_end:
            out.append(shifted(c, tile.start))
    return out


def read_sequence(fasta: str, seqid: str) -> str:
    """The sequence of ``seqid`` from a (gzipped) FASTA, case preserved."""
    from model.labels.admission import _iter_fasta

    for name, seq in _iter_fasta(fasta):
        if name == seqid:
            return seq
    raise ValueError(f"{seqid!r} not in {fasta}")


def predict_sequence(model, seqid: str, seq: str, *, code, structure, tables,
                     window: int, overlap: int = DEFAULT_OVERLAP,
                     decode_batch: int = 16, device=None, dtype=None,
                     clock: Optional[StageClock] = None,
                     strands: Sequence[str] = ("+", "-"),
                     stride: Optional[int] = None,
                     segments: Optional[int] = None) -> Tuple[List[str], dict]:
    """Decode ``seq`` on ``strands`` and return ``(gff3 rows, counts)``.

    ``counts``: ``oriented_bases`` (sum of decoded window lengths, overlap
    included), ``strand_bases`` (``len(seq)`` per strand), ``windows``,
    ``chains`` (reported), ``chains_decoded`` (before core de-duplication).

    ``stride`` is the encoder's pooling stride (default ``POOL_STRIDE``); each
    window's encoder input starts at the stride multiple at or before the tile
    so the pooling grid is the oriented chromosome's, and windows are flushed
    through the decoder in groups of ``decode_batch`` (at most that many
    emission tensors are alive at once).

    With ``segments`` (proposal 3.3's carried seams) each strand is instead
    cut into about that many *segments* that overlap by ``overlap`` (the
    containment guarantee then holds at segment seams only); every segment is
    encoded in ``window``-base tiles that do not overlap, and all segments of
    all strands are decoded in one batch by :func:`model.a.fast_viterbi.
    scan_segments`, which scans ``window`` bases at a time and carries the
    Viterbi state across the interior tile seams, so a gene crossing such a
    seam is decoded whole. ``counts["windows"]`` is then the number of
    segments and ``counts["tiles"]`` the number of encoder tiles; the
    oversampling is ``1 + (segments - 1) * overlap / n`` per strand instead of
    ``window / (window - overlap)``. Each tile is encoded when the scan
    reaches it, so the live state is one tile of emissions and operands for
    every segment, the packed back-pointers of every segment (23 bytes per
    oriented base at ``K = 24``, ``R = 3``;
    :class:`model.a.fast_viterbi.PackedBackPointers`), and one segment's
    dense pointers (120 bytes per base) during its traceback; the part that
    grows with the chromosome is the packed store.
    """
    import torch

    from . import pooled
    from .encoder import POOL_STRIDE
    from .fast_viterbi import viterbi_windows
    from .features import encode_sequence
    from .train import pad_window

    if decode_batch < 1:
        raise ValueError("decode_batch must be positive")
    stride = POOL_STRIDE if stride is None else stride
    clock = clock or StageClock()
    n = len(seq)
    rows: List[str] = []
    counts = dict(oriented_bases=0, strand_bases=n * len(strands), windows=0,
                  chains=0, chains_decoded=0)
    gene = 0

    def flush(strand, oriented, group):
        nonlocal gene
        if not group:
            return
        tiles_, xs, ems = zip(*group)
        with clock("decode"):
            decoded = viterbi_windows(list(xs), list(ems), codes=[code] * len(xs),
                                      duration=structure, tables=tables,
                                      batch_size=decode_batch)
        with clock("output"):
            for tile, (_score, chains) in zip(tiles_, decoded):
                counts["chains_decoded"] += len(chains)
                for chain in claimed(chains, tile):
                    gene += 1
                    parent = f"{seqid}.g{gene}"
                    rows.extend(gff3_rows(chain, seqid, n, strand, parent))
                    counts["chains"] += 1

    if segments is not None:
        return _predict_segments(model, seqid, seq, code=code, structure=structure,
                                 tables=tables, window=window, overlap=overlap,
                                 segments=segments, device=device, dtype=dtype,
                                 clock=clock, strands=strands, stride=stride,
                                 rows=rows, counts=counts)

    with torch.no_grad():
        for strand in strands:
            with clock("io"):
                oriented = seq if strand == "+" else reverse_complement(seq)
            group: List[Tuple[Tile, str, "torch.Tensor"]] = []
            for tile in tiles(n, window, overlap):
                # Anchor the pooling grid to the oriented origin: encode from
                # the stride multiple at or before the tile, crop afterwards.
                enc_start = tile.start - tile.start % stride
                x = oriented[tile.start:tile.end]
                with clock("preprocess"):
                    padded, available = pad_window(oriented[enc_start:tile.end], stride)
                    feats = encode_sequence(padded, available=available).unsqueeze(0).to(device)
                with clock("encoder"):
                    emissions = model.encoder(feats)[0][:, tile.start - enc_start:tile.end - enc_start].to(dtype)
                with clock("decode"):
                    emissions = emissions + pooled.motif_bias(x, model.decoder, dtype=dtype, device=device)
                group.append((tile, x, emissions))
                counts["oriented_bases"] += len(x)
                counts["windows"] += 1
                if len(group) == decode_batch:
                    flush(strand, oriented, group)
                    group = []
            flush(strand, oriented, group)
            del group
    return rows, counts


def segment_length(n: int, segments: int, overlap: int) -> int:
    """The window length that cuts ``[0, n)`` into about ``segments`` tiles
    overlapping by ``overlap`` (:func:`tiles`); at least ``overlap + 1``."""
    if segments < 1:
        raise ValueError("segments must be positive")
    return max(overlap + 1, -(-(n + (segments - 1) * overlap) // segments))


def _predict_segments(model, seqid, seq, *, code, structure, tables, window, overlap,
                      segments, device, dtype, clock, strands, stride, rows, counts):
    """The ``segments`` mode of :func:`predict_sequence`: one carried-state
    scan over every segment of every strand, fed one ``window``-base tile at
    a time; each tile is encoded when the scan asks for it, so the emissions
    alive at once are one tile of every segment, never a chromosome
    (engels-0090, stalin-0093)."""
    import torch

    from . import pooled
    from .encoder import POOL_STRIDE
    from .fast_viterbi import viterbi_segments
    from .features import encode_sequence
    from .train import pad_window

    stride = POOL_STRIDE if stride is None else stride
    n = len(seq)
    counts["tiles"] = 0
    seg_len = segment_length(n, segments, overlap)
    owners: List[Tuple[str, Tile]] = []
    xs: List[str] = []
    ems: List[callable] = []

    def emitter(oriented: str, seg: Tile, x: str):
        def emit(start: int, end: int) -> "torch.Tensor":
            # Encoder input anchored to the pooling grid as in the window
            # mode; the dinucleotide bias reads one base past a tile's end
            # and two before its start, within the segment.
            a, b = seg.start + start, seg.start + end
            enc_start = a - a % stride
            with clock("preprocess"):
                padded, available = pad_window(oriented[enc_start:b], stride)
                feats = encode_sequence(padded, available=available).unsqueeze(0).to(device)
            with clock("encoder"):
                em = model.encoder(feats)[0][:, a - enc_start:b - enc_start].to(dtype)
            lo, hi = max(0, start - 2), min(len(x), end + 1)
            bias = pooled.motif_bias(x[lo:hi], model.decoder, dtype=dtype, device=device)
            counts["tiles"] += 1
            return em + bias[:, start - lo:end - lo]
        return emit

    with torch.no_grad():
        for strand in strands:
            with clock("io"):
                oriented = seq if strand == "+" else reverse_complement(seq)
            for seg in tiles(n, seg_len, overlap):
                x = oriented[seg.start:seg.end]
                owners.append((strand, seg))
                xs.append(x)
                ems.append(emitter(oriented, seg, x))
                counts["oriented_bases"] += len(x)
                counts["windows"] += 1
        # The encoder runs inside the scan; its stages are clocked by the
        # emitter and taken back out of ``decode``.
        inner = {s: (clock.cpu.get(s, 0.0), clock.wall.get(s, 0.0)) for s in ("preprocess", "encoder")}
        with clock("decode"):
            decoded = viterbi_segments(xs, ems, tile=window, code=code, duration=structure,
                                       tables=tables, dtype=dtype, device=device)
        for s, (c0, w0) in inner.items():
            clock.cpu["decode"] -= clock.cpu.get(s, 0.0) - c0
            clock.wall["decode"] -= clock.wall.get(s, 0.0) - w0
        del ems
        gene = 0
        with clock("output"):
            for (strand, seg), (_score, chains) in zip(owners, decoded):
                counts["chains_decoded"] += len(chains)
                for chain in claimed(chains, seg):
                    gene += 1
                    rows.extend(gff3_rows(chain, seqid, n, strand, f"{seqid}.g{gene}"))
                    counts["chains"] += 1
    return rows, counts


def write_gff3(path: str, rows: Sequence[str]) -> int:
    """Write ``rows`` with a GFF3 header; returns bytes written."""
    text = "##gff-version 3\n" + "".join(r + "\n" for r in rows)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return len(text.encode("utf-8"))
