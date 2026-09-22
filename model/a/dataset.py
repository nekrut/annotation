"""Structured training-window loader for candidate A (proposal section 3.6).

The chain loss (`model.a.loss`) scores one admitted gold chain given the
encoder's eleven emission channels and the chain's oriented CDS/intron
intervals. This module produces those `(window, cds_ranges, intron_ranges)`
triples from a species' *checksummed* source GFF3 and FASTA, so a training
step can pair each window's featurized DNA (`model.a.features`) with its
support mask (`model.a.loss.numerator_scores`).

Design decisions, all following the admission audit (T-human-013) rather than
re-deriving it:

- **One admission authority.** Which transcripts are admitted, and with which
  boundary flags, is decided by `model.labels.admission.audit_species` exactly
  as the committed manifests were built. This loader calls it and keeps only
  the transcripts it marks `admitted`; it never re-implements the metadata or
  sequence audit. The manifests omit exon/site coordinates, so the windows are
  rejoined from the raw GFF3/FASTA (which the checksum gate below pins), not
  read back from the manifest.
- **The pinned summary is the loading interface.** `iter_windows` takes the
  species `*.summary.json`, hard-verifies the GFF3/FASTA digests against it
  (`verify_source`), and reads the audit's `m` and `table` from it, so the
  loader can never audit arbitrary inputs or run the audit with settings that
  disagree with the committed manifest. There is no way to yield a window
  without passing this gate.
- **Complete targets on clean windows only, for now.** `model.a.loss.chain_nll`
  is scoped to *complete* chains (no edge prior), and `numerator_scores` forces
  every non-CDS/intron base in the window to intergenic ``U``. That is only
  correct when the flanks are genuinely intergenic: a neighbouring gene's CDS
  falling in a flank must not be supervised as background, and a masked span
  must stay unconstrained (proposal section 3.6). This increment therefore
  restricts to *clean* windows -- no transcript of any other gene overlaps the
  focal window -- and skips (and counts) the rest. Adjacent-gene support built
  from the full window's annotations, and edge-partial boundary support, are
  later section-3.6 increments; `complete_only=False` receives edge-partial
  chains (with their flags) once the boundary loss lands.
- **Intergenic context, clipped at the neighbours (fit v2).** With
  ``context > 0`` a clean window is widened by up to ``context`` bases on each
  side, but never past the nearest span of a *different* gene (any transcript,
  admitted or not) or the sequence end, so every added base is annotated
  intergenic and the ``U`` supervision ``numerator_scores`` applies to it stays
  correct without a loss change. Fit v1 supervised only 3 % of its sampled
  bases as intergenic (10-base flanks plus background tiles) and its decoder
  fused genes through intergenic sequence (a-pilot 3.3); this is the
  minimal loader increment that adds real gene-adjacent intergenic sequence
  to every chain window. The clean check itself is unchanged (a gene inside
  the 10-base flank still skips the window).
- **No silent length cap.** A complete gene's oriented window is the CDS span
  plus flank; ~40% of admitted train chains are longer than the 3,072-base
  encoder core (engels-0047). Cropping such genes into the core with retained
  intron/codon state is the separate section-3.6 crop increment. Until then a
  caller may pass ``max_window`` to skip (and count) windows that would not fit
  a chosen chunk; the default keeps every window and never turns the encoder
  core into a gene-length cap.

Soft masking is carried through orientation: the window preserves lower-case,
so the featurizer's soft-mask channel (channel 5) observes it. This is why the
window is rebuilt from the raw FASTA slice here rather than reusing
``oriented_chain``'s upper-cased grammar-check window; ``oriented_chain`` stays
the single authority for the CDS/intron coordinates.

Standard library plus the admission module only; `WindowExample.features`
imports torch lazily, exactly as `model.a.features` does, so the loader and its
tests run on a host without torch. The loader is a *checkout-time* tool: the
admission audit loads `benchmark/score.py` by path, so it runs from a repository
checkout (where the source data and that filter live), not an installed wheel.
`model.a.__init__` therefore imports this module lazily.
"""
from __future__ import annotations

import json
import random
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from model.labels.admission import (
    audit_species,
    load_gff_rows,
    md5_file,
    metadata_flags,
    revcomp,
    _iter_fasta,
    _open,
)
from model.labels.numerator_check import FLANK, oriented_chain
from model.a.loss import numerator_scores

Range = Tuple[int, int]


class SourceMismatch(RuntimeError):
    """A source file's digest does not match the manifest summary."""


def verify_source(summary_path: str, gff: str, fasta: Optional[str] = None) -> Dict[str, str]:
    """Check ``gff`` (and ``fasta`` when given) against the digests a species'
    ``*.summary.json`` recorded when its manifest was built.

    Returns the verified ``{"gff_md5": ..., "fasta_md5": ...}`` on success and
    raises :class:`SourceMismatch` on the first disagreement. A summary that
    recorded a FASTA digest but is verified without a FASTA (or vice versa) is a
    mismatch: the audited labels depend on both files.
    """
    with open(summary_path, "r", encoding="utf-8") as fh:
        summary = json.load(fh)
    want_gff = summary.get("gff_md5")
    if want_gff is None:
        raise SourceMismatch(f"{summary_path} records no gff_md5")
    got_gff = md5_file(gff)
    if got_gff != want_gff:
        raise SourceMismatch(f"gff md5 {got_gff} != manifest {want_gff}")
    want_fasta = summary.get("fasta_md5")
    if bool(want_fasta) != bool(fasta):
        raise SourceMismatch(
            "fasta presence disagrees with the manifest: "
            f"manifest fasta_md5={want_fasta!r}, fasta={fasta!r}"
        )
    got_fasta = md5_file(fasta) if fasta else None
    if want_fasta and got_fasta != want_fasta:
        raise SourceMismatch(f"fasta md5 {got_fasta} != manifest {want_fasta}")
    return {"gff_md5": got_gff, "fasta_md5": got_fasta}


@dataclass
class WindowExample:
    """One admitted complete chain as an oriented training window.

    ``window`` is the oriented sequence (both strands processed as ``+``), with
    ``FLANK`` unrewarded bases on each side where the sequence allows; ``n`` is
    its length. Case is preserved through orientation, so a soft-masked
    (lower-case) base stays lower-case for the featurizer's soft-mask channel.
    ``cds_ranges`` and ``intron_ranges`` are the merged, oriented, half-open CDS
    and intron intervals in that window, the exact convention
    ``model.a.loss.numerator_scores`` consumes. ``table`` is the genetic-code id
    (from the pinned summary, not inferred). ``key`` is the source ``(seqid,
    strand, transcript-id)`` identity for provenance.
    """

    key: Tuple[str, str, str]
    window: str
    cds_ranges: List[Range]
    intron_ranges: List[Range]
    table: int
    partial_5: bool = False
    partial_3: bool = False
    flank: int = FLANK
    # Intergenic context actually added beyond ``flank`` on the oriented 5'
    # and 3' sides (0 when the loader ran without ``context`` or the window
    # was clipped at a neighbouring gene or the sequence end).
    context_5: int = 0
    context_3: int = 0

    @property
    def n(self) -> int:
        return len(self.window)

    @property
    def cds_bases(self) -> int:
        return sum(b - a for a, b in self.cds_ranges)

    @property
    def intron_bases(self) -> int:
        return sum(b - a for a, b in self.intron_ranges)

    @property
    def intergenic_bases(self) -> int:
        """Bases the support mask supervises as intergenic ``U``: everything
        outside the CDS/intron cover (flanks, context, or the whole window of
        a gene-free background tile)."""
        return self.n - self.cds_bases - self.intron_bases

    def support(self, base=None):
        """The hard support mask for this chain's numerator, over ``base``
        (the encoder emissions; zeros when ``None``). Delegates to
        :func:`model.a.loss.numerator_scores`."""
        return numerator_scores(base, self.n, self.cds_ranges, self.intron_ranges)

    def features(self):
        """The [8, n] feature tensor for the window (requires torch).

        The window is a real FASTA slice clamped to the sequence edges, so every
        position is real DNA (no invented padding); availability is therefore
        the default all-real mask. Case is preserved, so channel 5 (soft
        masking) is populated from the window's lower-case bases."""
        from model.a.features import encode_sequence

        return encode_sequence(self.window)


@dataclass
class LoaderStats:
    admitted: int = 0
    yielded: int = 0
    skipped_partial: int = 0
    skipped_too_long: int = 0
    skipped_neighbor: int = 0
    windows_by_seqid: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    # Gene-free (background) windows, counted apart from the chain windows.
    background_candidates: int = 0
    background_yielded: int = 0
    # Intergenic context bases added to yielded chain windows (``context``),
    # and how many of those windows were clipped short of the requested
    # context by a neighbouring gene or the sequence end on at least one side.
    context_bases: int = 0
    context_clipped: int = 0


@dataclass
class CoverageRow:
    """One species' admission/training coverage under the current loader scope.

    ``admitted`` is what the audit accepted; ``yielded`` is what this loader
    actually trains on given the complete-target, clean-window scope; the three
    ``skipped_*`` counts partition the difference. Reviewers asked for this
    accounting before any train-panel conclusion is drawn from A, because the
    clean-window and complete-only restrictions are temporary section-3.6
    scope, not the final training set (engels-0065, stalin-0068).
    """

    species: str
    admitted: int
    yielded: int
    skipped_partial: int
    skipped_too_long: int
    skipped_neighbor: int

    @property
    def yielded_fraction(self) -> float:
        return self.yielded / self.admitted if self.admitted else 0.0


def coverage_report(
    sources: Sequence[Tuple[str, str, str, str]],
    *,
    complete_only: bool = True,
    max_window: Optional[int] = None,
) -> List[CoverageRow]:
    """Run the loader over each species and return its training coverage.

    ``sources`` is a sequence of ``(species, summary, gff, fasta)``. Each entry
    is consumed exactly as :func:`iter_windows` would in training (same checksum
    gate, admission delegation, and ``complete_only``/``max_window`` scope), and
    only the counts are kept -- the windows themselves are discarded, so this is
    a cheap pass that quantifies how much of each species' admitted set the
    current scope actually reaches. It never re-derives the audit.
    """
    rows: List[CoverageRow] = []
    for species, summary, gff, fasta in sources:
        stats = LoaderStats()
        for _ in iter_windows(
            summary, gff, fasta,
            complete_only=complete_only, max_window=max_window, stats=stats,
        ):
            pass
        rows.append(CoverageRow(
            species=species,
            admitted=stats.admitted,
            yielded=stats.yielded,
            skipped_partial=stats.skipped_partial,
            skipped_too_long=stats.skipped_too_long,
            skipped_neighbor=stats.skipped_neighbor,
        ))
    return rows


def format_coverage(rows: Sequence[CoverageRow]) -> str:
    """Render :func:`coverage_report` rows as a TSV table with a TOTAL row.

    Columns: species, admitted, yielded, yielded_pct, skip_partial,
    skip_neighbor, skip_too_long. The same convention as
    ``docs/cost-baseline/measured.tsv``: tab-separated, one header line, a
    trailing ``TOTAL`` aggregating every species.
    """
    header = ("species", "admitted", "yielded", "yielded_pct",
              "skip_partial", "skip_neighbor", "skip_too_long")
    lines = ["\t".join(header)]
    tot = CoverageRow("TOTAL", 0, 0, 0, 0, 0)
    for r in rows:
        lines.append("\t".join((
            r.species, str(r.admitted), str(r.yielded),
            "%.1f" % (100.0 * r.yielded_fraction),
            str(r.skipped_partial), str(r.skipped_neighbor),
            str(r.skipped_too_long),
        )))
        tot.admitted += r.admitted
        tot.yielded += r.yielded
        tot.skipped_partial += r.skipped_partial
        tot.skipped_neighbor += r.skipped_neighbor
        tot.skipped_too_long += r.skipped_too_long
    lines.append("\t".join((
        tot.species, str(tot.admitted), str(tot.yielded),
        "%.1f" % (100.0 * tot.yielded_fraction),
        str(tot.skipped_partial), str(tot.skipped_neighbor),
        str(tot.skipped_too_long),
    )))
    return "\n".join(lines) + "\n"


def _oriented_window(t, seq: str, ws: Optional[int] = None,
                     we: Optional[int] = None) -> str:
    """The focal window's oriented sequence with case preserved.

    By default the same ``[ws, we]`` clamp as ``oriented_chain`` (which stays
    the authority for the CDS/intron coordinates); this only avoids that
    helper's ``.upper()`` so the soft-mask channel survives. ``revcomp``
    preserves case on the minus strand. Explicit 1-based inclusive ``ws``/``we``
    select a wider (context) slice.
    """
    s, e = t.span
    if ws is None:
        ws = max(1, s - FLANK)
    if we is None:
        we = min(len(seq), e + FLANK)
    window = seq[ws - 1:we]
    if t.strand == "-":
        window = revcomp(window)
    return window


def _context_bounds(ws: int, we: int, seq_len: int, context: int, gid: str,
                    others) -> Tuple[int, int]:
    """Widen the clean window ``[ws, we]`` by up to ``context`` bases per side,
    stopping at the sequence ends and at the nearest span of any *other* gene
    (so the added bases are annotated intergenic). Returns the 1-based
    inclusive bounds; equal to the input when ``context`` is 0."""
    if context <= 0:
        return ws, we
    left_gene_end = 0
    right_gene_start = seq_len + 1
    for os_, oe, g in others:
        if g == gid:
            continue
        if oe < ws and oe > left_gene_end:
            left_gene_end = oe
        elif os_ > we and os_ < right_gene_start:
            right_gene_start = os_
    return (max(1, ws - context, left_gene_end + 1),
            min(seq_len, we + context, right_gene_start - 1))


def iter_windows(
    summary: str,
    gff: str,
    fasta: str,
    *,
    complete_only: bool = True,
    max_window: Optional[int] = None,
    stats: Optional[LoaderStats] = None,
    context: int = 0,
) -> Iterator[WindowExample]:
    """Yield one :class:`WindowExample` per admitted representative chain.

    ``summary`` is the species' ``*.summary.json``. Its GFF3/FASTA digests are
    hard-verified (:func:`verify_source`, raising :class:`SourceMismatch` on any
    disagreement) and its audit ``m``/``table`` are read from it, before any
    window is yielded, so the loader cannot audit unpinned inputs or use
    settings that disagree with the committed manifest.

    Admission is delegated to ``audit_species`` (with the summary's ``m`` and
    ``table``), so the yielded set is exactly the admitted representatives. Each
    is rejoined to the checksummed FASTA and oriented; ``oriented_chain`` gives
    the merged CDS/intron coordinates and :func:`_oriented_window` the
    case-preserving window.

    Only *clean* windows are yielded: a window is skipped (and counted in
    ``stats.skipped_neighbor``) when any transcript of a *different* gene
    overlaps it, because ``numerator_scores`` would otherwise force that
    neighbour's coding (or masked) bases to intergenic ``U``. When
    ``complete_only`` (the default), edge-partial admitted chains are also
    skipped and counted in ``stats.skipped_partial``. ``max_window`` skips (and
    counts) any window longer than it (after context); ``None`` keeps every
    window.

    ``context`` widens each clean window by up to that many bases per side of
    annotated intergenic sequence, clipped at the nearest other gene's span
    and the sequence ends (:func:`_context_bounds`); the CDS/intron
    coordinates are shifted accordingly and the example records the bases
    actually added as ``context_5``/``context_3`` (oriented).
    """
    if stats is None:
        stats = LoaderStats()
    if context < 0:
        raise ValueError("context must not be negative")

    verify_source(summary, gff, fasta)
    with open(summary, "r", encoding="utf-8") as fh:
        meta = json.load(fh)
    m = int(meta.get("m", 20))
    table = int(meta.get("table", 1))

    res = audit_species("_loader_", gff, fasta, m=m, table=table, log=lambda *a, **k: None)
    admitted = {
        (r["seqid"], r["strand"], r["transcript"])
        for r in res.manifest_rows
        if r["status"] == "admitted"
    }
    stats.admitted = len(admitted)

    transcripts, _gene_biotype, _shim = load_gff_rows(gff)
    # Every gene's span on each sequence, for the clean-window neighbour check.
    # Non-representative and masked transcripts count too; isoforms of the same
    # gene do not (same gene id), so a multi-isoform locus is not falsely dirty.
    genes_by_seqid: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for t in transcripts:
        s, e = t.span
        genes_by_seqid[t.seqid].append((s, e, t.gid))

    by_seqid: Dict[str, List] = defaultdict(list)
    for t in transcripts:
        if t.key in admitted:
            by_seqid[t.seqid].append(t)

    for name, seq in _iter_fasta(fasta):
        others = genes_by_seqid.get(name, ())
        for t in by_seqid.get(name, ()):  # only admitted transcripts on this seqid
            mt = metadata_flags(t, m, table)
            if complete_only and (mt.five or mt.three):
                stats.skipped_partial += 1
                continue
            s, e = t.span
            ws, we = max(1, s - FLANK), min(len(seq), e + FLANK)
            if any(gid != t.gid and not (oe < ws or os_ > we)
                   for os_, oe, gid in others):
                stats.skipped_neighbor += 1
                continue
            _up_window, cds, introns = oriented_chain(t, seq)
            ws2, we2 = _context_bounds(ws, we, len(seq), context, t.gid, others)
            window = _oriented_window(t, seq, ws2, we2)
            # Oriented offset of the flank-only window inside the wider one.
            if t.strand == "-":
                shift, ctx5, ctx3 = we2 - we, we2 - we, ws - ws2
            else:
                shift, ctx5, ctx3 = ws - ws2, ws - ws2, we2 - we
            if window.upper()[shift:shift + len(_up_window)] != _up_window:
                raise RuntimeError(
                    f"cased window disagrees with oriented_chain for {t.key}")
            if shift:
                cds = [(a + shift, b + shift) for a, b in cds]
                introns = [(a + shift, b + shift) for a, b in introns]
            if max_window is not None and len(window) > max_window:
                stats.skipped_too_long += 1
                continue
            stats.yielded += 1
            stats.windows_by_seqid[t.seqid] += 1
            stats.context_bases += ctx5 + ctx3
            if context and (ctx5 < context or ctx3 < context):
                stats.context_clipped += 1
            yield WindowExample(
                key=t.key,
                window=window,
                cds_ranges=cds,
                intron_ranges=introns,
                table=table,
                partial_5=bool(mt.five),
                partial_3=bool(mt.three),
                context_5=ctx5,
                context_3=ctx3,
            )


# GFF3 feature types whose rows block background. Everything the source
# annotates as a gene, a transcript, or a piece of one counts (protein-coding
# or not, admitted or not): ``gene``/``pseudogene`` and their Ensembl variants
# (``ncRNA_gene``, ``*_gene_segment``), every ``*RNA`` transcript type, and
# the exon/CDS/UTR rows themselves so an orphan transcript without a parent
# gene row is still blocked. Landmarks (``region``, ``chromosome``) and
# non-gene features (``centromere``, ``origin_of_replication``,
# ``long_terminal_repeat``, ``mobile_genetic_element``, ...) are not genes
# and stay eligible.
_GENE_TYPES = frozenset({
    "gene", "pseudogene", "exon", "CDS", "transcript", "primary_transcript",
    "pseudogenic_transcript", "pseudogenic_exon",
    "five_prime_UTR", "three_prime_UTR", "UTR", "start_codon", "stop_codon",
})


def is_gene_feature(ftype: str) -> bool:
    """True when a GFF3 ``type`` column value marks a gene, transcript, or a
    part of one (see ``_GENE_TYPES``)."""
    return (ftype in _GENE_TYPES or ftype.endswith("gene")
            or ftype.endswith("gene_segment") or ftype.endswith("RNA"))


def annotated_gene_spans(gff: str) -> Dict[str, List[Tuple[int, int]]]:
    """1-based inclusive ``[start, end]`` spans per seqid of every raw GFF3
    row whose type :func:`is_gene_feature`, regardless of strand, parent,
    or whether the admission parser kept it.

    This is deliberately broader than ``load_gff_rows``: that parser keeps
    only transcripts with CDS rows and ``Transcript.span`` covers the CDS
    endpoints, so a UTR, an ncRNA, or a pseudogene would be invisible to it
    and could be tiled as background (engels-0096).
    """
    spans: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with _open(gff) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or not is_gene_feature(f[2]):
                continue
            try:
                s, e = int(f[3]), int(f[4])
            except ValueError:
                continue
            spans[f[0]].append((min(s, e), max(s, e)))
    return spans


def gene_free_intervals(seq_len: int, spans: Sequence[Tuple[int, int]],
                        margin: int = FLANK) -> List[Range]:
    """Zero-based half-open intervals of a sequence that no gene touches.

    ``spans`` are 1-based inclusive gene spans (:func:`annotated_gene_spans`);
    each is widened by ``margin`` on both sides so a background window never
    abuts a gene closer than a chain window's own flank. Every annotated gene
    and transcript counts, coding or not, admitted or not, masked or not: a
    base is background only if the source annotation has no gene feature
    there on either strand.
    """
    if seq_len <= 0:
        return []
    blocked = sorted((max(0, s - 1 - margin), min(seq_len, e + margin)) for s, e in spans)
    out: List[Range] = []
    cursor = 0
    for a, b in blocked:
        if a > cursor:
            out.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < seq_len:
        out.append((cursor, seq_len))
    return out


def iter_background_windows(
    summary: str,
    gff: str,
    fasta: str,
    *,
    length: int,
    count: int,
    seed: int = 0,
    stats: Optional[LoaderStats] = None,
) -> Iterator[WindowExample]:
    """Yield up to ``count`` gene-free windows of ``length`` bases (proposal
    section 3, "gene-free/background windows").

    The source pair is checksum-gated exactly as :func:`iter_windows`, and
    only sequences carrying at least one *admitted* chain are sampled (the
    audit's own inventory, so a mitochondrion or other sequence the admission
    excluded -- a different genetic code, no admitted representative -- is
    never background for the nuclear model). Every gene-free interval of such
    a sequence (:func:`gene_free_intervals` over :func:`annotated_gene_spans`:
    every raw gene, pseudogene, transcript, exon, CDS or UTR row of the GFF3
    on both strands, including UTRs and non-coding genes the admission parser
    never sees, widened by ``FLANK``) is tiled with non-overlapping
    ``length``-base candidates from its left edge; a deterministic
    ``random.Random(seed)`` draw without replacement picks ``count`` of them
    (all, when fewer exist) and gives each a strand, so the minus-strand half
    is reverse-complemented with case preserved. The yielded example has empty
    ``cds_ranges``/``intron_ranges``: its support mask is the all-intergenic
    path (``numerator_scores`` with no chain). ``key`` is
    ``(seqid, strand, "background:<start>-<end>")`` in zero-based half-open
    forward coordinates, so the dev split by seqid applies unchanged.
    """
    if length <= 0:
        raise ValueError("background window length must be positive")
    if count < 0:
        raise ValueError("background window count must not be negative")
    if stats is None:
        stats = LoaderStats()

    verify_source(summary, gff, fasta)
    with open(summary, "r", encoding="utf-8") as fh:
        meta = json.load(fh)
    m = int(meta.get("m", 20))
    table = int(meta.get("table", 1))

    res = audit_species("_loader_", gff, fasta, m=m, table=table, log=lambda *a, **k: None)
    admitted_seqids = {r["seqid"] for r in res.manifest_rows if r["status"] == "admitted"}

    spans_by_seqid = annotated_gene_spans(gff)

    candidates: List[Tuple[str, int]] = []
    sequences: Dict[str, str] = {}
    for name, seq in _iter_fasta(fasta):
        if name not in admitted_seqids:
            continue
        free = gene_free_intervals(len(seq), spans_by_seqid.get(name, ()))
        starts = [a for lo, hi in free for a in range(lo, hi - length + 1, length)]
        if starts:
            sequences[name] = seq
            candidates.extend((name, a) for a in starts)
    stats.background_candidates = len(candidates)

    rng = random.Random(seed)
    chosen = candidates if count >= len(candidates) else rng.sample(candidates, count)
    for name, a in chosen:
        strand = "+" if rng.random() < 0.5 else "-"
        window = sequences[name][a:a + length]
        if strand == "-":
            window = revcomp(window)
        stats.background_yielded += 1
        yield WindowExample(
            key=(name, strand, f"background:{a}-{a + length}"),
            window=window,
            cds_ranges=[],
            intron_ranges=[],
            table=table,
        )
