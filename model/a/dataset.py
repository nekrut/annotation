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

    @property
    def n(self) -> int:
        return len(self.window)

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


def _oriented_window(t, seq: str) -> str:
    """The focal window's oriented sequence with case preserved.

    Same ``[ws, we]`` clamp as ``oriented_chain`` (which stays the authority for
    the CDS/intron coordinates); this only avoids that helper's ``.upper()`` so
    the soft-mask channel survives. ``revcomp`` preserves case on the minus
    strand.
    """
    s, e = t.span
    ws, we = max(1, s - FLANK), min(len(seq), e + FLANK)
    window = seq[ws - 1:we]
    if t.strand == "-":
        window = revcomp(window)
    return window


def iter_windows(
    summary: str,
    gff: str,
    fasta: str,
    *,
    complete_only: bool = True,
    max_window: Optional[int] = None,
    stats: Optional[LoaderStats] = None,
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
    counts) any window longer than it; ``None`` keeps every window.
    """
    if stats is None:
        stats = LoaderStats()

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
            window = _oriented_window(t, seq)
            if window.upper() != _up_window:
                raise RuntimeError(
                    f"cased window disagrees with oriented_chain for {t.key}")
            if max_window is not None and len(window) > max_window:
                stats.skipped_too_long += 1
                continue
            stats.yielded += 1
            stats.windows_by_seqid[t.seqid] += 1
            yield WindowExample(
                key=t.key,
                window=window,
                cds_ranges=cds,
                intron_ranges=introns,
                table=table,
                partial_5=bool(mt.five),
                partial_3=bool(mt.three),
            )
