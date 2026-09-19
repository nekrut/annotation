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
- **Complete targets only, for now.** `model.a.loss.chain_nll` is scoped to
  *complete* chains (no edge prior); a chain whose own 5' or 3' end is declared
  partial needs the section-3.6 boundary-support increment. `iter_windows`
  therefore skips edge-partial admitted chains by default and counts them, so
  the training set is exactly the complete-target oracle's domain. Set
  ``complete_only=False`` to receive them too (with their flags) once the
  boundary loss lands.
- **No silent length cap.** A complete gene's oriented window is the CDS span
  plus flank; ~40% of admitted train chains are longer than the 3,072-base
  encoder core (engels-0047). Cropping such genes into the core with retained
  intron/codon state is the separate section-3.6 crop increment. Until then a
  caller may pass ``max_window`` to skip (and count) windows that would not fit
  a chosen chunk; the default keeps every window and never turns the encoder
  core into a gene-length cap.

Checksum verification (`verify_source`) is a hard precondition: a window built
from a FASTA/GFF3 that does not match the manifest summary's recorded digests
is not the audited label, so this refuses to yield from mismatched inputs.

Standard library plus the admission module only; `WindowExample.features`
imports torch lazily, exactly as `model.a.features` does, so the loader and its
tests run on a host without torch.
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
    its length. ``cds_ranges`` and ``intron_ranges`` are the merged, oriented,
    half-open CDS and intron intervals in that window, the exact convention
    ``model.a.loss.numerator_scores`` consumes. ``table`` is the genetic-code id
    (input metadata, not inferred). ``key`` is the source ``(seqid, strand,
    transcript-id)`` identity for provenance.
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

        ``oriented_chain`` cuts the window from real FASTA and clamps to the
        sequence edges, so every position is real DNA (no invented padding);
        availability is therefore the default all-real mask."""
        from model.a.features import encode_sequence

        return encode_sequence(self.window)


@dataclass
class LoaderStats:
    admitted: int = 0
    yielded: int = 0
    skipped_partial: int = 0
    skipped_too_long: int = 0
    windows_by_seqid: Dict[str, int] = field(default_factory=lambda: defaultdict(int))


def iter_windows(
    gff: str,
    fasta: str,
    *,
    m: int = 20,
    table: int = 1,
    complete_only: bool = True,
    max_window: Optional[int] = None,
    stats: Optional[LoaderStats] = None,
) -> Iterator[WindowExample]:
    """Yield one :class:`WindowExample` per admitted representative chain.

    Admission is delegated to ``audit_species`` (same ``m`` and ``table`` as the
    manifest), so the yielded set is exactly the admitted representatives. Each
    is rejoined to the checksummed FASTA and oriented with
    ``model.labels.numerator_check.oriented_chain`` (the flank-padded, merged,
    strand-corrected windowing the numerator check already verifies).

    When ``complete_only`` (the default), edge-partial admitted chains are
    skipped and counted in ``stats.skipped_partial``; they need the boundary
    support the current loss does not yet score. ``max_window`` skips (and
    counts) any window longer than it, for a caller that trains on a fixed
    chunk; ``None`` keeps every window.
    """
    if stats is None:
        stats = LoaderStats()

    res = audit_species("_loader_", gff, fasta, m=m, table=table, log=lambda *a, **k: None)
    admitted = {
        (r["seqid"], r["strand"], r["transcript"])
        for r in res.manifest_rows
        if r["status"] == "admitted"
    }
    stats.admitted = len(admitted)

    transcripts, _gene_biotype, _shim = load_gff_rows(gff)
    by_seqid: Dict[str, List] = defaultdict(list)
    for t in transcripts:
        if t.key in admitted:
            by_seqid[t.seqid].append(t)

    for name, seq in _iter_fasta(fasta):
        for t in by_seqid.get(name, ()):  # only admitted transcripts on this seqid
            mt = metadata_flags(t, m, table)
            if complete_only and (mt.five or mt.three):
                stats.skipped_partial += 1
                continue
            window, cds, introns = oriented_chain(t, seq)
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
