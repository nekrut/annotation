#!/usr/bin/env python3
"""Score a predicted GFF3 against a reference annotation.

Implements section 4 of ``docs/benchmark.md`` for one species and emits one
JSON object carrying every field section 4.8 asks for that can be computed
from annotation alone: nucleotide, CDS-exon, splice-site, transcript, locus,
and start/stop-codon metrics with their stratifications.  BUSCO/OMArk (4.6)
and cost (4.7) are measured by other tools and merged in by ``report.py``.

Standard library only.  Nothing is written except the JSON on stdout (or
``--out``).  Reference and prediction may be plain or gzipped GFF3.

Everything is scored on the **CDS**, not the transcript: UTRs are out of
scope for this charter and a model that does not predict them must not be
penalized (section 4).  A reference CDS chain is assumed to include its stop
codon, which is the RefSeq and Ensembl GFF3 convention; pass
``--stop-outside-cds`` if the prediction follows the other one.

Refuses to run without the submission declaration of section 3.3
(``--declaration``), because a run without it is not comparable.  The
exception is ``--self-test``, which scores built-in fixtures and checks the
answers.

Usage:
    python3 score.py --reference REF.gff.gz --prediction PRED.gff3 \\
        --species Homo_sapiens --declaration decl.yaml [--genome REF.fna.gz]
    python3 score.py --self-test
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import io
import json
import math
import os
import sys
import tempfile
from collections import defaultdict, deque

# Section 3.3.  A declaration missing any of these is not a submission.
DECLARATION_KEYS = [
    "model",
    "training_species",
    "pretraining_corpus",
    "heldout_seen_in_pretraining",
    "protein_db",
    "alignment",
    "alignment_rows_dropped",
    "informants",
    "rnaseq",
    "hardware",
]

# Section 4: "no unplaced scaffolds under 10 kb".
MIN_SEQ_LEN = 10000
# Section 4.3: local GC is measured in a 200 bp window centred on the site.
GC_WINDOW = 200
# Not every gap between consecutive CDS blocks is an intron.  RefSeq encodes a
# programmed ribosomal frameshift as two CDS blocks separated by 1 bp, which is
# 47 of the 343 CDS gaps in the S. cerevisiae reference alone.  The shortest
# known spliceosomal introns are around 30 nt, so gaps under this threshold are
# counted and reported separately instead of being scored as splice sites.
MIN_INTRON = 20
GC_BINS = [0.30, 0.40, 0.50, 0.60]  # five bins: <30, 30-40, 40-50, 50-60, >=60

# Section 4.5, stop codons by NCBI translation table, used only by the
# genome-based convention check below.  Only the tables the panel uses are
# listed: table 1 for nineteen species, and table 6 for Tetrahymena
# thermophila, where TAA and TAG code for glutamine and TGA is the only stop.
# An unlisted table is an error rather than a silent fall back to table 1,
# because falling back is exactly how a code-6 genome gets judged as if two
# thirds of its stops were real.
STOP_CODONS = {1: ("TAA", "TAG", "TGA"), 6: ("TGA",)}
# How many predicted chains the genome check needs before it will return a
# verdict.  Three is enough for the fixtures; a real prediction brings
# thousands, and the count is reported so a thin verdict is visible.
MIN_CONVENTION_CHAINS = 3

TX_TYPES = {"mRNA", "transcript", "CDS_transcript"}

# Section 4 scores the nuclear primary assembly.  Organelle genomes are out of
# charter scope and use their own genetic codes (human chrM is code 2, which
# would make every stop-codon check on it wrong), so they are dropped.  RefSeq
# GFF3 states this in the ``genome=`` attribute of the ``region`` feature; the
# values below are the ones NCBI uses.
ORGANELLE_GENOMES = {
    "mitochondrion", "chloroplast", "plastid", "apicoplast", "kinetoplast",
    "chromoplast", "cyanelle", "leucoplast", "proplastid",
}
# Annotations without RefSeq ``region`` features (Ensembl, and most predictor
# output) carry no ``genome=``, so organelles there are recognised by the
# sequence names the major providers use.  Compared case-insensitively after
# stripping a leading "chr".
ORGANELLE_SEQIDS = {
    "m", "mt", "mtdna", "mito", "mitochondrion", "mitochondrion_genome",
    "pt", "plastid", "chloroplast", "apicoplast", "api",
}
# Ensembl names alt loci and patch scaffolds with these suffixes.  Heuristic,
# unlike the RefSeq attribute rule: use --seqids to override it.
ALT_SEQID_SUFFIXES = ("_patch", "_alt", "_ctg1", "_hap1")


# ---------------------------------------------------------------- GFF3 input


def _open(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _attr(field, key):
    # GFF3 attributes: key=value;key=value.
    for part in field.split(";"):
        part = part.strip()
        if part.startswith(key) and part[len(key):len(key) + 1] == "=":
            return part[len(key) + 1:]
    return None


class Annotation:
    """CDS chains keyed by transcript, with their loci and sequence lengths."""

    def __init__(self):
        self.cds = defaultdict(list)      # tid -> [(start, end), ...]
        self.stops = defaultdict(list)    # tid -> [(start, end), ...] stop_codon
        self.where = {}                   # tid -> (seqid, strand)
        self.gene_of = {}                 # tid -> gene id (may be the tid)
        self.seq_len = {}                 # seqid -> length
        self.region_attrs = {}            # seqid -> attribute string
        self.explicit_genes = False
        # Section 4.  Reference rows that are not a complete protein-coding
        # gene: a pseudogene is truth for "there is no gene here", and a
        # partial CDS has no correct start or stop to hit.
        self.pseudo = set()               # tids with pseudo=true, or a
                                          # pseudogene parent
        self.biotype = {}                 # tid -> gene_biotype, when stated
        self.partial5 = set()             # tids whose CDS 5' end is incomplete
        self.partial3 = set()             # tids whose CDS 3' end is incomplete
        # Transcript ids that appear on more than one sequence or strand.
        # Concatenating per-chromosome predictor output produces exactly this,
        # because AUGUSTUS restarts its gene numbering at g1 in every run.
        # ``load_gff`` keys every chain by (sequence, strand, id) so the
        # collision no longer welds unrelated chains into one transcript; this
        # set records the raw ids that had to be disambiguated, because a
        # submission that needs it is not a valid GFF3 and its author should
        # know.
        self.id_conflicts = set()
        # True when partial ends were read off missing start_codon/stop_codon
        # features rather than off reference range attributes.
        self.partial_from_missing_codon_feature = False

    def chains(self):
        """tid -> (seqid, strand, ((s, e), ...)) with CDS sorted by position."""
        out = {}
        for tid, blocks in self.cds.items():
            seqid, strand = self.where[tid]
            out[tid] = (seqid, strand, tuple(sorted(blocks)))
        return out

    def stop_codon_convention(self):
        """Which stop-codon convention this file actually uses.

        GFF3 says the CDS includes the stop codon and that is what RefSeq,
        Ensembl and the panel references do.  GTF says the opposite, and every
        tool with a GTF lineage inherits it: AUGUSTUS ships
        ``stopCodonExcludedFromCDS true`` in its species parameter files and
        warns about it on stderr, and BRAKER, GeneMark and SNAP output is in
        the same family.  Getting this wrong is not a rounding error -- it
        moves the 3' end of every CDS chain by 3 bp, which is every terminal
        exon, every single-exon gene, every stop codon and every exact
        transcript match -- so it is detected from the file rather than
        trusted to a flag.

        Returns ``(verdict, counts)`` where verdict is ``inside``, ``outside``,
        ``mixed`` or ``unknown``.  ``unknown`` means there were no
        ``stop_codon`` features to judge by, which is the normal case for a
        reference.
        """
        inside = outside = 0
        for tid, stops in self.stops.items():
            blocks = self.cds.get(tid)
            if not blocks or not stops:
                continue
            hit = any(s <= be and bs <= e for (s, e) in stops
                      for (bs, be) in blocks)
            if hit:
                inside += 1
            else:
                outside += 1
        total = inside + outside
        counts = {"transcripts_with_stop_codon_feature": total,
                  "stop_inside_cds": inside, "stop_outside_cds": outside}
        if not total:
            return "unknown", counts
        if outside >= 0.9 * total:
            return "outside", counts
        if inside >= 0.9 * total:
            return "inside", counts
        return "mixed", counts

    def include_stop_codons(self, seq_len=None):
        """Bring a stop-codon-excluding annotation into the GFF3 convention.

        Where a ``stop_codon`` feature exists it is unioned into the CDS
        chain, which is correct even when the stop is split across an intron.
        Where none exists -- a bare GTF-style CDS set -- the last block in the
        direction of translation is extended by 3 bp instead, clipped to the
        sequence, and counted separately so the guess is visible.

        A 3'-partial gene also has no ``stop_codon`` feature and must *not* be
        extended: it stops because the contig does.  Those two cases are the
        same absence and are told apart by the rest of the file -- a
        prediction that uses ``stop_codon`` features at all is stating, by
        omitting one, that the gene is incomplete, which is what
        ``partial3`` holds.  A file with no ``stop_codon`` feature anywhere
        leaves ``partial3`` empty and every chain is extended, as before.
        """
        seq_len = self.seq_len if seq_len is None else seq_len
        merged = extended = not_extended = 0
        for tid, blocks in list(self.cds.items()):
            stops = self.stops.get(tid)
            if stops and not any(s <= be and bs <= e for (s, e) in stops
                                 for (bs, be) in blocks):
                self.cds[tid] = merge_intervals(blocks + stops)
                merged += 1
            elif not stops and tid in self.partial3:
                not_extended += 1
            elif not stops:
                seqid, strand = self.where[tid]
                limit = seq_len.get(seqid)
                b = sorted(blocks)
                if strand == "-":
                    new = (max(1, b[0][0] - 3), b[0][1])
                    b[0] = new
                else:
                    end = b[-1][1] + 3
                    if limit:
                        end = min(end, limit)
                    b[-1] = (b[-1][0], end)
                self.cds[tid] = merge_intervals(b)
                extended += 1
        return {"transcripts_stop_merged_from_feature": merged,
                "transcripts_stop_extended_by_3bp": extended,
                "transcripts_stop_not_extended_partial": not_extended}


def merge_intervals(blocks):
    """Sorted, non-overlapping, adjacency-joined 1-based inclusive blocks."""
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return [tuple(b) for b in out]


def load_gff(path):
    """Parse a GFF3 into CDS chains.

    Tolerant on purpose: predictors emit CDS rows with ``Parent`` pointing at
    an mRNA that may or may not be declared, or with only ``transcript_id``.
    A CDS with no usable parent becomes its own transcript.
    """
    ann = Annotation()
    tx_gene = {}
    gene_biotype = {}   # gene id -> gene_biotype (or "pseudogene" by feature)
    tx_pseudo = set()   # tids marked pseudo on their own transcript row
    max_end = defaultdict(int)
    # A GFF3 id is only unique within the file that declared it, and a
    # prediction assembled by concatenating one AUGUSTUS run per scaffold is
    # 1,158 files with a ``g1.t1`` in each.  A CDS chain can lie on exactly one
    # sequence and one strand, so keying every id by (sequence, strand) is
    # always safe and resolves the collision instead of merging across it.
    # ``raw_where`` keeps the first place each raw id was seen so the file can
    # still be told it was not unique.
    raw_where = {}
    with_start = set()  # tids carrying a start_codon feature
    with_stop = set()   # tids carrying a stop_codon feature

    def _ns(seqid, strand, ident, record=True):
        if ident is None:
            return None
        if record:
            seen = raw_where.get(ident)
            if seen is None:
                raw_where[ident] = (seqid, strand)
            elif seen != (seqid, strand):
                ann.id_conflicts.add(ident)
        return "%s\x00%s\x00%s" % (seqid, strand, ident)
    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                if line.startswith("##sequence-region"):
                    f = line.split()
                    if len(f) == 4:
                        ann.seq_len[f[1]] = int(f[3]) - int(f[2]) + 1
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, ftype, attrs = f[0], f[2], f[8]
            if ftype == "region" and f[3] == "1":
                ann.region_attrs[seqid] = attrs
                ann.seq_len.setdefault(seqid, int(f[4]))
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            if end > max_end[seqid]:
                max_end[seqid] = end
            if ftype in ("gene", "pseudogene"):
                gid = _ns(seqid, f[6], _attr(attrs, "ID")
                          or _attr(attrs, "gene_id"), record=False)
                if gid:
                    # RefSeq writes ``gene_biotype``, GENCODE ``gene_type``,
                    # Ensembl ``biotype``.  All three name the same thing, and
                    # reading only the first silently disables section 4's
                    # filter on the other two: GENCODE 50 states
                    # ``gene_type=IG_V_gene`` on 421 CDS-bearing transcripts
                    # that RefSeq drops as ``V_segment`` and friends, so a
                    # GENCODE-shaped submission was charged 421 false
                    # positives the reference was not allowed to answer.
                    gene_biotype[gid] = (_attr(attrs, "gene_biotype")
                                         or _attr(attrs, "gene_type")
                                         or _attr(attrs, "biotype")
                                         or ("pseudogene" if ftype == "pseudogene"
                                             else None))
            elif ftype not in ("CDS", "stop_codon", "start_codon", "exon"):
                # Anything else with an ID may be the parent of a CDS.  Not
                # only mRNA: RefSeq gives immunoglobulin and T-cell receptor
                # segments their own feature types (``V_gene_segment``,
                # ``C_gene_segment``), and their CDS rows would otherwise have
                # no path to the gene row that says what they are.
                tid = _ns(seqid, f[6], _attr(attrs, "ID")
                          or _attr(attrs, "transcript_id"))
                gid = _ns(seqid, f[6], (_attr(attrs, "Parent")
                                        or _attr(attrs, "gene_id")),
                          record=False)
                if tid:
                    tx_gene[tid] = gid or tid
                    if gid and ftype in TX_TYPES:
                        ann.explicit_genes = True
                    if _attr(attrs, "pseudo") == "true":
                        tx_pseudo.add(tid)
            elif ftype == "CDS":
                tid = (_attr(attrs, "Parent") or _attr(attrs, "transcript_id")
                       or _attr(attrs, "ID"))
                if not tid:
                    tid = "%s:%d:%s" % (seqid, start, f[6])
                # A CDS may list several parents; the first is the transcript.
                tid = _ns(seqid, f[6], tid.split(",")[0])
                ann.cds[tid].append((start, end))
                ann.where[tid] = (seqid, f[6])
                if _attr(attrs, "pseudo") == "true":
                    ann.pseudo.add(tid)
                # RefSeq marks an incomplete CDS end with start_range= or
                # end_range=, in genome coordinates, so which *biological* end
                # is missing depends on the strand.  ``partial=true`` with
                # neither range attribute says only that something is missing;
                # both ends are then treated as unknown.  These attributes sit
                # on the CDS row: on pombe's mRNA rows 788 of 5,166 carry
                # partial=true, describing incomplete UTRs, while only 6 CDS
                # rows do.
                sr = _attr(attrs, "start_range") is not None
                er = _attr(attrs, "end_range") is not None
                if f[6] == "-":
                    sr, er = er, sr
                if sr:
                    ann.partial5.add(tid)
                if er:
                    ann.partial3.add(tid)
                if not sr and not er and _attr(attrs, "partial") == "true":
                    ann.partial5.add(tid)
                    ann.partial3.add(tid)
            elif ftype == "stop_codon":
                tid = (_attr(attrs, "Parent") or _attr(attrs, "transcript_id")
                       or _attr(attrs, "ID"))
                if tid:
                    tid = _ns(seqid, f[6], tid.split(",")[0])
                    ann.stops[tid].append((start, end))
                    with_stop.add(tid)
            elif ftype == "start_codon":
                # Not scored directly -- the start position comes from the CDS
                # chain -- but its *absence* is how a GTF-lineage predictor
                # says a gene runs off the end of a contig.
                tid = (_attr(attrs, "Parent") or _attr(attrs, "transcript_id")
                       or _attr(attrs, "ID"))
                if tid:
                    with_start.add(_ns(seqid, f[6], tid.split(",")[0]))
    # Section 4.5.  A reference states an incomplete CDS end with
    # ``start_range=``/``end_range=``; a *predictor* has no such convention and
    # none writes those attributes, so ``predicted_partial_*`` read 0 on the
    # first prediction that actually contained partial genes (AUGUSTUS
    # ``--genemodel=partial`` on the 1,158 scaffolds of *T. thermophila*: 250
    # transcripts with no start codon, 69 with no stop).  What a GTF-lineage
    # predictor does instead is emit no ``start_codon``/``stop_codon`` feature
    # for the end that runs off the contig, so the omission is the statement,
    # and it is only readable in a file that uses those features at all --
    # Helixer and Tiberius emit neither and must not be read as all-partial.
    # The first CDS block's ``phase`` is the other candidate signal and is
    # weaker: 69 of those 250 truncated genes happen to resume in frame and
    # carry phase 0.
    if with_start:
        ann.partial5.update(tid for tid in ann.cds if tid not in with_start)
        ann.partial_from_missing_codon_feature = True
    if with_stop:
        ann.partial3.update(tid for tid in ann.cds if tid not in with_stop)
        ann.partial_from_missing_codon_feature = True
    for tid in ann.cds:
        gid = tx_gene.get(tid, tid)
        ann.gene_of[tid] = gid
        bt = gene_biotype.get(gid)
        if bt:
            ann.biotype[tid] = bt
        if tid in tx_pseudo or (bt and bt.endswith("pseudogene")):
            ann.pseudo.add(tid)
    for seqid, end in max_end.items():
        ann.seq_len.setdefault(seqid, end)
    return ann


def select_seqids(ref, whitelist=None, min_len=MIN_SEQ_LEN, reasons=None):
    """Sequences the run is scored on (section 4: nuclear primary assembly).

    A whitelist is taken verbatim, intersected with what the reference has.

    Otherwise, in order:

    * drop anything shorter than ``min_len``;
    * where the reference has RefSeq ``region`` features, drop organelles
      (``genome=`` in `ORGANELLE_GENOMES`) and drop alt loci and patch
      scaffolds, which are the ``genome=genomic`` regions that carry a
      cytogenetic ``map=`` band;
    * where it does not (Ensembl, predictor output), fall back to the
      sequence-name heuristics in `ORGANELLE_SEQIDS` and `ALT_SEQID_SUFFIXES`.

    The ``map=`` test, not ``chromosome=``, is what separates an alt locus from
    an unplaced scaffold.  Unlocalized scaffolds carry ``map=unlocalized``, and
    unplaced ones carry ``chromosome=Unknown`` with no ``map=`` at all: in
    GRCh38.p14 all 680 ``genome=genomic`` regions have a ``chromosome=``, and
    so do all 675 of maize's and all 32 of Nematostella's, so keying on
    ``chromosome=`` drops every unplaced scaffold in the panel.

    ``reasons``, if given, is a dict that receives ``seqid -> reason`` for
    every sequence *not* kept.
    """
    if reasons is None:
        reasons = {}
    if whitelist is not None:
        keep = set(whitelist) & set(ref.seq_len)
        for seqid in ref.seq_len:
            if seqid not in keep:
                reasons[seqid] = "not in --seqids"
        return keep
    keep = set()
    for seqid, length in ref.seq_len.items():
        if length < min_len:
            reasons[seqid] = "shorter than %d bp" % min_len
            continue
        attrs = ref.region_attrs.get(seqid)
        if attrs is not None:
            genome = _attr(attrs, "genome")
            if genome in ORGANELLE_GENOMES:
                reasons[seqid] = "organelle (genome=%s)" % genome
                continue
            band = _attr(attrs, "map")
            if genome == "genomic" and band and band != "unlocalized":
                reasons[seqid] = "alt locus or patch (map=%s)" % band
                continue
        else:
            bare = seqid.lower()
            if bare.startswith("chr"):
                bare = bare[3:]
            if bare in ORGANELLE_SEQIDS:
                reasons[seqid] = "organelle (sequence name)"
                continue
            if bare.endswith(ALT_SEQID_SUFFIXES):
                reasons[seqid] = "alt locus or patch (sequence name)"
                continue
        keep.add(seqid)
    return keep


def select_transcripts(chains, ann, keep_all=False, reasons=None):
    """Transcripts the run is scored on (section 4: complete protein-coding).

    A reference CDS row is not automatically a protein-coding gene.  RefSeq
    annotates pseudogenes with real CDS blocks carrying ``pseudo=true`` --
    32 transcripts in *S. pombe*, one of them described as "malic enzyme with
    2 frameshifts", 314 in GRCh38.p14 -- and a predictor that correctly
    declines to call a frameshifted pseudogene is otherwise charged a false
    negative at nucleotide, exon, locus and transcript level.  Immunoglobulin
    and T-cell receptor segments (``gene_biotype=V_segment`` and friends, 891
    transcripts in GRCh38.p14) are the other case: they are gene *fragments*
    with no start or stop codon of their own.

    Both are dropped by default and counted by reason.  The filter is applied
    to the prediction as well as the reference, so an identity run still
    scores 1.0, and it only ever acts on evidence the file states: a
    prediction without ``pseudo=`` or ``gene_biotype=`` loses nothing.

    ``--score-all-transcripts`` restores the old behaviour for auditing.

    ``reasons``, if given, is a dict that receives ``tid -> reason``.
    """
    if reasons is None:
        reasons = {}
    if keep_all:
        return dict(chains)
    keep = {}
    for tid, chain in chains.items():
        if tid in ann.pseudo:
            reasons[tid] = "pseudogene"
            continue
        bt = ann.biotype.get(tid)
        if bt is not None and bt != "protein_coding":
            reasons[tid] = "biotype=%s" % bt
            continue
        keep[tid] = chain
    return keep


# ------------------------------------------------------------ interval maths


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + 1:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def total(intervals):
    return sum(e - s + 1 for s, e in intervals)


def isect_len(a, b):
    """Length of the intersection of two merged, sorted interval lists."""
    i = j = n = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e >= s:
            n += e - s + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def prf(tp, fp, fn):
    sens = tp / (tp + fn) if tp + fn else None
    prec = tp / (tp + fp) if tp + fp else None
    if sens and prec:
        f1 = 2 * sens * prec / (sens + prec)
    else:
        f1 = 0.0 if (tp + fp + fn) else None
    return {"tp": tp, "fp": fp, "fn": fn,
            "sensitivity": _r(sens), "precision": _r(prec), "f1": _r(f1)}


def _r(x, nd=5):
    return None if x is None else round(x, nd)


# ------------------------------------------------------------------- metrics


def cds_intervals_by_strand(chains, seqids):
    by = defaultdict(list)
    for seqid, strand, blocks in chains.values():
        if seqid in seqids:
            by[(seqid, strand)].extend(blocks)
    return {k: merge(v) for k, v in by.items()}


def nucleotide(ref_chains, pred_chains, seqids, scored_bp):
    """Section 4.1.  Strand-aware; overlapping genes unioned; MCC reported."""
    r = cds_intervals_by_strand(ref_chains, seqids)
    p = cds_intervals_by_strand(pred_chains, seqids)
    tp = 0
    for key in set(r) & set(p):
        tp += isect_len(r[key], p[key])
    ref_bp = sum(total(v) for v in r.values())
    pred_bp = sum(total(v) for v in p.values())
    fp, fn = pred_bp - tp, ref_bp - tp
    out = prf(tp, fp, fn)
    tn = 2 * scored_bp - tp - fp - fn      # two strands
    denom = math.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out["mcc"] = _r((tp * tn - fp * fn) / denom) if denom > 0 else None
    out["reference_cds_bp"] = ref_bp
    out["predicted_cds_bp"] = pred_bp
    out["scored_bp"] = scored_bp
    return out


def exon_type(i, n):
    if n == 1:
        return "single"
    if i == 0:
        return "initial"
    if i == n - 1:
        return "terminal"
    return "internal"


TYPE_RANK = {"single": 0, "initial": 1, "terminal": 2, "internal": 3}


def exon_sets(chains, seqids):
    """(seqid, start, end, strand) -> exon type, over all transcripts.

    An exon shared by isoforms in different roles is labelled with the
    highest-priority role it plays (single > initial > terminal > internal),
    so that stratification is a partition.
    """
    out = {}
    for seqid, strand, blocks in chains.values():
        if seqid not in seqids:
            continue
        n = len(blocks)
        order = blocks if strand == "+" else tuple(reversed(blocks))
        for i, (s, e) in enumerate(order):
            key = (seqid, s, e, strand)
            t = exon_type(i, n)
            if key not in out or TYPE_RANK[t] < TYPE_RANK[out[key]]:
                out[key] = t
    return out


def exons(ref_chains, pred_chains, seqids):
    """Section 4.2.  Both boundaries exact; stratified by exon type."""
    ref = exon_sets(ref_chains, seqids)
    pred = exon_sets(pred_chains, seqids)
    tp_keys = set(ref) & set(pred)
    out = {"all": prf(len(tp_keys), len(pred) - len(tp_keys), len(ref) - len(tp_keys))}
    by_type = {}
    for t in ("initial", "internal", "terminal", "single"):
        r = {k for k, v in ref.items() if v == t}
        p = {k for k, v in pred.items() if v == t}
        tp = len(r & p)
        by_type[t] = prf(tp, len(p) - tp, len(r) - tp)
    out["by_type"] = by_type

    # "Roughly right, structurally wrong": predicted exons that overlap an
    # annotated exon but match neither boundary (section 4.2).
    starts = defaultdict(list)
    for seqid, s, e, strand in ref:
        starts[(seqid, strand)].append((s, e))
    for key in starts:
        starts[key].sort()
    overlap_only = 0
    for key in pred:
        seqid, s, e, strand = key
        if key in ref:
            continue
        arr = starts.get((seqid, strand))
        if not arr:
            continue
        i = bisect.bisect_right(arr, (e, float("inf")))
        hit = False
        # Reference exons are short relative to the genome; walking back a
        # bounded number of entries finds any overlap without an interval tree.
        j = i - 1
        while j >= 0 and arr[j][0] >= s - 5_000_000:
            rs, re_ = arr[j]
            if re_ >= s and rs <= e and rs != s and re_ != e:
                hit = True
                break
            j -= 1
        if hit:
            overlap_only += 1
    out["predicted_overlap_no_boundary"] = overlap_only
    out["predicted_overlap_no_boundary_frac"] = _r(
        overlap_only / len(pred) if pred else None)
    return out


def introns_of(chains, seqids, min_len=MIN_INTRON):
    """(seqid, start, end, strand) of every CDS-chain intron, de-duplicated.

    Returns ``(introns, short_gaps)``: gaps shorter than ``min_len`` are not
    splice junctions (see ``MIN_INTRON``) and are returned separately so the
    count is reported rather than silently dropped.
    """
    out, short = set(), set()
    for seqid, strand, blocks in chains.values():
        if seqid not in seqids:
            continue
        for i in range(len(blocks) - 1):
            s, e = blocks[i][1] + 1, blocks[i + 1][0] - 1
            if e < s:
                continue
            (out if e - s + 1 >= min_len else short).add((seqid, s, e, strand))
    return out, short


def sites_of(introns):
    """Donor and acceptor positions, strand-aware, with their intron length.

    One position can carry several intron lengths: alternative splicing shares
    a donor between a short and a long intron.  On the panel's references that
    is 0.00% of *S. cerevisiae* donors but 12.35% of human ones, with a spread
    of over 1 Mb between the shortest and the longest intron at one site, so
    which length is kept decides that site's §4.3 decile.  ``introns`` is a
    set, so taking whichever came last made the stratification depend on the
    interpreter's hash seed: the same file scored twice gave different
    per-decile counts.  The shortest intron at the site is kept, which is
    arbitrary but fixed; ``sites_multiple_intron_lengths`` in the result says
    how many sites the choice was made for.
    """
    donors, acceptors = {}, {}
    multi = [set(), set()]
    for seqid, s, e, strand in introns:
        length = e - s + 1
        d, a = (s, e) if strand == "+" else (e, s)
        for i, (m, key) in enumerate(((donors, (seqid, d, strand)),
                                      (acceptors, (seqid, a, strand)))):
            prev = m.get(key)
            if prev is None:
                m[key] = length
            else:
                if prev != length:
                    multi[i].add(key)
                m[key] = min(prev, length)
    return donors, acceptors, (len(multi[0]), len(multi[1]))


def deciles(values):
    v = sorted(values)
    if not v:
        return []
    return [v[min(len(v) - 1, int(round(q / 10 * (len(v) - 1))))] for q in range(1, 10)]


def splice(ref_chains, pred_chains, seqids, seqfetch=None):
    """Section 4.3.  Donors and acceptors separately, stratified."""
    ref_i, ref_short = introns_of(ref_chains, seqids)
    pred_i, pred_short = introns_of(pred_chains, seqids)
    cuts = deciles([e - s + 1 for _, s, e, _ in ref_i])

    def bucket(length):
        return bisect.bisect_left(cuts, length) if cuts else 0

    out = {"reference_introns": len(ref_i), "predicted_introns": len(pred_i),
           "reference_cds_gaps_below_min": len(ref_short),
           "predicted_cds_gaps_below_min": len(pred_short),
           "min_intron": MIN_INTRON,
           "intron_length_decile_cuts": cuts}
    ref_d, ref_a, ref_multi = sites_of(ref_i)
    pred_d, pred_a, pred_multi = sites_of(pred_i)
    for i, (name, r, p) in enumerate((("donor", ref_d, pred_d),
                                      ("acceptor", ref_a, pred_a))):
        tp = set(r) & set(p)
        out[name] = prf(len(tp), len(p) - len(tp), len(r) - len(tp))
        # Sites whose decile was decided by the shortest-intron rule of
        # ``sites_of``; the totals above do not depend on it, the
        # stratification below does.
        out[name]["reference_sites_multiple_intron_lengths"] = ref_multi[i]
        out[name]["predicted_sites_multiple_intron_lengths"] = pred_multi[i]
        strat = []
        for b in range(10):
            rb = {k for k in r if bucket(r[k]) == b}
            pb = {k for k in p if bucket(p[k]) == b}
            t = len(rb & pb)
            strat.append(prf(t, len(pb) - t, len(rb) - t))
        out[name]["by_intron_length_decile"] = strat

    if seqfetch is not None:
        out["by_dinucleotide"] = _splice_by_dinuc(ref_i, pred_i, seqfetch)
        out["by_local_gc"] = _splice_by_gc(ref_d, pred_d, seqfetch)
    else:
        out["by_dinucleotide"] = None
        out["by_local_gc"] = None
    return out


COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _revcomp(s):
    return s.translate(COMP)[::-1]


def dinuc_class(seqfetch, seqid, s, e, strand):
    d = seqfetch(seqid, s, s + 1)
    a = seqfetch(seqid, e - 1, e)
    if d is None or a is None:
        return "unknown"
    if strand == "-":
        d, a = _revcomp(a), _revcomp(d)
    pair = (d + "-" + a).upper()
    return pair if pair in ("GT-AG", "GC-AG", "AT-AC") else "other"


def _splice_by_dinuc(ref_i, pred_i, seqfetch):
    classes = defaultdict(lambda: [0, 0, 0])  # tp, fp, fn
    pred = set(pred_i)
    for iv in ref_i:
        c = dinuc_class(seqfetch, *iv)
        if iv in pred:
            classes[c][0] += 1
        else:
            classes[c][2] += 1
    for iv in pred - ref_i:
        classes[dinuc_class(seqfetch, *iv)][1] += 1
    return {c: prf(*v) for c, v in sorted(classes.items())}


def gc_bin(seqfetch, seqid, pos):
    half = GC_WINDOW // 2
    s = seqfetch(seqid, max(1, pos - half), pos + half - 1)
    if not s:
        return None
    s = s.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    if acgt == 0:
        return None
    gc = (s.count("G") + s.count("C")) / acgt
    return bisect.bisect_left(GC_BINS, gc)


def _splice_by_gc(ref_d, pred_d, seqfetch):
    bins = defaultdict(lambda: [0, 0, 0])
    pred = set(pred_d)
    for k in ref_d:
        b = gc_bin(seqfetch, k[0], k[1])
        if b is None:
            continue
        bins[b][0 if k in pred else 2] += 1
    for k in pred - set(ref_d):
        b = gc_bin(seqfetch, k[0], k[1])
        if b is not None:
            bins[b][1] += 1
    return {str(b): prf(*v) for b, v in sorted(bins.items())}


def loci_of(ann, chains, seqids):
    """gene id -> (seqid, strand, merged CDS, [transcript ids]).

    Uses the GFF3 gene grouping when the file has one; otherwise clusters
    transcripts whose CDS overlaps on the same strand, so a predictor that
    emits no gene features is still scored at the locus level.
    """
    groups = defaultdict(list)
    for tid, (seqid, strand, blocks) in chains.items():
        if seqid in seqids:
            groups[ann.gene_of.get(tid, tid)].append(tid)
    if not ann.explicit_genes:
        groups = _cluster_by_overlap(chains, seqids)
    out = {}
    for gid, tids in groups.items():
        seqid, strand, _ = chains[tids[0]]
        blocks = []
        for tid in tids:
            blocks.extend(chains[tid][2])
        out[gid] = (seqid, strand, merge(blocks), tids)
    return out


def _cluster_by_overlap(chains, seqids):
    by = defaultdict(list)
    for tid, (seqid, strand, blocks) in chains.items():
        if seqid in seqids:
            by[(seqid, strand)].append((blocks[0][0], blocks[-1][1], tid))
    groups = {}
    n = 0
    for key, items in by.items():
        items.sort()
        cur, cur_end = [], -1
        for s, e, tid in items:
            if cur and s > cur_end:
                n += 1
                groups["cluster%d" % n] = cur
                cur, cur_end = [], -1
            cur.append(tid)
            cur_end = max(cur_end, e)
        if cur:
            n += 1
            groups["cluster%d" % n] = cur
    return groups


def _overlap_pairs(ref_loci, pred_loci):
    """Overlapping (ref, pred) locus pairs with their shared CDS bases."""
    index = defaultdict(list)
    for pid, (seqid, strand, blocks, _) in pred_loci.items():
        index[(seqid, strand)].append((blocks[0][0], blocks[-1][1], pid, blocks))
    for key in index:
        index[key].sort()
    pairs = []
    for rid, (seqid, strand, rblocks, _) in ref_loci.items():
        arr = index.get((seqid, strand))
        if not arr:
            continue
        rs, re = rblocks[0][0], rblocks[-1][1]
        i = bisect.bisect_right(arr, (re, float("inf"), "", ()))
        j = i - 1
        while j >= 0:
            ps, pe, pid, pblocks = arr[j]
            if pe >= rs:
                n = isect_len(rblocks, pblocks)
                if n > 0:
                    pairs.append((n, rid, pid))
            if ps < rs - 10_000_000:
                break
            j -= 1
    return pairs


def _components(loci):
    """locus id -> component id of the same-strand CDS-overlap graph.

    Two annotated genes whose CDS overlap are not separable at CDS level in
    that annotation, so no prediction should be blamed for touching both.
    Transitivity matters and pairwise tests are not enough: in GRCh38.p14 there
    are 82 places where gene B overlaps both A and C while A and C are disjoint,
    and a pairwise rule scores a *perfect* prediction of B as a fusion of A and
    C.  Grouping by connected component makes an identity run score exactly
    zero fusions and zero splits, on every panel species.
    """
    parent = {lid: lid for lid in loci}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    by = defaultdict(list)
    for lid, (seqid, strand, blocks, _) in loci.items():
        by[(seqid, strand)].append((blocks[0][0], blocks[-1][1], lid, blocks))
    for items in by.values():
        items.sort()
        for i, (s, e, lid, blocks) in enumerate(items):
            for j in range(i + 1, len(items)):
                if items[j][0] > e:
                    break
                if isect_len(blocks, items[j][3]) > 0:
                    union(lid, items[j][2])
    return {lid: find(lid) for lid in loci}


def _overlapping_pairs(loci):
    """Same-strand CDS-overlapping locus pairs, a property of the reference."""
    by = defaultdict(list)
    for lid, (seqid, strand, blocks, _) in loci.items():
        by[(seqid, strand)].append((blocks[0][0], blocks[-1][1], lid, blocks))
    n = 0
    for items in by.values():
        items.sort()
        for i, (s, e, lid, blocks) in enumerate(items):
            for j in range(i + 1, len(items)):
                if items[j][0] > e:
                    break
                if isect_len(blocks, items[j][3]) > 0:
                    n += 1
    return n


def loci(ref_ann, pred_ann, ref_chains, pred_chains, seqids):
    """Section 4.4: locus level with fusion and split counts, then transcripts."""
    ref_loci = loci_of(ref_ann, ref_chains, seqids)
    pred_loci = loci_of(pred_ann, pred_chains, seqids)
    pairs = _overlap_pairs(ref_loci, pred_loci)

    per_ref, per_pred = defaultdict(set), defaultdict(set)
    for _, rid, pid in pairs:
        per_ref[rid].add(pid)
        per_pred[pid].add(rid)
    # A prediction that overlaps two annotated genes is only a fusion if the
    # reference itself keeps those genes apart.  Yeast alone has 91 same-strand
    # CDS-overlapping RefSeq gene pairs and a correct prediction of one member
    # necessarily touches the other, so overlap is charged to the annotation,
    # not to the predictor: fusions are counted across connected components of
    # the reference's own overlap graph.  Splits are the mirror image, over the
    # prediction's overlap graph.
    ref_comp = _components(ref_loci)
    pred_comp = _components(pred_loci)
    fusion = sum(1 for pid, rids in per_pred.items()
                 if len({ref_comp[r] for r in rids}) > 1)
    split = sum(1 for rid, pids in per_ref.items()
                if len({pred_comp[p] for p in pids}) > 1)

    # One-to-one, greedy by shared CDS bases: each reference locus keeps the
    # prediction that overlaps it best, and no prediction is used twice.
    matched_r, matched_p, matches = set(), set(), []
    for n, rid, pid in sorted(pairs, key=lambda x: (-x[0], x[1], x[2])):
        if rid in matched_r or pid in matched_p:
            continue
        matched_r.add(rid)
        matched_p.add(pid)
        matches.append((rid, pid))
    locus = prf(len(matches), len(pred_loci) - len(matches), len(ref_loci) - len(matches))
    locus["fusion"] = fusion
    locus["split"] = split
    locus["reference_loci"] = len(ref_loci)
    locus["predicted_loci"] = len(pred_loci)
    locus["reference_overlapping_locus_pairs"] = _overlapping_pairs(ref_loci)

    tx = _transcripts(ref_loci, pred_loci, ref_chains, pred_chains, matches)
    tx["unmatched_reference_loci"] = len(ref_loci) - len(matches)
    return locus, tx


def _transcripts(ref_loci, pred_loci, ref_chains, pred_chains, matches):
    """Section 4.4 transcript exact match, under the isoform rule.

    Within a matched locus, predicted transcripts are matched one-to-one to
    annotated isoforms, exact chain matches first and the rest by shared CDS
    bases; unmatched annotated isoforms are not counted as false negatives, so
    deeply annotated species are not penalized against one-transcript-per-gene
    species.

    Exactness has to outrank shared bases, not just break ties inside it.  An
    isoform that contains the prediction's whole CDS shares exactly as many
    bases with it as the isoform the prediction *equals*, so ordering on
    shared bases alone leaves the winner to the id tie-break -- and the
    accession that sorts first is not the one that matches.  Fugu is the
    demonstration: the prediction of `rab44` is `rna-XM_029826178.1` base for
    base, but `rna-XM_011613896.2` is that chain plus 51 more bases at one
    end, so both share 10,749, `011` sorts before `029`, and an exact hit is
    scored as a miss.  343 of 5,907 exact fugu matches (5.8%) were lost that
    way.  Ranking exactness first is also the only order that keeps the metric
    a property of the two annotations rather than of their accession strings.
    """
    tp = fp = fn = 0
    for rid, pid in matches:
        r_tids = ref_loci[rid][3]
        p_tids = pred_loci[pid][3]
        cand = []
        for pt in p_tids:
            pb = pred_chains[pt][2]
            for rt in r_tids:
                n = isect_len(list(pb), list(ref_chains[rt][2]))
                if n > 0:
                    exact = ref_chains[rt] == pred_chains[pt]
                    cand.append((0 if exact else 1, -n, rt, pt))
        used_r, used_p = set(), set()
        for _, _, rt, pt in sorted(cand):
            if rt in used_r or pt in used_p:
                continue
            used_r.add(rt)
            used_p.add(pt)
            if ref_chains[rt] == pred_chains[pt]:
                tp += 1
            else:
                fp += 1
                fn += 1
        fp += len(p_tids) - len(used_p)
    fn += len(ref_loci) - len(matches)
    for pid, (_, _, _, p_tids) in pred_loci.items():
        if pid not in {p for _, p in matches}:
            fp += len(p_tids)
    return prf(tp, fp, fn)


def codons(ref_chains, pred_chains, seqids, stop_in_cds=True,
           detected="unknown", detected_counts=None, merge_stats=None,
           genome_detected="unknown", genome_counts=None, source="assumed",
           ref_ann=None, pred_ann=None):
    """Section 4.5.  Start and stop codon positions, reported separately.

    A CDS annotated as incomplete at one end contributes no codon at that end.
    RefSeq states this with ``start_range=``/``end_range=`` on the CDS row:
    1,415 human and 194 maize transcripts are partial this way.  A 5'-partial
    gene has no annotated start codon to hit, so counting its first base as a
    reference start codon charges every predictor a false negative it cannot
    avoid; the 3' end is the mirror image.  The two ends are excluded
    independently.

    Dropping the reference position alone would only move the charge: the
    predictor's own start, somewhere inside that gene, would then be a false
    positive.  So a predicted position inside the span of a reference
    transcript that is partial at that end is dropped from the denominator
    too -- unless it coincides with a reference position that survived, which
    is a real match and stays a true positive.
    """
    r5 = getattr(ref_ann, "partial5", set()) if ref_ann else set()
    r3 = getattr(ref_ann, "partial3", set()) if ref_ann else set()
    p5 = getattr(pred_ann, "partial5", set()) if pred_ann else set()
    p3 = getattr(pred_ann, "partial3", set()) if pred_ann else set()

    def ends(chains, skip5, skip3):
        starts, stops = set(), set()
        for tid, (seqid, strand, blocks) in chains.items():
            if seqid not in seqids or not blocks:
                continue
            if tid in skip5 and tid in skip3:
                continue
            if strand == "+":
                first, last = blocks[0][0], blocks[-1][1]
            else:
                first, last = blocks[-1][1], blocks[0][0]
            if tid not in skip5:
                starts.add((seqid, first, strand))
            if tid not in skip3:
                stops.add((seqid, last, strand))
        return starts, stops

    def spans(tids):
        """(seqid, strand) -> sorted list of (lo, hi) CDS spans."""
        out = defaultdict(list)
        for tid in tids:
            chain = ref_chains.get(tid)
            if not chain or not chain[2]:
                continue
            seqid, strand, blocks = chain
            out[(seqid, strand)].append((blocks[0][0], blocks[-1][1]))
        # Merged, so one bisect answers "is this position inside a span".
        return {k: merge_intervals(v) for k, v in out.items()}

    def outside(positions, unscorable, keep):
        """Predicted positions not inside an unscorable reference span."""
        out = set()
        for pos in positions:
            if pos in keep:
                out.add(pos)
                continue
            seqid, coord, strand = pos
            spanlist = unscorable.get((seqid, strand))
            if not spanlist:
                out.add(pos)
                continue
            i = bisect.bisect_right(spanlist, (coord, float("inf"))) - 1
            if i < 0 or spanlist[i][1] < coord:
                out.add(pos)
        return out

    r_start, r_stop = ends(ref_chains, r5, r3)
    p_start, p_stop = ends(pred_chains, p5, p3)
    p_start = outside(p_start, spans(r5), r_start)
    p_stop = outside(p_stop, spans(r3), r_stop)
    out = {}
    for name, r, p in (("start", r_start, p_start), ("stop", r_stop, p_stop)):
        tp = len(r & p)
        out[name] = prf(tp, len(p) - tp, len(r) - tp)
    # What the run was told, what the prediction file actually shows, and
    # what was done about it.
    out["reference_partial_5prime"] = len(r5 & set(ref_chains))
    out["reference_partial_3prime"] = len(r3 & set(ref_chains))
    out["predicted_partial_5prime"] = len(p5 & set(pred_chains))
    out["predicted_partial_3prime"] = len(p3 & set(pred_chains))
    # Where those two numbers came from: ``range_attribute`` is the reference
    # convention (``start_range=``/``end_range=``), ``missing_codon_feature``
    # is the predictor one (no ``start_codon``/``stop_codon`` row for the end
    # that runs off the contig) -- which reads that way with both counts zero
    # when the file uses codon features and declares every gene complete --
    # and ``none`` means the file cannot state it at all, so every chain end
    # was taken as real.
    out["predicted_partial_source"] = (
        "missing_codon_feature"
        if getattr(pred_ann, "partial_from_missing_codon_feature", False)
        else ("range_attribute" if (p5 or p3) else "none"))
    out["stop_codon_inside_cds"] = bool(stop_in_cds)
    out["stop_codon_convention_detected"] = detected
    # Where the convention actually came from: a stop_codon feature, the
    # genome, the --stop-outside-cds flag, or nothing at all.  "assumed" means
    # the run had no evidence either way and took the GFF3 default; it is the
    # value to look for before believing a terminal-exon or stop-codon score.
    out["stop_codon_convention_from_genome"] = genome_detected
    out["stop_codon_convention_source"] = source
    out.update(detected_counts or {})
    out.update(genome_counts or {})
    out.update(merge_stats or {})
    return out


# ---------------------------------------------------- optional genome access


class WindowFetcher:
    """One streaming pass over a FASTA, serving a fixed set of windows.

    Random access into a 3 Gb gzipped FASTA is not worth an index here: every
    window the scorer needs is known before the pass starts, so the sequence
    is streamed once with a rolling buffer no larger than the gap between
    consecutive requested windows.
    """

    def __init__(self, path, windows):
        self.cache = {}
        want = defaultdict(list)
        for seqid, s, e in windows:
            want[seqid].append((s, e))
        for seqid in want:
            want[seqid] = sorted(set(want[seqid]))
        self._fill(path, want)

    def _fill(self, path, want):
        with _open(path) as fh:
            seqid = None
            pending = deque()
            buf, buf_start, pos = [], 1, 0
            for line in fh:
                if line.startswith(">"):
                    self._flush(seqid, pending, buf, buf_start, pos)
                    seqid = line[1:].split()[0]
                    pending = deque(want.get(seqid, []))
                    buf, buf_start, pos = [], 1, 0
                    continue
                if seqid is None or not pending:
                    continue
                seq = line.strip()
                if not seq:
                    continue
                buf.append(seq)
                pos += len(seq)
                joined = None
                while pending and pending[0][1] <= pos:
                    s, e = pending.popleft()
                    if joined is None:
                        joined = "".join(buf)
                    if s >= buf_start:
                        self.cache[(seqid, s, e)] = joined[s - buf_start:e - buf_start + 1]
                if pending:
                    # Never trim past what has actually been read: buf_start
                    # must stay equal to pos - len(buffer) + 1.
                    keep = min(pending[0][0], pos + 1)
                    if keep > buf_start:
                        if joined is None:
                            joined = "".join(buf)
                        buf = [joined[keep - buf_start:]]
                        buf_start = keep
            self._flush(seqid, pending, buf, buf_start, pos)

    def _flush(self, seqid, pending, buf, buf_start, pos):
        """Serve the windows still queued when a FASTA record ends.

        A window whose end runs past the end of the sequence -- a local-GC
        window on a splice site within ``GC_WINDOW // 2`` of a scaffold end --
        is never satisfied by the ``pending[0][1] <= pos`` test above.  Because
        ``pending`` is ordered by window *start*, that unsatisfiable window
        also sits in front of every later window on the same sequence and
        blocks all of them: before this flush, one such intron on a 43 kb
        *T. rubripes* scaffold silently cost both its acceptor dinucleotide and
        its GC bin, and on a more fragmented assembly it would cost every
        window behind it too.  Windows are dropped, not clipped, at the *start*
        of a sequence by ``max(1, ...)`` in ``needed_windows``; the end needs
        this.  The buffer is never trimmed past ``pending[0][0]``, so every
        queued window is still fully in ``buf`` here, and the sequence is
        truncated rather than the window discarded.
        """
        if seqid is None or not pending:
            return
        joined = "".join(buf)
        while pending:
            s, e = pending.popleft()
            if s >= buf_start and s <= pos:
                self.cache[(seqid, s, e)] = joined[s - buf_start:e - buf_start + 1]

    def __call__(self, seqid, s, e):
        return self.cache.get((seqid, s, e))


def needed_windows(ref_chains, pred_chains, seqids):
    """Every window the dinucleotide and GC stratifications will ask for."""
    wins = set()
    half = GC_WINDOW // 2
    for chains in (ref_chains, pred_chains):
        for seqid, s, e, strand in introns_of(chains, seqids)[0]:
            wins.add((seqid, s, s + 1))
            wins.add((seqid, e - 1, e))
            d = s if strand == "+" else e
            wins.add((seqid, max(1, d - half), d + half - 1))
    return wins


def _convention_probes(chains, seqids):
    """(tid, inside window, outside window) for the stop-convention check.

    ``inside`` is the last codon of the CDS chain in the direction of
    translation; ``outside`` is the codon immediately after it on the genome.
    A chain whose terminal block is shorter than 3 bp is skipped rather than
    walked back across the intron: the split-codon case is rare and a wrong
    window would be counted as evidence.
    """
    out = []
    for tid, (seqid, strand, blocks) in chains.items():
        if seqid not in seqids or not blocks:
            continue
        if strand == "-":
            s = blocks[0][0]
            if blocks[0][1] - s + 1 < 3 or s < 4:
                continue
            out.append((tid, seqid, strand, (s, s + 2), (s - 3, s - 1)))
        else:
            e = blocks[-1][1]
            if e - blocks[-1][0] + 1 < 3:
                continue
            out.append((tid, seqid, strand, (e - 2, e), (e + 1, e + 3)))
    return out


def stop_convention_windows(chains, seqids):
    wins = set()
    for _, seqid, _, ins, outs in _convention_probes(chains, seqids):
        wins.add((seqid, ins[0], ins[1]))
        wins.add((seqid, outs[0], outs[1]))
    return wins


def stop_convention_from_genome(chains, seqids, fetch, stops):
    """Which stop-codon convention the prediction uses, read off the genome.

    ``Annotation.stop_codon_convention`` can only answer when the file carries
    ``stop_codon`` features.  Helixer and Tiberius carry none -- their GFF3 is
    gene/mRNA/exon/CDS/UTR and nothing else -- so for them the answer was
    previously the *default of a flag*, and being wrong about it moves the 3'
    end of every chain by 3 bp while leaving nucleotide and locus numbers
    looking healthy.  With ``--genome`` the question is decidable without any
    feature: either the last codon of the chain is a stop or the codon just
    past it is.

    Returns ``(verdict, counts)`` with the same vocabulary as the feature-based
    check.  ``unknown`` means too few usable chains, or neither position looks
    like a stop often enough to call it.
    """
    stops = set(stops)
    inside = outside = neither = 0
    for _, seqid, strand, ins, outs in _convention_probes(chains, seqids):
        a = fetch(seqid, ins[0], ins[1])
        b = fetch(seqid, outs[0], outs[1])
        if a is None or b is None or len(a) != 3 or len(b) != 3:
            continue
        if strand == "-":
            a, b = _revcomp(a), _revcomp(b)
        a, b = a.upper(), b.upper()
        if a in stops:
            inside += 1
        elif b in stops:
            outside += 1
        else:
            neither += 1
    total = inside + outside + neither
    counts = {"stop_convention_genome_chains": total,
              "stop_convention_genome_inside": inside,
              "stop_convention_genome_outside": outside,
              "stop_convention_genome_neither": neither}
    if total < MIN_CONVENTION_CHAINS:
        return "unknown", counts
    if inside >= 0.9 * total:
        return "inside", counts
    if outside >= 0.9 * total:
        return "outside", counts
    return "unknown", counts


# ----------------------------------------------------------- the declaration


def parse_declaration(path):
    """Minimal YAML: top-level ``key: value`` lines, ``#`` comments.

    Deliberately not a YAML parser.  It checks that every key section 3.3
    requires is present and non-empty and records a hash of the file, so the
    scored run is traceable to the declaration it was scored under.
    """
    raw = open(path, "rb").read()
    seen = {}
    for line in raw.decode("utf-8", "replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#") or line[0] in " \t-":
            continue
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        seen[k.strip()] = v.strip()
    missing = [k for k in DECLARATION_KEYS if not seen.get(k)]
    return seen, missing, hashlib.sha256(raw).hexdigest()


# -------------------------------------------------------------------- driver


def selection_summary(ref, keep, dropped):
    """What section 4's sequence filter did, so a scored run is auditable.

    Reasons are grouped, not listed per sequence: a human run drops 638
    sequences for one reason and naming them would swamp the report.
    """
    by_reason = defaultdict(lambda: [0, 0])
    for seqid, reason in dropped.items():
        kind = reason.split(" (")[0]
        by_reason[kind][0] += 1
        by_reason[kind][1] += ref.seq_len.get(seqid, 0)
    return {
        "kept": len(keep),
        "dropped": len(dropped),
        "dropped_by_reason": {k: {"sequences": v[0], "bp": v[1]}
                              for k, v in sorted(by_reason.items())},
        "reference_has_region_features": bool(ref.region_attrs),
    }


def transcript_selection_summary(dropped, kept):
    by_reason = defaultdict(int)
    for reason in dropped.values():
        by_reason[reason] += 1
    return {"kept": kept, "dropped": len(dropped),
            "dropped_by_reason": dict(sorted(by_reason.items()))}


def score(reference, prediction, species, genome=None, seqids=None,
          stop_in_cds=True, all_transcripts=False, genetic_code=1):
    ref = load_gff(reference)
    pred = load_gff(prediction)
    # Section 4.5.  Detected from the file, not taken on trust: a prediction
    # whose CDS excludes the stop codon is 3 bp short at the 3' end of every
    # chain, which silently zeroes terminal exons, single-exon genes, stop
    # codons and exact transcript matches while leaving the nucleotide and
    # locus numbers looking healthy.
    pred_convention, convention_counts = pred.stop_codon_convention()
    dropped = {}
    keep = select_seqids(ref, whitelist=seqids, reasons=dropped)
    ref_chains = {t: c for t, c in ref.chains().items() if c[0] in keep}
    pred_chains = {t: c for t, c in pred.chains().items() if c[0] in keep}

    # One pass over the FASTA serves both the splice strata and the
    # convention probes.  The probes are planned from the seqid-filtered
    # chains rather than the transcript-filtered ones, which is a superset and
    # costs a few unused windows; planning them after the 3 bp extension is
    # not possible, because the extension is what the probes decide.  Introns
    # are unaffected by that extension -- it moves the 3' end of the terminal
    # block and adds no junction -- so the splice windows are the same either
    # way.
    fetch = None
    planned = set()
    genome_convention, genome_counts = "unknown", {}
    if genome:
        planned = (needed_windows(ref_chains, pred_chains, keep)
                   | stop_convention_windows(pred_chains, keep))
        fetch = WindowFetcher(genome, planned)
        genome_convention, genome_counts = stop_convention_from_genome(
            pred_chains, keep, fetch, STOP_CODONS[genetic_code])

    merge_stats = {"transcripts_stop_merged_from_feature": 0,
                   "transcripts_stop_extended_by_3bp": 0,
                   "transcripts_stop_not_extended_partial": 0}
    convention_source = "assumed"
    if not stop_in_cds:
        merge_stats = pred.include_stop_codons(seq_len=ref.seq_len)
        convention_source = "flag"
    elif pred_convention == "unknown" and genome_convention == "outside":
        # No stop_codon feature to read and the genome says the stop sits
        # just past the chain.  Taking the default here would score every
        # terminal exon, single-exon gene, stop codon and exact transcript
        # match 3 bp short, so the evidence wins over the default.  It never
        # overrides an explicit --stop-outside-cds; that case is the branch
        # above.
        merge_stats = pred.include_stop_codons(seq_len=ref.seq_len)
        convention_source = "genome"
    elif genome_convention != "unknown":
        convention_source = "genome"
    elif pred_convention != "unknown":
        convention_source = "stop_codon_feature"
    if any(merge_stats.values()):
        pred_chains = {t: c for t, c in pred.chains().items() if c[0] in keep}
        # The comment above is right that the 3 bp extension adds no junction,
        # but merging a ``stop_codon`` feature can: AUGUSTUS splits a stop
        # codon across an intron, and unioning that feature into the chain
        # creates a junction the window plan never saw.  One predicted intron
        # in the *S. pombe* cross-parameter run is exactly this, and without a
        # second pass its dinucleotide was reported as `unknown` rather than
        # the GT-AG it is.  A second pass runs only for the windows the first
        # one did not plan, so the common case still reads the FASTA once.
        if fetch is not None:
            extra = needed_windows(ref_chains, pred_chains, keep) - planned
            if extra:
                fetch.cache.update(WindowFetcher(genome, extra).cache)
                planned |= extra
    # Section 4: pseudogenes and gene fragments are not scored as truth, and
    # the same filter is applied to the prediction so that an identity run is
    # still exactly 1.0.
    ref_dropped_tx, pred_dropped_tx = {}, {}
    ref_chains = select_transcripts(ref_chains, ref, all_transcripts,
                                    ref_dropped_tx)
    pred_chains = select_transcripts(pred_chains, pred, all_transcripts,
                                     pred_dropped_tx)
    scored_bp = sum(ref.seq_len[s] for s in keep)

    locus, tx = loci(ref, pred, ref_chains, pred_chains, keep)
    return {
        "species": species,
        "reference": os.path.basename(reference),
        "prediction": os.path.basename(prediction),
        "scored_sequences": len(keep),
        "scored_bp": scored_bp,
        "sequence_selection": selection_summary(ref, keep, dropped),
        "reference_transcripts": len(ref_chains),
        "predicted_transcripts": len(pred_chains),
        # Reference rows that carry a CDS but are not a complete
        # protein-coding gene, grouped by why they were dropped.
        "reference_transcript_selection":
            transcript_selection_summary(ref_dropped_tx, len(ref_chains)),
        "predicted_transcript_selection":
            transcript_selection_summary(pred_dropped_tx, len(pred_chains)),
        # Predictions on sequences the filter excluded, or on sequence names
        # the reference does not have at all, are not scored.  A large count
        # here means the prediction was made against a different assembly or
        # a different naming convention and the run is not comparable.
        # Rows the section 4 biotype filter dropped are *not* counted here --
        # they are on a scored sequence and are reported, by reason, in
        # ``predicted_transcript_selection`` -- because pooling the two makes
        # a GENCODE-shaped submission look like an assembly mismatch.
        "predicted_transcripts_not_scored":
            len(pred.chains()) - len(pred_chains) - len(pred_dropped_tx),
        "predicted_sequences_absent_from_reference":
            len({c[0] for c in pred.chains().values()} - set(ref.seq_len)),
        # Non-zero means the prediction reuses a transcript id across
        # sequences or strands, so those chains are a merge of unrelated
        # genes and nothing downstream of here is meaningful for them.
        "predicted_conflicting_transcript_ids": len(pred.id_conflicts),
        "reference_conflicting_transcript_ids": len(ref.id_conflicts),
        "nucleotide": nucleotide(ref_chains, pred_chains, keep, scored_bp),
        "exon": exons(ref_chains, pred_chains, keep),
        "splice": splice(ref_chains, pred_chains, keep, fetch),
        "transcript": tx,
        "locus": locus,
        "codon": codons(ref_chains, pred_chains, keep, stop_in_cds,
                        pred_convention, convention_counts, merge_stats,
                        genome_convention, genome_counts, convention_source,
                        ref, pred),
    }


# ----------------------------------------------------------------- self-test

REF_FIXTURE = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tregion\t1\t20000\t.\t+\t.\tID=chr1;genome=chromosome
chr1\t.\tgene\t1000\t3000\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t1000\t3000\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tID=c1;Parent=t1
chr1\t.\tCDS\t2000\t2100\t.\t+\t0\tID=c2;Parent=t1
chr1\t.\tCDS\t2900\t3000\t.\t+\t0\tID=c3;Parent=t1
chr1\t.\tgene\t6000\t6500\t.\t-\t.\tID=g2
chr1\t.\tmRNA\t6000\t6500\t.\t-\t.\tID=t2;Parent=g2
chr1\t.\tCDS\t6000\t6200\t.\t-\t0\tID=c4;Parent=t2
chr1\t.\tCDS\t6400\t6500\t.\t-\t0\tID=c5;Parent=t2
chr1\t.\tgene\t9000\t9300\t.\t+\t.\tID=g3
chr1\t.\tmRNA\t9000\t9300\t.\t+\t.\tID=t3;Parent=g3
chr1\t.\tCDS\t9000\t9300\t.\t+\t0\tID=c6;Parent=t3
"""

# Four predicted loci exercising four outcomes: pg1 reproduces t1 exactly;
# pg2 shifts one boundary of t2's minus-strand CDS, which moves the acceptor
# and not the donor; pg4 overlaps g3 sharing neither boundary; pg3 overlaps
# nothing annotated.
PRED_FIXTURE = """##gff-version 3
chr1\t.\tgene\t1000\t3000\t.\t+\t.\tID=pg1
chr1\t.\tmRNA\t1000\t3000\t.\t+\t.\tID=p1;Parent=pg1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2000\t2100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2900\t3000\t.\t+\t0\tParent=p1
chr1\t.\tgene\t6000\t6500\t.\t-\t.\tID=pg2
chr1\t.\tmRNA\t6000\t6500\t.\t-\t.\tID=p2;Parent=pg2
chr1\t.\tCDS\t6000\t6250\t.\t-\t0\tParent=p2
chr1\t.\tCDS\t6400\t6500\t.\t-\t0\tParent=p2
chr1\t.\tgene\t9050\t9250\t.\t+\t.\tID=pg4
chr1\t.\tmRNA\t9050\t9250\t.\t+\t.\tID=p4;Parent=pg4
chr1\t.\tCDS\t9050\t9250\t.\t+\t0\tParent=p4
chr1\t.\tgene\t15000\t15200\t.\t+\t.\tID=pg3
chr1\t.\tmRNA\t15000\t15200\t.\t+\t.\tID=p3;Parent=pg3
chr1\t.\tCDS\t15000\t15200\t.\t+\t0\tParent=p3
"""


# Two annotated neighbours predicted as one locus (a fusion), and one
# annotated gene predicted as two (a split).
FUSION_REF = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tgene\t1000\t1500\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t1000\t1500\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t1000\t1500\t.\t+\t0\tParent=t1
chr1\t.\tgene\t1700\t2200\t.\t+\t.\tID=g2
chr1\t.\tmRNA\t1700\t2200\t.\t+\t.\tID=t2;Parent=g2
chr1\t.\tCDS\t1700\t2200\t.\t+\t0\tParent=t2
chr1\t.\tgene\t5000\t6000\t.\t+\t.\tID=g3
chr1\t.\tmRNA\t5000\t6000\t.\t+\t.\tID=t3;Parent=g3
chr1\t.\tCDS\t5000\t6000\t.\t+\t0\tParent=t3
"""

FUSION_PRED = """##gff-version 3
chr1\t.\tgene\t1000\t2200\t.\t+\t.\tID=pg1
chr1\t.\tmRNA\t1000\t2200\t.\t+\t.\tID=p1;Parent=pg1
chr1\t.\tCDS\t1000\t1500\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t1700\t2200\t.\t+\t0\tParent=p1
chr1\t.\tgene\t5000\t5400\t.\t+\t.\tID=pg2
chr1\t.\tmRNA\t5000\t5400\t.\t+\t.\tID=p2;Parent=pg2
chr1\t.\tCDS\t5000\t5400\t.\t+\t0\tParent=p2
chr1\t.\tgene\t5600\t6000\t.\t+\t.\tID=pg3
chr1\t.\tmRNA\t5600\t6000\t.\t+\t.\tID=p3;Parent=pg3
chr1\t.\tCDS\t5600\t6000\t.\t+\t0\tParent=p3
"""


# The same three transcripts as REF_FIXTURE, written the way AUGUSTUS writes
# them: every CDS chain stops 3 bp short and the stop codon is its own
# feature.  On the minus strand the missing 3 bp are at the low coordinate.
# With --stop-outside-cds this must score exactly 1.0 everywhere; without it,
# nothing whose 3' end is a stop codon can match.
# A GTF-lineage predictor run with ``--genemodel=partial`` on a fragmented
# assembly: complete genes carry both codon features, a gene truncated at the
# contig start carries no ``start_codon``, and one truncated at the contig end
# carries no ``stop_codon``.  That omission is the only statement of
# partiality such a file makes -- no predictor writes RefSeq's
# ``start_range=``/``end_range=``.
PARTIAL_PRED = """##gff-version 3
chr1\t.\tgene\t1000\t3000\t.\t+\t.\tID=pg1
chr1\t.\ttranscript\t1000\t3000\t.\t+\t.\tID=p1;Parent=pg1
chr1\t.\tstart_codon\t1000\t1002\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2000\t2100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2900\t2997\t.\t+\t0\tParent=p1
chr1\t.\tstop_codon\t2998\t3000\t.\t+\t0\tParent=p1
chr1\t.\tgene\t6000\t6450\t.\t-\t.\tID=pg2
chr1\t.\ttranscript\t6000\t6450\t.\t-\t.\tID=p2;Parent=pg2
chr1\t.\tstop_codon\t6000\t6002\t.\t-\t0\tParent=p2
chr1\t.\tCDS\t6003\t6200\t.\t-\t0\tParent=p2
chr1\t.\tCDS\t6400\t6450\t.\t-\t0\tParent=p2
chr1\t.\tgene\t9000\t9297\t.\t+\t.\tID=pg3
chr1\t.\ttranscript\t9000\t9297\t.\t+\t.\tID=p3;Parent=pg3
chr1\t.\tstart_codon\t9000\t9002\t.\t+\t0\tParent=p3
chr1\t.\tCDS\t9000\t9297\t.\t+\t0\tParent=p3
"""

STOP_PRED = """##gff-version 3
chr1\t.\tgene\t1000\t3000\t.\t+\t.\tID=pg1
chr1\t.\ttranscript\t1000\t3000\t.\t+\t.\tID=p1;Parent=pg1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2000\t2100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2900\t2997\t.\t+\t0\tParent=p1
chr1\t.\tstop_codon\t2998\t3000\t.\t+\t0\tParent=p1
chr1\t.\tgene\t6000\t6500\t.\t-\t.\tID=pg2
chr1\t.\ttranscript\t6000\t6500\t.\t-\t.\tID=p2;Parent=pg2
chr1\t.\tstop_codon\t6000\t6002\t.\t-\t0\tParent=p2
chr1\t.\tCDS\t6003\t6200\t.\t-\t0\tParent=p2
chr1\t.\tCDS\t6400\t6500\t.\t-\t0\tParent=p2
chr1\t.\tgene\t9000\t9300\t.\t+\t.\tID=pg3
chr1\t.\ttranscript\t9000\t9300\t.\t+\t.\tID=p3;Parent=pg3
chr1\t.\tCDS\t9000\t9297\t.\t+\t0\tParent=p3
chr1\t.\tstop_codon\t9298\t9300\t.\t+\t0\tParent=p3
"""

# The same again with the stop_codon rows removed, which is what a GTF-derived
# prediction looks like after a naive conversion: the convention can no longer
# be detected and the 3 bp have to be guessed back.
# A prediction whose ``stop_codon`` feature sits on the far side of an intron
# from its last CDS block: unioning it into the chain creates a junction that
# did not exist in the file, so the window plan has to be extended after the
# merge rather than before it.  AUGUSTUS emits this whenever a stop codon is
# split by an intron.
SPLIT_STOP_REF = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tregion\t1\t20000\t.\t+\t.\tID=chr1;genome=chromosome
chr1\t.\tgene\t1000\t2502\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t1000\t2502\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tID=c1;Parent=t1
chr1\t.\tCDS\t2000\t2502\t.\t+\t0\tID=c2;Parent=t1
"""

SPLIT_STOP_PRED = """##gff-version 3
chr1\t.\tgene\t1000\t2502\t.\t+\t.\tID=pg1
chr1\t.\ttranscript\t1000\t2502\t.\t+\t.\tID=p1;Parent=pg1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2000\t2101\t.\t+\t0\tParent=p1
chr1\t.\tstop_codon\t2500\t2502\t.\t+\t0\tParent=p1
"""


STOP_PRED_BARE = "\n".join(
    l for l in STOP_PRED.splitlines() if "\tstop_codon\t" not in l) + "\n"


# Section 4 transcript selection, in the shapes RefSeq actually writes:
# a normal gene; a pseudogene whose CDS carries pseudo=true under a
# gene_biotype=pseudogene parent; an immunoglobulin V segment; a 5'-partial
# CDS on the plus strand (start_range=); and a partial CDS on the minus
# strand whose start_range= is its *3'* end.
PSEUDO_REF = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tregion\t1\t20000\t.\t+\t.\tID=chr1;genome=chromosome
chr1\t.\tgene\t1000\t2100\t.\t+\t.\tID=g1;gene_biotype=protein_coding
chr1\t.\tmRNA\t1000\t2100\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tID=c1;Parent=t1
chr1\t.\tCDS\t2000\t2100\t.\t+\t0\tID=c2;Parent=t1
chr1\t.\tpseudogene\t4000\t4300\t.\t+\t.\tID=g2;gene_biotype=pseudogene
chr1\t.\tmRNA\t4000\t4300\t.\t+\t.\tID=t2;Parent=g2
chr1\t.\tCDS\t4000\t4300\t.\t+\t0\tID=c3;Parent=t2;pseudo=true
chr1\t.\tgene\t5000\t5200\t.\t+\t.\tID=g3;gene_biotype=V_segment
chr1\t.\tV_gene_segment\t5000\t5200\t.\t+\t.\tID=t3;Parent=g3
chr1\t.\tCDS\t5000\t5200\t.\t+\t0\tID=c4;Parent=t3
chr1\t.\tgene\t7000\t7300\t.\t+\t.\tID=g4;gene_biotype=protein_coding
chr1\t.\tmRNA\t7000\t7300\t.\t+\t.\tID=t4;Parent=g4
chr1\t.\tCDS\t7000\t7300\t.\t+\t0\tID=c5;Parent=t4;partial=true;start_range=.,7000
chr1\t.\tgene\t9000\t9300\t.\t-\t.\tID=g5;gene_biotype=protein_coding
chr1\t.\tmRNA\t9000\t9300\t.\t-\t.\tID=t5;Parent=g5
chr1\t.\tCDS\t9000\t9300\t.\t-\t0\tID=c6;Parent=t5;partial=true;start_range=.,9000
"""

# A predictor's view of the same locus set: it calls the two complete genes
# exactly, calls each partial one with an end it invented 60 bp inside the
# annotated span -- the 5' end of t4 and, on the minus strand, the 3' end of
# t5 -- and does not call the pseudogene or the V segment at all.
PSEUDO_PRED = """##gff-version 3
chr1\t.\tgene\t1000\t2100\t.\t+\t.\tID=pg1
chr1\t.\tmRNA\t1000\t2100\t.\t+\t.\tID=p1;Parent=pg1
chr1\t.\tCDS\t1000\t1100\t.\t+\t0\tParent=p1
chr1\t.\tCDS\t2000\t2100\t.\t+\t0\tParent=p1
chr1\t.\tgene\t7060\t7300\t.\t+\t.\tID=pg2
chr1\t.\tmRNA\t7060\t7300\t.\t+\t.\tID=p2;Parent=pg2
chr1\t.\tCDS\t7060\t7300\t.\t+\t0\tParent=p2
chr1\t.\tgene\t9060\t9300\t.\t-\t.\tID=pg3
chr1\t.\tmRNA\t9060\t9300\t.\t-\t.\tID=p3;Parent=pg3
chr1\t.\tCDS\t9060\t9300\t.\t-\t0\tParent=p3
"""

# The same prediction as PSEUDO_PRED, plus the V segment and a pseudogene,
# declared the way GENCODE (``gene_type``) and Ensembl (``biotype``) declare
# them rather than the way RefSeq does (``gene_biotype``).  Section 4 applies
# the same filter to the prediction as to the reference, so all three
# spellings must be read; otherwise a GENCODE-shaped submission is charged
# false positives for the very rows the reference is forbidden to score.
PSEUDO_PRED_SYNONYM = PSEUDO_PRED + (
    "chr1\t.\tgene\t5000\t5200\t.\t+\t.\tID=pg4;gene_type=IG_V_gene\n"
    "chr1\t.\tmRNA\t5000\t5200\t.\t+\t.\tID=p4;Parent=pg4\n"
    "chr1\t.\tCDS\t5000\t5200\t.\t+\t0\tParent=p4\n"
    "chr1\t.\tgene\t4000\t4300\t.\t+\t.\tID=pg5;biotype=pseudogene\n"
    "chr1\t.\tmRNA\t4000\t4300\t.\t+\t.\tID=p5;Parent=pg5\n"
    "chr1\t.\tCDS\t4000\t4300\t.\t+\t0\tParent=p5\n")


# Section 4 sequence selection.  Every line is a shape taken from a real panel
# reference: a chromosome, an alt locus (map= band), an unlocalized scaffold
# (map=unlocalized), an unplaced scaffold (chromosome=Unknown, no map=, the
# shape maize and Nematostella use), a Tetrahymena-style scaffold with neither,
# a mitochondrion, a chloroplast, and a scaffold under the 10 kb floor.
SHARED_DONOR_REF = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tregion\t1\t20000\t.\t+\t.\tID=chr1;genome=chromosome
chr1\t.\tgene\t1000\t9000\t.\t+\t.\tID=g1
chr1\t.\tmRNA\t1000\t2200\t.\t+\t.\tID=t1;Parent=g1
chr1\t.\tCDS\t1000\t1099\t.\t+\t0\tID=c1;Parent=t1
chr1\t.\tCDS\t1200\t2200\t.\t+\t2\tID=c2;Parent=t1
chr1\t.\tmRNA\t1000\t9000\t.\t+\t.\tID=t2;Parent=g1
chr1\t.\tCDS\t1000\t1099\t.\t+\t0\tID=c3;Parent=t2
chr1\t.\tCDS\t8000\t9000\t.\t+\t2\tID=c4;Parent=t2
"""


# The fugu case in miniature (see ``_transcripts``): the prediction equals
# isoform ``rna-t9`` base for base, while ``rna-t1`` is that same chain with
# 51 extra bases at the 3' end, so both share every one of the prediction's
# bases and the accession that sorts first is the one that does not match.
CONTAINED_ISOFORM_REF = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tregion\t1\t20000\t.\t+\t.\tID=chr1;genome=chromosome
chr1\t.\tgene\t1000\t2251\t.\t+\t.\tID=g1;gene_biotype=protein_coding
chr1\t.\tmRNA\t1000\t2251\t.\t+\t.\tID=rna-t1;Parent=g1
chr1\t.\tCDS\t1000\t1099\t.\t+\t0\tID=c1;Parent=rna-t1
chr1\t.\tCDS\t1200\t2251\t.\t+\t2\tID=c2;Parent=rna-t1
chr1\t.\tmRNA\t1000\t2200\t.\t+\t.\tID=rna-t9;Parent=g1
chr1\t.\tCDS\t1000\t1099\t.\t+\t0\tID=c3;Parent=rna-t9
chr1\t.\tCDS\t1200\t2200\t.\t+\t2\tID=c4;Parent=rna-t9
"""


CONTAINED_ISOFORM_PRED = """##gff-version 3
##sequence-region chr1 1 20000
chr1\t.\tgene\t1000\t2200\t.\t+\t.\tID=pg1
chr1\t.\tmRNA\t1000\t2200\t.\t+\t.\tID=pt1;Parent=pg1
chr1\t.\tCDS\t1000\t1099\t.\t+\t0\tID=pc1;Parent=pt1
chr1\t.\tCDS\t1200\t2200\t.\t+\t2\tID=pc2;Parent=pt1
"""


SELECT_FIXTURE = """##gff-version 3
chr1\t.\tregion\t1\t248956422\t.\t+\t.\tID=r1;chromosome=1;genome=chromosome
alt1\t.\tregion\t1\t200000\t.\t+\t.\tID=r2;chromosome=1;genome=genomic;map=19q13.42
unloc1\t.\tregion\t1\t175055\t.\t+\t.\tID=r3;chromosome=1;genome=genomic;map=unlocalized
unplaced1\t.\tregion\t1\t109197\t.\t+\t.\tID=r4;chromosome=Unknown;genome=genomic
bare1\t.\tregion\t1\t56000\t.\t+\t.\tID=r5;genome=genomic
chrM\t.\tregion\t1\t16569\t.\t+\t.\tID=r6;genome=mitochondrion
chrC\t.\tregion\t1\t154478\t.\t+\t.\tID=r7;genome=chloroplast
tiny1\t.\tregion\t1\t1564\t.\t+\t.\tID=r8;genome=genomic
"""

# The same assembly as Ensembl writes it: no region features at all, so only
# the sequence-name fallback can act.
SELECT_FIXTURE_ENSEMBL = """##gff-version 3
##sequence-region 1 1 248956422
##sequence-region MT 1 16569
##sequence-region MtDNA 1 13794
##sequence-region Pt 1 154478
##sequence-region HG1012_PATCH 1 200000
##sequence-region GL000195.1 1 182896
##sequence-region KI270713.1 1 40745
"""


def select_test(check):
    d = tempfile.mkdtemp(prefix="score-select-")
    try:
        for name, text, want, want_reasons in (
            ("refseq", SELECT_FIXTURE,
             {"chr1", "unloc1", "unplaced1", "bare1"},
             {"organelle": 2, "alt locus or patch": 1,
              "shorter than 10000 bp": 1}),
            ("ensembl", SELECT_FIXTURE_ENSEMBL,
             {"1", "GL000195.1", "KI270713.1"},
             {"organelle": 3, "alt locus or patch": 1}),
        ):
            path = os.path.join(d, name + ".gff3")
            open(path, "w").write(text)
            ann = load_gff(path)
            reasons = {}
            got = select_seqids(ann, reasons=reasons)
            check("select %s kept" % name, got, want)
            counts = defaultdict(int)
            for r in reasons.values():
                counts[r.split(" (")[0]] += 1
            check("select %s reasons" % name, dict(counts), want_reasons)
            os.unlink(path)
        # An explicit whitelist wins over every rule, including the organelle
        # and length ones: it is the documented escape hatch.
        path = os.path.join(d, "w.gff3")
        open(path, "w").write(SELECT_FIXTURE)
        ann = load_gff(path)
        check("select whitelist", select_seqids(ann, whitelist=["chrM", "tiny1"]),
              {"chrM", "tiny1"})
        os.unlink(path)
    finally:
        os.rmdir(d)


def self_test():
    d = tempfile.mkdtemp(prefix="score-selftest-")
    ref = os.path.join(d, "ref.gff3")
    pred = os.path.join(d, "pred.gff3")
    open(ref, "w").write(REF_FIXTURE)
    open(pred, "w").write(PRED_FIXTURE)
    r = score(ref, pred, "fixture")

    fails = []

    n_checks = [0]

    def check(name, got, want):
        n_checks[0] += 1
        if got != want:
            fails.append("%s: got %r, want %r" % (name, got, want))

    # Nucleotide: reference CDS is 303 + 302 + 301 = 906 bp, prediction is
    # 303 + 352 + 201 + 201 = 1057 bp.  t1 is exact, all 302 of t2's bases are
    # covered by a prediction 50 bp too long, and 201 of g3's 301 are hit.
    check("nt tp", r["nucleotide"]["tp"], 303 + 302 + 201)
    check("nt fn", r["nucleotide"]["fn"], 50 + 50)
    check("nt fp", r["nucleotide"]["fp"], 50 + 201)
    check("scored_bp", r["scored_bp"], 20000)
    # Exons: 3 of t1 plus the unshifted t2 exon match; 6 reference, 7 predicted.
    check("exon tp", r["exon"]["all"]["tp"], 4)
    check("exon fp", r["exon"]["all"]["fp"], 3)
    check("exon fn", r["exon"]["all"]["fn"], 2)
    check("single exons", r["exon"]["by_type"]["single"]["fn"], 1)
    # Only pg4 overlaps a reference exon while matching neither boundary; the
    # shifted t2 exon still shares its start and does not count here.
    check("overlap-no-boundary", r["exon"]["predicted_overlap_no_boundary"], 1)
    # Splice sites: t1's two introns are exact.  t2 is on the minus strand, so
    # its donor is the intron's high coordinate, which the shift left alone,
    # and its acceptor is the low one, which the shift moved.
    check("donor tp", r["splice"]["donor"]["tp"], 3)
    check("acceptor tp", r["splice"]["acceptor"]["tp"], 2)
    check("ref introns", r["splice"]["reference_introns"], 3)
    check("short gaps", r["splice"]["reference_cds_gaps_below_min"], 0)
    # Loci: all three annotated genes are hit, pg3 is spurious, and a
    # one-to-one matching means no fusion and no split.
    check("locus tp", r["locus"]["tp"], 3)
    check("locus fn", r["locus"]["fn"], 0)
    check("locus fp", r["locus"]["fp"], 1)
    check("fusion", r["locus"]["fusion"], 0)
    check("split", r["locus"]["split"], 0)
    # Transcripts: only t1 is an exact chain match; the two inexact matches are
    # each both a false positive and a false negative, and pg3 adds one more FP.
    check("tx tp", r["transcript"]["tp"], 1)
    check("tx fp", r["transcript"]["fp"], 3)
    check("tx fn", r["transcript"]["fn"], 2)
    # Codons: t1 and t2 have both right; pg4 misplaces g3's start and stop.
    check("start tp", r["codon"]["start"]["tp"], 2)
    check("stop tp", r["codon"]["stop"]["tp"], 2)

    # Fusion and split are counted, and neither inflates the locus TP count:
    # pg1 covers g1 and g2 (one fusion), pg2 and pg3 both sit in g3 (one
    # split), so three annotated loci get two matches and one FP is left over.
    fref = os.path.join(d, "fref.gff3")
    fpred = os.path.join(d, "fpred.gff3")
    open(fref, "w").write(FUSION_REF)
    open(fpred, "w").write(FUSION_PRED)
    fr = score(fref, fpred, "fixture")
    check("fusion count", fr["locus"]["fusion"], 1)
    check("split count", fr["locus"]["split"], 1)
    check("fusion locus tp", fr["locus"]["tp"], 2)
    check("fusion locus fn", fr["locus"]["fn"], 1)
    check("fusion locus fp", fr["locus"]["fp"], 1)
    os.unlink(fref)
    os.unlink(fpred)

    # A prediction identical to the reference must score 1.0 everywhere.
    perfect = score(ref, ref, "fixture")
    for path in (("nucleotide", "f1"), ("exon", "all", "f1"),
                 ("splice", "donor", "f1"), ("transcript", "f1"),
                 ("locus", "f1"), ("codon", "start", "f1")):
        v = perfect
        for k in path:
            v = v[k]
        if v != 1.0:
            fails.append("self-comparison %s: got %r, want 1.0" % ("/".join(path), v))

    # Section 4.5, the stop-codon convention.  AUGUSTUS-shaped output must
    # score 1.0 everywhere with --stop-outside-cds and must be *detected* as
    # such whether or not the flag is given, because without the flag the
    # metrics anchored on the 3' end all collapse while nucleotide and locus
    # stay high enough to look like a real result.
    spred = os.path.join(d, "spred.gff3")
    open(spred, "w").write(STOP_PRED)
    sr = score(ref, spred, "fixture", stop_in_cds=False)
    check("stop detected", sr["codon"]["stop_codon_convention_detected"], "outside")
    check("stop merged", sr["codon"]["transcripts_stop_merged_from_feature"], 3)
    check("stop extended", sr["codon"]["transcripts_stop_extended_by_3bp"], 0)
    for path in (("nucleotide", "f1"), ("exon", "all", "f1"),
                 ("transcript", "f1"), ("locus", "f1"),
                 ("codon", "start", "f1"), ("codon", "stop", "f1")):
        v = sr
        for k in path:
            v = v[k]
        if v != 1.0:
            fails.append("stop-outside %s: got %r, want 1.0" % ("/".join(path), v))
    # The same file scored without the flag: still detected, and the damage is
    # exactly the terminal and single exons and every stop codon.
    sw = score(ref, spred, "fixture")
    check("stop detected without flag",
          sw["codon"]["stop_codon_convention_detected"], "outside")
    check("stop tp without flag", sw["codon"]["stop"]["tp"], 0)
    # t1's two 5' exons and t2's initial exon survive: on the minus strand it
    # is the *first* block that carries the stop codon.
    check("stop exon tp without flag", sw["exon"]["all"]["tp"], 3)
    check("stop tx tp without flag", sw["transcript"]["tp"], 0)
    check("stop start tp without flag", sw["codon"]["start"]["tp"], 3)
    os.unlink(spred)

    # No stop_codon features at all: the convention cannot be detected, the
    # 3 bp are guessed, and the guess is reported.
    bpred = os.path.join(d, "bpred.gff3")
    open(bpred, "w").write(STOP_PRED_BARE)
    br = score(ref, bpred, "fixture", stop_in_cds=False)
    check("bare stop detected",
          br["codon"]["stop_codon_convention_detected"], "unknown")
    check("bare stop extended",
          br["codon"]["transcripts_stop_extended_by_3bp"], 3)
    check("bare stop tp", br["codon"]["stop"]["tp"], 3)
    check("bare exon f1", br["exon"]["all"]["f1"], 1.0)
    check("bare source", br["codon"]["stop_codon_convention_source"], "flag")
    # ... and without the flag it is scored 3 bp short and says only that it
    # assumed, which is the state Helixer and Tiberius output was in before
    # the genome check below.
    bw = score(ref, bpred, "fixture")
    check("bare no flag source",
          bw["codon"]["stop_codon_convention_source"], "assumed")
    check("bare no flag stop tp", bw["codon"]["stop"]["tp"], 0)

    # The same file with --genome: no stop_codon feature exists, but the
    # sequence settles it.  chr1 is filled with C so that no filler codon and
    # no reverse complement of one is a stop, and a stop is placed at exactly
    # the three positions REF_FIXTURE's chains end on: 2998-3000 for t1,
    # 6000-6002 (TTA, read TAA on the minus strand) for t2, and 9298-9300 for
    # t3.  STOP_PRED_BARE stops 3 bp short of each, so its stops are outside.
    gfa = os.path.join(d, "conv.fa")
    cseq = list("C" * 20000)
    for start, codon in ((2998, "TAA"), (6000, "TTA"), (9298, "TAA")):
        cseq[start - 1:start + 2] = list(codon)
    cseq = "".join(cseq)
    with open(gfa, "w") as fh:
        fh.write(">chr1 convention fixture\n")
        for i in range(0, len(cseq), 60):
            fh.write(cseq[i:i + 60] + "\n")
    bg = score(ref, bpred, "fixture", genome=gfa)
    check("genome convention outside",
          bg["codon"]["stop_codon_convention_from_genome"], "outside")
    check("genome convention source",
          bg["codon"]["stop_codon_convention_source"], "genome")
    check("genome chains", bg["codon"]["stop_convention_genome_chains"], 3)
    check("genome outside count",
          bg["codon"]["stop_convention_genome_outside"], 3)
    check("genome inside count",
          bg["codon"]["stop_convention_genome_inside"], 0)
    # The payoff: the 3 bp are put back without the flag, so the run that was
    # scoring 0.00 stop-codon F1 now scores the same as the flagged one.
    check("genome extended", bg["codon"]["transcripts_stop_extended_by_3bp"], 3)
    check("genome stop tp", bg["codon"]["stop"]["tp"], 3)
    check("genome exon f1", bg["exon"]["all"]["f1"], 1.0)
    check("genome tx f1", bg["transcript"]["f1"], 1.0)
    # A prediction that already includes its stop codons is read as such and
    # is left alone.
    ig = score(ref, ref, "fixture", genome=gfa)
    check("genome convention inside",
          ig["codon"]["stop_codon_convention_from_genome"], "inside")
    check("genome inside not extended",
          ig["codon"]["transcripts_stop_extended_by_3bp"], 0)
    check("genome inside f1", ig["codon"]["stop"]["f1"], 1.0)
    # The translation table is consulted rather than hard-coded: under table 6
    # (Tetrahymena) TAA is glutamine, so the same three codons are no longer
    # evidence and the scorer declines to guess instead of extending wrongly.
    t6 = score(ref, bpred, "fixture", genome=gfa, genetic_code=6)
    check("code 6 convention",
          t6["codon"]["stop_codon_convention_from_genome"], "unknown")
    check("code 6 neither", t6["codon"]["stop_convention_genome_neither"], 3)
    check("code 6 not extended",
          t6["codon"]["transcripts_stop_extended_by_3bp"], 0)
    os.unlink(gfa)
    os.unlink(bpred)

    # A reference has no stop_codon features either, so an ordinary run must
    # report "unknown" and change nothing.
    check("ref convention", r["codon"]["stop_codon_convention_detected"], "unknown")
    check("ref not merged", r["codon"]["transcripts_stop_merged_from_feature"], 0)

    # An id reused across sequences must be counted *and* resolved.  This is
    # not hypothetical: AUGUSTUS restarts its gene numbering at g1 in every
    # run, so the only way to parallelise it over a fragmented assembly --
    # one run per scaffold, concatenated -- produces a ``g1.t1`` per scaffold.
    # Merging those chains does not merely lose the colliding transcripts, it
    # welds one chain per sequence into a single chain spanning the genome,
    # and every metric collapses.  The fixture therefore checks the fix, not
    # the count: the colliding file must score exactly what the same
    # prediction scores after ``--uniqueGeneId=true`` would have renamed it.
    cref = os.path.join(d, "cref.gff3")
    open(cref, "w").write(REF_FIXTURE + REF_FIXTURE
                          .replace("##gff-version 3\n", "")
                          .replace("chr1", "chr2")
                          .replace("=t", "=u").replace("=g", "=h")
                          .replace("=c", "=d"))
    cpred = os.path.join(d, "cpred.gff3")
    open(cpred, "w").write(PRED_FIXTURE + PRED_FIXTURE
                           .replace("##gff-version 3\n", "")
                           .replace("chr1", "chr2"))
    upred = os.path.join(d, "upred.gff3")
    open(upred, "w").write(PRED_FIXTURE + PRED_FIXTURE
                           .replace("##gff-version 3\n", "")
                           .replace("chr1", "chr2")
                           .replace("=p", "=q").replace("=pg", "=qg"))
    cr = score(cref, cpred, "fixture")
    ur = score(cref, upred, "fixture")
    check("id conflicts", cr["predicted_conflicting_transcript_ids"], 4)
    check("ref id conflicts", cr["reference_conflicting_transcript_ids"], 0)
    check("unique-id run has no conflicts",
          ur["predicted_conflicting_transcript_ids"], 0)
    check("both sequences scored", cr["scored_sequences"], 2)
    check("colliding chains not merged", cr["predicted_transcripts"], 8)
    check("colliding == renamed transcripts",
          cr["predicted_transcripts"], ur["predicted_transcripts"])
    for block in ("nucleotide", "exon", "splice", "transcript", "locus",
                  "codon"):
        if cr[block] != ur[block]:
            fails.append("colliding ids change %s: %r != %r"
                         % (block, cr[block], ur[block]))
    check("colliding nucleotide f1", cr["nucleotide"]["f1"],
          ur["nucleotide"]["f1"])
    os.unlink(cpred)
    os.unlink(upred)
    os.unlink(cref)

    # Section 4.5, a prediction that is itself partial.  RefSeq's range
    # attributes are a *reference* convention; a predictor states an
    # incomplete end by omitting the codon feature, and if that is not read
    # the truncated end is charged as a wrong codon and, worse, the blind 3 bp
    # extension puts a stop codon onto a gene that ran off the contig -- which
    # can turn a miss into a spurious true positive.
    apred = os.path.join(d, "apred.gff3")
    open(apred, "w").write(PARTIAL_PRED)
    ar = score(ref, apred, "fixture", stop_in_cds=False)
    check("partial pred 5prime", ar["codon"]["predicted_partial_5prime"], 1)
    check("partial pred 3prime", ar["codon"]["predicted_partial_3prime"], 1)
    check("partial pred source", ar["codon"]["predicted_partial_source"],
          "missing_codon_feature")
    check("partial merged", ar["codon"]["transcripts_stop_merged_from_feature"],
          2)
    check("partial not extended",
          ar["codon"]["transcripts_stop_extended_by_3bp"], 0)
    check("partial extension skipped",
          ar["codon"]["transcripts_stop_not_extended_partial"], 1)
    check("partial start tp", ar["codon"]["start"]["tp"], 2)
    check("partial start fp", ar["codon"]["start"]["fp"], 0)
    check("partial start fn", ar["codon"]["start"]["fn"], 1)
    check("partial stop tp", ar["codon"]["stop"]["tp"], 2)
    check("partial stop fp", ar["codon"]["stop"]["fp"], 0)
    check("partial stop fn", ar["codon"]["stop"]["fn"], 1)
    os.unlink(apred)

    # Section 4 transcript selection: pseudogenes and gene fragments are not
    # truth, and the ends of a partial CDS are not scorable at either side.
    pref = os.path.join(d, "pref.gff3")
    ppred = os.path.join(d, "ppred.gff3")
    open(pref, "w").write(PSEUDO_REF)
    open(ppred, "w").write(PSEUDO_PRED)
    pr = score(pref, ppred, "fixture")
    check("pseudo ref kept", pr["reference_transcripts"], 3)
    check("pseudo ref dropped",
          pr["reference_transcript_selection"]["dropped_by_reason"],
          {"biotype=V_segment": 1, "pseudogene": 1})
    check("pseudo pred dropped",
          pr["predicted_transcript_selection"]["dropped"], 0)
    # Neither the pseudogene nor the V segment is a false negative now.
    check("pseudo locus fn", pr["locus"]["fn"], 0)
    check("pseudo locus tp", pr["locus"]["tp"], 3)
    # t4 has no annotated start and t5 no annotated stop, so each contributes
    # to one denominator and not the other; the predictor's invented start
    # inside t4 is not charged as a false positive either.
    check("partial 5prime", pr["codon"]["reference_partial_5prime"], 1)
    check("partial 3prime", pr["codon"]["reference_partial_3prime"], 1)
    check("partial start tp", pr["codon"]["start"]["tp"], 2)
    check("partial start fp", pr["codon"]["start"]["fp"], 0)
    check("partial start fn", pr["codon"]["start"]["fn"], 0)
    check("partial stop tp", pr["codon"]["stop"]["tp"], 2)
    check("partial stop fp", pr["codon"]["stop"]["fp"], 0)
    # An identical file must still be exactly 1.0 with the filter on, and the
    # audit flag must put the dropped rows back.
    pid = score(pref, pref, "fixture")
    check("pseudo identity locus f1", pid["locus"]["f1"], 1.0)
    check("pseudo identity tx f1", pid["transcript"]["f1"], 1.0)
    check("pseudo identity start f1", pid["codon"]["start"]["f1"], 1.0)
    pall = score(pref, ppred, "fixture", all_transcripts=True)
    check("all-transcripts ref kept", pall["reference_transcripts"], 5)
    check("all-transcripts locus fn", pall["locus"]["fn"], 2)

    # GENCODE's ``gene_type`` and Ensembl's ``biotype`` must drop the same
    # rows RefSeq's ``gene_biotype`` does.  Before the synonyms were read,
    # these two predicted loci were false positives against a reference that
    # is not allowed to contain their answer.
    spred = os.path.join(d, "spred.gff3")
    open(spred, "w").write(PSEUDO_PRED_SYNONYM)
    sr = score(pref, spred, "fixture")
    check("synonym pred dropped",
          sr["predicted_transcript_selection"]["dropped_by_reason"],
          {"biotype=IG_V_gene": 1, "pseudogene": 1})
    # Declaring them costs the submission nothing: every scored count is the
    # one PSEUDO_PRED gets without those two rows.
    check("synonym locus", sr["locus"], pr["locus"])
    check("synonym transcript", sr["transcript"], pr["transcript"])
    check("synonym nucleotide", sr["nucleotide"], pr["nucleotide"])
    # ``--score-all-transcripts`` puts both sides back, and then the two rows
    # match the reference's own V segment and pseudogene.
    sall = score(pref, spred, "fixture", all_transcripts=True)
    check("synonym audit pred kept", sall["predicted_transcripts"], 5)
    check("synonym audit locus fp", sall["locus"]["fp"], 0)
    # A biotype-dropped row is on a scored sequence, so it must not be
    # reported as if the prediction named a sequence the reference lacks.
    check("synonym not scored", sr["predicted_transcripts_not_scored"], 0)
    os.unlink(spred)
    os.unlink(pref)
    os.unlink(ppred)

    # A donor shared by two introns of different lengths must land in a fixed
    # decile.  ``introns`` is a set, so before ``sites_of`` took the minimum
    # the length attached to a shared site depended on the hash seed and the
    # §4.3 stratification of a real reference was not reproducible: on human
    # 12.35% of donors are shared this way, with over 1 Mb between the
    # shortest and the longest intron at one site.  Here g1's two isoforms
    # share the donor at 1099 with a 100 bp and a 6,900 bp intron; the site is
    # counted once, in the decile of the shorter, and reported as ambiguous.
    sdref = os.path.join(d, "sdref.gff3")
    open(sdref, "w").write(SHARED_DONOR_REF)
    sd = score(sdref, sdref, "fixture")
    check("shared donor sites", sd["splice"]["donor"]["tp"], 1)
    check("shared acceptor sites", sd["splice"]["acceptor"]["tp"], 2)
    check("shared donor ambiguous",
          sd["splice"]["donor"]["reference_sites_multiple_intron_lengths"], 1)
    check("shared acceptor ambiguous",
          sd["splice"]["acceptor"]["reference_sites_multiple_intron_lengths"], 0)
    check("shared donor decile cuts", sd["splice"]["intron_length_decile_cuts"],
          [100, 100, 100, 100, 100, 6900, 6900, 6900, 6900])
    # The shorter intron's decile, not the longer one's: bucket 0.
    check("shared donor decile", sd["splice"]["donor"]["by_intron_length_decile"][0]["tp"], 1)
    check("shared donor top decile",
          sd["splice"]["donor"]["by_intron_length_decile"][9]["tp"], 0)
    os.unlink(sdref)

    # A stop codon split across an intron: merging the feature into the chain
    # creates a junction the pre-merge window plan never asked for.  Its
    # dinucleotide has to be read anyway, so the plan is extended after the
    # merge; before that it came back `unknown`, which is how one predicted
    # intron of the S. pombe cross-parameter run was classified.
    ssref = os.path.join(d, "ssref.gff3")
    sspred = os.path.join(d, "sspred.gff3")
    open(ssref, "w").write(SPLIT_STOP_REF)
    open(sspred, "w").write(SPLIT_STOP_PRED)
    ssfa = os.path.join(d, "ss.fa")
    sseq = list("C" * 20000)
    sseq[2101:2103] = list("GT")       # donor of the junction the merge makes
    sseq[2497:2499] = list("AG")       # its acceptor
    sseq = "".join(sseq)
    with open(ssfa, "w") as fh:
        fh.write(">chr1 split-stop fixture\n")
        for i in range(0, len(sseq), 60):
            fh.write(sseq[i:i + 60] + "\n")
    ss = score(ssref, sspred, "fixture", stop_in_cds=False, genome=ssfa)
    ssd = ss["splice"]["by_dinucleotide"]
    check("split stop merged",
          ss["codon"]["transcripts_stop_merged_from_feature"], 1)
    check("split stop predicted introns", ss["splice"]["predicted_introns"], 2)
    check("split stop no unknown class", "unknown" in ssd, False)
    check("split stop junction classified", ssd.get("GT-AG", {}).get("fp"), 1)
    check("split stop reference intron matched", ssd.get("other", {}).get("tp"), 1)
    os.unlink(ssref)
    os.unlink(sspred)
    os.unlink(ssfa)

    # An annotated isoform that contains the prediction's whole CDS shares
    # exactly as many bases with it as the isoform the prediction equals, so
    # ordering the within-locus pairing on shared bases alone leaves the
    # winner to the id tie-break, and `rna-t1` sorts before `rna-t9`.  That
    # scored 343 exact fugu matches as misses.  Exactness outranks overlap.
    ciref = os.path.join(d, "ciref.gff3")
    cipred = os.path.join(d, "cipred.gff3")
    open(ciref, "w").write(CONTAINED_ISOFORM_REF)
    open(cipred, "w").write(CONTAINED_ISOFORM_PRED)
    ci = score(ciref, cipred, "fixture")
    check("contained isoform tx tp", ci["transcript"]["tp"], 1)
    check("contained isoform tx fp", ci["transcript"]["fp"], 0)
    check("contained isoform tx fn", ci["transcript"]["fn"], 0)
    check("contained isoform tx f1", ci["transcript"]["f1"], 1.0)
    check("contained isoform locus tp", ci["locus"]["tp"], 1)
    # The identity run must stay 1.0: every prediction now has an exact
    # partner, and the contained isoform must not steal the other's.
    ciid = score(ciref, ciref, "fixture")
    check("contained isoform identity tx f1", ciid["transcript"]["f1"], 1.0)
    check("contained isoform identity tx tp", ciid["transcript"]["tp"], 2)
    os.unlink(ciref)
    os.unlink(cipred)

    select_test(check)

    # The window fetcher must return the same bases as a plain read.
    fa = os.path.join(d, "g.fa")
    seq = "ACGT" * 500
    with open(fa, "w") as fh:
        fh.write(">chr1 test\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    wf = WindowFetcher(fa, [("chr1", 1, 4), ("chr1", 100, 103), ("chr1", 1997, 2000)])
    check("window 1-4", wf("chr1", 1, 4), seq[0:4])
    check("window 100-103", wf("chr1", 100, 103), seq[99:103])
    check("window 1997-2000", wf("chr1", 1997, 2000), seq[1996:2000])

    # A window running past the end of a record must be clipped, and must not
    # block the windows queued behind it -- ``pending`` is ordered by window
    # start, so an unsatisfiable window with an early start hides every later
    # one on the same sequence.  A splice site within GC_WINDOW // 2 of a
    # scaffold end produces exactly this, on chr1 and again on the last record
    # in the file.
    fa2 = os.path.join(d, "g2.fa")
    s1, s2 = "ACGT" * 25, "GGCC" * 25  # 100 bp each
    with open(fa2, "w") as fh:
        for name, s in (("chr1", s1), ("chr2", s2)):
            fh.write(">%s test\n" % name)
            for i in range(0, len(s), 30):
                fh.write(s[i:i + 30] + "\n")
    wf2 = WindowFetcher(fa2, [("chr1", 90, 150), ("chr1", 95, 96),
                              ("chr2", 90, 150), ("chr2", 95, 96),
                              ("chr3", 1, 4)])
    check("overrun clipped chr1", wf2("chr1", 90, 150), s1[89:100])
    check("behind overrun chr1", wf2("chr1", 95, 96), s1[94:96])
    check("overrun clipped last record", wf2("chr2", 90, 150), s2[89:100])
    check("behind overrun last record", wf2("chr2", 95, 96), s2[94:96])
    check("absent sequence", wf2("chr3", 1, 4), None)

    for f in (ref, pred, fa, fa2):
        os.unlink(f)
    os.rmdir(d)
    if fails:
        for f in fails:
            print("FAIL " + f, file=sys.stderr)
        return 1
    print("self-test: %d checks passed" % n_checks[0])
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", help="reference GFF3 (.gz ok)")
    ap.add_argument("--prediction", help="predicted GFF3 (.gz ok)")
    ap.add_argument("--species", help="panel species name, e.g. Homo_sapiens")
    ap.add_argument("--declaration", help="submission declaration, docs/benchmark.md 3.3")
    ap.add_argument("--genome", help="reference FASTA; enables the splice "
                                     "dinucleotide and local-GC strata")
    ap.add_argument("--seqids", help="file of sequence ids to score, one per "
                                     "line; default is every reference sequence "
                                     "over 10 kb that is not an alt locus or patch")
    ap.add_argument("--stop-outside-cds", action="store_true",
                    help="the prediction excludes the stop codon from its CDS")
    ap.add_argument("--genetic-code", type=int, default=1,
                    choices=sorted(STOP_CODONS),
                    help="NCBI translation table of the species, used only to "
                         "read the stop-codon convention off --genome "
                         "(panel.tsv column genetic_code; 6 for Tetrahymena)")
    ap.add_argument("--score-all-transcripts", action="store_true",
                    help="score pseudogene and gene-fragment CDS rows as "
                         "protein-coding truth (section 4; for auditing only)")
    ap.add_argument("--out", help="write JSON here instead of stdout")
    ap.add_argument("--self-test", action="store_true",
                    help="score built-in fixtures and check the answers")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    for req in ("reference", "prediction", "species"):
        if not getattr(args, req):
            ap.error("--%s is required" % req)
    if not args.declaration:
        ap.error("--declaration is required: a run without the section 3.3 "
                 "block is not comparable and is not scored (--self-test is "
                 "the only exception)")
    decl, missing, decl_sha = parse_declaration(args.declaration)
    if missing:
        print("declaration %s is missing required keys: %s"
              % (args.declaration, ", ".join(missing)), file=sys.stderr)
        return 2

    seqids = None
    if args.seqids:
        seqids = [l.strip() for l in open(args.seqids) if l.strip()]

    row = score(args.reference, args.prediction, args.species,
                genome=args.genome, seqids=seqids,
                stop_in_cds=not args.stop_outside_cds,
                all_transcripts=args.score_all_transcripts,
                genetic_code=args.genetic_code)
    row["declaration"] = {"path": os.path.basename(args.declaration),
                          "sha256": decl_sha,
                          "model": decl.get("model")}
    # A prediction in a different sequence naming convention scores 0.0
    # everywhere and looks exactly like a bad predictor.  Ensembl calls
    # C. elegans chromosome I "I" and RefSeq calls it NC_003279.8, so this is
    # one wrong download away at all times.  Say so on stderr; the counts are
    # in the JSON either way.
    # The stop-codon convention is the other silent way a correct prediction
    # scores like a broken one.  AUGUSTUS, BRAKER, GeneMark and SNAP all
    # exclude the stop codon from the CDS; the panel references all include
    # it.  Without --stop-outside-cds an AUGUSTUS run scores 0.00 exon F1 and
    # 0.00 stop-codon F1 on S. cerevisiae while nucleotide F1 stays at 0.96,
    # which reads as a real result.
    detected = row["codon"]["stop_codon_convention_detected"]
    if detected == "outside" and not args.stop_outside_cds:
        print("warning: the prediction's stop_codon features lie OUTSIDE its "
              "CDS (%d of %d transcripts with a stop_codon feature), but the "
              "run was not given --stop-outside-cds. Every terminal exon, "
              "single-exon gene, stop codon and exact transcript match is "
              "being scored 3 bp short. Re-run with --stop-outside-cds."
              % (row["codon"]["stop_outside_cds"],
                 row["codon"]["transcripts_with_stop_codon_feature"]),
              file=sys.stderr)
    elif detected == "inside" and args.stop_outside_cds:
        print("warning: --stop-outside-cds was given but the prediction's "
              "stop_codon features lie INSIDE its CDS (%d of %d). The CDS "
              "chains were left alone where a stop_codon feature said so; "
              "check the flag."
              % (row["codon"]["stop_inside_cds"],
                 row["codon"]["transcripts_with_stop_codon_feature"]),
              file=sys.stderr)
    elif detected == "mixed":
        print("warning: the prediction mixes stop-codon conventions (%d "
              "inside, %d outside of %d transcripts with a stop_codon "
              "feature)."
              % (row["codon"]["stop_inside_cds"],
                 row["codon"]["stop_outside_cds"],
                 row["codon"]["transcripts_with_stop_codon_feature"]),
              file=sys.stderr)
    elif detected == "unknown" and args.stop_outside_cds:
        print("warning: --stop-outside-cds was given but the prediction has "
              "no stop_codon features, so %d CDS chains were extended by a "
              "blind 3 bp. A 3'-partial gene has no stop codon and is "
              "extended wrongly by this."
              % row["codon"]["transcripts_stop_extended_by_3bp"],
              file=sys.stderr)

    # The genome settles what the features cannot.  Helixer and Tiberius emit
    # no stop_codon feature at all, so before this check their convention was
    # the default of a flag; --genome makes it a measurement.
    from_genome = row["codon"]["stop_codon_convention_from_genome"]
    if row["codon"]["stop_codon_convention_source"] == "genome" \
            and row["codon"]["transcripts_stop_extended_by_3bp"]:
        print("note: the prediction has no stop_codon features, but %d of %d "
              "chains carry a stop codon in the 3 bp *after* the CDS and only "
              "%d carry one inside it, so the CDS chains were extended by 3 bp "
              "as if --stop-outside-cds had been given. Pass "
              "--stop-outside-cds to make that explicit."
              % (row["codon"]["stop_convention_genome_outside"],
                 row["codon"]["stop_convention_genome_chains"],
                 row["codon"]["stop_convention_genome_inside"]),
              file=sys.stderr)
    elif from_genome == "inside" and args.stop_outside_cds:
        print("warning: --stop-outside-cds was given, but %d of %d predicted "
              "chains already end in a stop codon on the genome. The blind "
              "3 bp extension has moved the 3' end of every chain past its "
              "own stop; drop the flag."
              % (row["codon"]["stop_convention_genome_inside"],
                 row["codon"]["stop_convention_genome_chains"]),
              file=sys.stderr)
    elif from_genome == "unknown" and detected == "unknown" and args.genome:
        print("warning: the prediction has no stop_codon features and the "
              "genome does not settle the convention either (%d chains: %d "
              "end in a stop, %d are followed by one, %d neither). Terminal "
              "exon, single-exon, stop-codon and exact-transcript scores are "
              "resting on the GFF3 default."
              % (row["codon"]["stop_convention_genome_chains"],
                 row["codon"]["stop_convention_genome_inside"],
                 row["codon"]["stop_convention_genome_outside"],
                 row["codon"]["stop_convention_genome_neither"]),
              file=sys.stderr)
    elif row["codon"]["stop_codon_convention_source"] == "assumed":
        print("warning: nothing in this run settles the stop-codon "
              "convention: the prediction has no stop_codon features and no "
              "--genome was given. Pass --genome so the convention is "
              "measured rather than assumed.", file=sys.stderr)

    if row["predicted_conflicting_transcript_ids"]:
        print("warning: %d ids in the prediction appear on more than one "
              "sequence or strand, so the file is not valid GFF3. They were "
              "disambiguated by (sequence, strand) before scoring, which is "
              "what concatenating one predictor run per scaffold needs; "
              "AUGUSTUS's --uniqueGeneId=true does it at the source."
              % row["predicted_conflicting_transcript_ids"], file=sys.stderr)

    unscored = row["predicted_transcripts_not_scored"]
    total_pred = (row["predicted_transcripts"] + unscored
                  + row["predicted_transcript_selection"]["dropped"])
    if total_pred and unscored > total_pred / 2:
        print("warning: %d of %d predicted transcripts (%.0f%%) are not on a "
              "scored sequence; %d predicted sequence names are absent from "
              "the reference. Check that the prediction uses the reference's "
              "sequence names."
              % (unscored, total_pred, 100.0 * unscored / total_pred,
                 row["predicted_sequences_absent_from_reference"]),
              file=sys.stderr)

    text = json.dumps(row, indent=1, sort_keys=True)
    if args.out:
        open(args.out, "w").write(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
