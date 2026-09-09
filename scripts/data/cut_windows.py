#!/usr/bin/env python3
"""Cut the output of ``fetch_window.py`` into fixed-length training examples.

This file *is* the window-cutting convention (docs/data-sources.md, section
6.3).  Standard library only; it writes NumPy ``.npz`` archives by hand so
that a fresh clone can run it, and PyTorch or NumPy can read the result.

One fetched locus (``<stem>.maf``, ``.fa``, ``.nh``, ``.annotation.json``,
``.conservation.json``, ``.manifest.json``) becomes one or more examples
``<out>/<stem>.w<k>.npz`` plus a sidecar ``<out>/<stem>.w<k>.json``.  Every
array is on the reference + strand, 0-based, one entry per reference base;
insertions relative to the reference are summarised, not expanded, so the
example has a fixed length ``L`` whatever the alignment did.

Arrays in each ``.npz`` (``L`` = ``--length``, ``K`` = number of informants):

  ``ref``    uint8 [L]     A=0 C=1 G=2 T=3, N or other=4, padding=4
  ``mask``   uint8 [L]     1 where the position is a real reference base, 0 padding
  ``label``  uint8 [L]     0 intergenic, 1 intron, 2 UTR exon, 3 CDS
                           (CDS > UTR > intron over the transcripts that
                           paint: all isoforms by default, see --isoforms)
  ``frame``  int8  [L]     codon position 0,1,2 of a CDS base counted from the
                           start codon on the transcript's strand; -1 elsewhere
  ``strand`` int8  [L]     +1 / -1 strand of the transcript that owns the
                           base, 0 for intergenic.  The owner is the transcript
                           giving the base its highest label class; ties go
                           to the shortest span (a gene nested in another
                           gene's intron owns its own exons and introns), then
                           annotation order.  ``frame`` is the owner's too.
  ``bound``  uint8 [L]     0 none, 1 first base of start codon, 2 last base of
                           stop codon, 3 donor (first intron base),
                           4 acceptor (last intron base); all placed at their
                           + strand coordinate, strand given by ``strand``.
                           A CDS end the annotation leaves incomplete gets no
                           start or stop mark (see --partial-ends)
  ``inf``    uint8 [K, L]  informant base opposite each reference base:
                           A=0 C=1 G=2 T=3, gap=4, unaligned=5 (no block covers
                           the position), N or other=6
  ``ins``    uint8 [K, L]  length of the informant insertion immediately after
                           this reference base, clipped at 255
  ``dist``   float32 [K]   patristic distance reference -> informant on the
                           track tree (substitutions per site); NaN if the
                           informant is not on the tree
  ``cons``   float32 [L]   per-base conservation score (phyloP / phastCons)
                           when the fetch had one, else NaN

The sidecar JSON carries the coordinates, the informant names in row order,
how that order was chosen (``row_order``), which rows are inferred ancestors
and which leaves each ancestor summarises, which species and ancestors were
dropped for leakage (``--drop-species``), the alignment class the leakage
rule was applied under, the transcripts painted and the transcripts
excluded by type, the distinct start codons, stop codons, donors and
acceptors the annotation states over all isoforms against those the
painted isoforms keep (``distinct_sites``, so the cost of an isoform
policy in splice-site diversity is on record), which CDS ends were judged
incomplete and on what evidence (``cds_ends``), how alignment blocks were
chosen where the MAF said two things about one base (``block_selection``:
overlapping blocks and what first-wins lost, minus-strand reference rows
flipped, duplicate rows for one species and the policy that picked one,
see --duplicate-rows), the source manifest's SHA-256 and the tool version.

Conventions worth stating once:

* Examples are cut on the fetched window with length ``L`` and stride ``S``
  (defaults 4096 / 2048); the last example is padded, never dropped, and
  ``mask`` says where the padding starts.  ``--length 0`` emits the whole
  fetched window as one unpadded example.
* Incomplete CDS ends.  An annotation states a truncated CDS end either
  explicitly (genePred ``cdsStartStat``/``cdsEndStat`` = ``incmpl``,
  GENCODE ``cds_start_NF``/``cds_end_NF`` tags, a non-zero frame on the
  first coding exon; the fetcher carries these as ``cds_start_status``,
  ``cds_end_status`` and ``cds_start_frame``) or only by omission, and some
  sources (Ensembl REST, UCSC's GENCODE ``cdsStartStat`` columns, which are
  ``none`` throughout) state nothing.  ``--partial-ends both`` (default)
  takes the declaration where there is one and otherwise reads the
  reference: a CDS whose first codon is not ATG is 5' incomplete, one whose
  last codon is not a stop (nor the codon after it, for a stop-excluded
  convention) is 3' incomplete, under ``--genetic-code``.  An incomplete
  end gets no start/stop mark in ``bound``, is not counted in
  ``distinct_sites``, and a declared 5' frame offsets ``frame``.
  ``declared`` and ``sequence`` use one signal only; ``none`` paints every
  end as before version 0.5.  The sidecar's ``cds_ends`` records the
  declared and sequence evidence per window, the transcripts judged
  incomplete at each end, and every disagreement between the two signals.
* Row order.  For a track with one tree (multiz, Cactus on UCSC) rows are
  the tree's leaves in depth-first order, so ``K`` is fixed per track.
  Ensembl returns one tree per block, so rows there are the species set's
  members in sorted order (``species_set_members`` in the fetcher's
  manifest; every member gets a row even when it does not align in the
  window), then any observed species outside the set, then the inferred
  ancestral rows sorted by name.  The leaf part is fixed per species set;
  the ancestor part varies by window, and ``--drop-ancestors`` removes it.
* Reverse-strand genes are *not* flipped here.  ``--both-strands`` writes a
  second example ``.rc.npz`` with every array reverse-complemented (labels
  and boundaries recomputed on the flipped coordinates), so a model can be
  trained orientation-symmetric; the default keeps one orientation and
  leaves augmentation to the loader.  An insertion recorded after the last
  reference base of an example has no preceding base to move to on the
  flipped strand and is dropped from ``.rc``; at most one position per
  informant per example.
* ``--transcript-types`` decides which transcripts paint labels.  The
  default ``benchmark`` follows docs/benchmark.md section 4: pseudogenes and
  immunoglobulin / T-cell receptor segments are excluded, recognised only
  from the type the source states (Ensembl ``biotype``, GenArk ``geneType``;
  a transcript without a type always paints).  ``all`` paints everything;
  a comma-separated list keeps exactly those types.  Excluded transcripts
  are listed in the sidecar.
* Stop codons are inside the CDS, as in RefSeq genePred and GFF3 from NCBI;
  ``bound == 2`` marks the last base of that codon.
* ``--drop-species`` implements the benchmark leakage rule for training
  alignments (docs/benchmark.md 3.2, docs/data-sources.md section 7): a
  held-out species' row is removed before the tensors are built, and so is
  every ancestral row whose clade (read from the block trees) contains a
  dropped species.  For a reference-anchored alignment (multiz, pairwise
  chain/net) this is exact; for a jointly inferred one (Cactus, EPO, PECAN)
  it removes the sequence but not its influence on the alignment, which the
  sidecar records as ``"dropped_rows_only": true`` together with the
  ``alignment_class`` it decided that from.  The class comes from an
  explicit table of track names; a track the table does not know is treated
  as jointly inferred (the strict branch of docs/benchmark.md 3.2) unless
  ``--reference-anchored`` says otherwise, which the sidecar then records as
  the operator's claim.  The sidecar's ``"dropped_species"`` list is the
  value a run declares as ``alignment_rows_dropped`` (docs/benchmark.md
  3.3).
* Rows for species on the tree but absent from every block are kept as
  all-``unaligned``, so ``K`` is fixed per track.

Usage::

    python3 scripts/data/cut_windows.py --stem /tmp/win/HBB/HBB_470 --out /tmp/ex
    python3 scripts/data/cut_windows.py --stem /tmp/win/Adh/Adh_124 --out /tmp/ex \
        --length 2048 --stride 1024 --drop-species apiMel4 --both-strands
    python3 scripts/data/cut_windows.py --self-test
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import shutil
import struct
import sys
import zipfile

TOOL_VERSION = "0.11"

# Length floors the sidecar counts features under (``feature_lengths``).
# They are the ones engels measured in Helixer's decoder (relay note
# 20260909T202745Z-engels-0020): the HMM's shortest intron path is 30
# bases (U12 AT-AC) and the U2 GT-AG / GC-AG paths 50, so a shorter intron
# cannot be decoded at all; a transcript needs a genic run above 80 bases
# to pass the default candidate gate (window 100, peak 0.8) before the
# decoder sees it; and the published minimum coding length is 60.  Any
# decoder with hard minimum durations has floors of this kind, and the
# labels below them are the ones such a model can never reproduce.
LENGTH_FLOORS = {"intron": (30, 50), "cds": (60,), "span": (81,)}
# The scorer's floor: benchmark/score.py ``MIN_INTRON`` treats a gap between
# consecutive CDS blocks as an intron only if it is at least 20 bases and
# counts the shorter ones apart (``reference_cds_gaps_below_min``; lenin,
# relay note 20260909T211910Z-lenin-0019).  The cutter follows the same
# floor: an exon gap shorter than ``--min-intron`` is painted ``short_gap``,
# gets no donor or acceptor mark, and is listed in the sidecar's
# ``short_gaps`` record with its flanking sequence, so what the scorer
# leaves out of its splice denominator is not a splice site in the labels
# either.  Twenty is a policy, not a biological fact: RefSeq encodes a
# programmed ribosomal frameshift as a 1-base gap between two CDS blocks
# (the Ty ORFs on sacCer3 chrIV, ``cdh-4`` on ce11 chrIII), and Stentor
# coeruleus has 15- and 16-base spliceosomal introns (stalin, relay note
# 20260909T214050Z-stalin-0021: 8,806 of them in one annotation), which
# this floor would file under short_gap; ``--min-intron 15`` restores them
# and 0 turns the floor off (the pre-0.9 behaviour).
MIN_INTRON = 20
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
# Stop codons by NCBI genetic code table number, for the sequence check of
# CDS ends.  1 standard; 4 mold/protozoan mitochondrial and Mycoplasma; 6
# ciliate, dasycladacean and hexamita (T. thermophila: TAA and TAG read
# Gln); 10 euplotid (TGA reads Cys); 12 alternative yeast; 25 candidate
# division SR1; 26 Pachysolen; 29-31 Mesodinium, Peritrich, Blastocrithidia.
GENETIC_CODE_STOPS = {1: {"TAA", "TAG", "TGA"}, 4: {"TAA", "TAG"}, 6: {"TGA"}, 10: {"TAA", "TAG"},
                      12: {"TAA", "TAG", "TGA"}, 25: {"TAA", "TAG"}, 26: {"TAA", "TAG", "TGA"},
                      29: {"TGA"}, 30: {"TGA"}, 31: {"TGA"}}
# Near-cognate initiation codons; a CDS opening on one is treated as 5'
# incomplete by the sequence rule (docs/data-sources.md 6.3: on GENCODE
# hg38 chr21, 3,664 of 3,669 untagged CDSs open on ATG, none of the 117
# cds_start_NF ones do, and 15 of those open on a near-cognate codon), and
# the sidecar counts them apart from ATG.
NEAR_COGNATE_STARTS = {"CTG", "GTG", "TTG", "ACG", "ATA", "ATC", "ATT", "AAG", "AGG"}
PARTIAL_END_MODES = ("both", "declared", "sequence", "none")
INF_GAP, INF_UNALIGNED, INF_OTHER = 4, 5, 6
LABEL = {"intergenic": 0, "intron": 1, "utr": 2, "cds": 3, "short_gap": 4}
# precedence when transcripts overlap: a base takes its highest-ranked
# class; short_gap (an exon gap under the intron floor) ranks above intron
# and below UTR, so an isoform that reads through the gap owns it
LABEL_RANK = {"intergenic": 0, "intron": 1, "short_gap": 2, "utr": 3, "cds": 4}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Alignment classes for the benchmark leakage rule (docs/benchmark.md 3.2).
# Reference-anchored: rows are independent given the reference, so dropping a
# row at cut time is a filter.  Jointly inferred: columns and ancestors were
# estimated over the whole taxon set, so a dropped row still shaped the
# alignment and the rebuild rule applies.  A track in neither table is
# "unknown" and handled as jointly inferred.
REFERENCE_ANCHORED_TRACKS = re.compile(r"^(multiz\d+way\w*|(chain|net|lastz)\w+)$", re.I)
JOINTLY_INFERRED_TRACKS = re.compile(r"^(cactus\d+way\w*|hprc\d+way\w*)$", re.I)
# Ensembl EPO ancestral rows are named after their descendants: Ggal-Mgal[2].
ANCESTOR_RE = re.compile(r"^[A-Za-z]+(-[A-Za-z]+)*\[\d+\]$")
# Transcript types the benchmark removes from truth (docs/benchmark.md 4):
# pseudogenes (any biotype containing the word) and Ig/TCR segments in
# Ensembl (IG_V_gene ...) or NCBI (V_segment, C_region) spelling.
BENCHMARK_EXCLUDED_TYPES = re.compile(r"pseudogene|^(IG|TR)_[VDJC]_gene$|^[VDJ]_segment$|^C_region$", re.I)


def classify_alignment(manifest: dict, reference_anchored: bool = False) -> tuple[str, str]:
    """Return (class, how): class is reference_anchored, jointly_inferred or
    unknown; how says whether the table or the operator decided."""
    if reference_anchored:
        return "reference_anchored", "--reference-anchored"
    if manifest.get("source") == "ensembl":
        if str(manifest.get("method") or "").upper().startswith("LASTZ"):
            return "reference_anchored", "table"
        return "jointly_inferred", "table"
    track = str(manifest.get("alignment_track") or "")
    if REFERENCE_ANCHORED_TRACKS.match(track):
        return "reference_anchored", "table"
    if JOINTLY_INFERRED_TRACKS.match(track):
        return "jointly_inferred", "table"
    return "unknown", "default-strict"


def select_transcripts(transcripts: list[dict], mode: str) -> tuple[list[dict], list[dict]]:
    """Split transcripts into (painted, excluded) by ``--transcript-types``.
    A transcript without a stated type is always painted, the same
    "recognised only from what the file states" rule the benchmark uses."""
    keep, out = [], []
    allowed = None if mode in ("benchmark", "all") else {t.strip() for t in mode.split(",") if t.strip()}
    for t in transcripts:
        ty = str(t.get("type") or "")
        if mode == "all" or not ty:
            keep.append(t)
        elif mode == "benchmark" and BENCHMARK_EXCLUDED_TYPES.search(ty):
            out.append({"id": t.get("id"), "type": ty})
        elif allowed is not None and ty not in allowed:
            out.append({"id": t.get("id"), "type": ty})
        else:
            keep.append(t)
    return keep, out


def _strip_version(tid: str) -> str:
    return re.sub(r"\.\d+$", "", str(tid))


def group_loci(transcripts: list[dict]) -> list[list[int]]:
    """Group transcript indices into loci.  Transcripts that state a gene
    share a locus by that name; the rest are clustered by overlapping span
    on the same strand (single linkage), which is the coverage tables' rule
    for an annotation that carries no gene names.  Loci keep annotation
    order of their first member."""
    named: dict[str, list[int]] = {}
    unnamed: list[int] = []
    order: list[tuple[int, str]] = []
    for i, t in enumerate(transcripts):
        g = t.get("gene")
        if g:
            key = "gene:" + str(g)
            if key not in named:
                named[key] = []
                order.append((i, key))
            named[key].append(i)
        else:
            unnamed.append(i)
    # single-linkage clustering of the unnamed transcripts by strand and span
    clusters: list[dict] = []
    for i in sorted(unnamed, key=lambda j: (transcripts[j].get("strand", "+"), transcripts[j]["start"])):
        t = transcripts[i]
        st = t.get("strand", "+")
        if clusters and clusters[-1]["strand"] == st and t["start"] < clusters[-1]["end"]:
            clusters[-1]["members"].append(i)
            clusters[-1]["end"] = max(clusters[-1]["end"], t["end"])
        else:
            clusters.append({"strand": st, "start": t["start"], "end": t["end"], "members": [i]})
    groups = {key: idx for key, idx in named.items()}
    for c in clusters:
        key = f"span:{c['strand']}{c['start']}-{c['end']}"
        groups[key] = sorted(c["members"])
        order.append((min(c["members"]), key))
    order.sort()
    return [groups[key] for _, key in order]


def locus_name(transcripts: list[dict], members: list[int]) -> str:
    g = transcripts[members[0]].get("gene")
    if g:
        return str(g)
    lo = min(transcripts[i]["start"] for i in members)
    hi = max(transcripts[i]["end"] for i in members)
    return f"{transcripts[members[0]].get('strand', '+')}{lo}-{hi}"


def select_isoforms(transcripts: list[dict], policy: str = "union",
                    representatives: set[str] | None = None) -> tuple[list[dict], list[dict], list[str]]:
    """Choose which isoforms of each locus paint labels.

    ``union``: every transcript paints (the label is a union over isoforms).
    ``longest-cds``: one transcript per locus, the one with the most CDS
    bases; ties go to the longest exonic span, then annotation order; a
    locus with no coding transcript keeps its longest exonic span.
    ``representative``: only the transcripts named in ``representatives``
    paint (ids compared with and without a trailing version); a locus in
    which none is named falls back to ``longest-cds`` and is reported.

    Returns (kept, dropped, fallback_loci) with dropped entries of the form
    {"id", "locus", "reason"}.  The benchmark scores any of these targets
    (relay note 20260909T141456Z-lenin-0014): the choice is a training
    decision, and the sidecar records which one was made."""
    if policy == "union":
        return list(transcripts), [], []
    if policy not in ("longest-cds", "representative"):
        raise SystemExit(f"unknown --isoforms policy {policy!r}")
    reps: set[str] = set()
    for r in representatives or ():
        reps.add(str(r)); reps.add(_strip_version(r))
    keep_idx: set[int] = set()
    dropped: list[dict] = []
    fallback: list[str] = []

    def cds_len(t: dict) -> int:
        return sum(b - a for a, b in t.get("cds", []))

    def exon_len(t: dict) -> int:
        return sum(b - a for a, b in t.get("exons", []))

    for members in group_loci(transcripts):
        name = locus_name(transcripts, members)
        chosen: list[int] = []
        reason = "shorter-cds"
        if policy == "representative":
            chosen = [i for i in members
                      if str(transcripts[i].get("id")) in reps or _strip_version(transcripts[i].get("id")) in reps]
            reason = "not-representative"
            if not chosen:
                fallback.append(name)
        if not chosen:
            best = max(members, key=lambda i: (cds_len(transcripts[i]), exon_len(transcripts[i]), -i))
            chosen = [best]
            if policy == "representative":
                reason = "shorter-cds"
        keep_idx.update(chosen)
        for i in members:
            if i not in keep_idx:
                dropped.append({"id": transcripts[i].get("id"), "locus": name, "reason": reason})
    kept = [t for i, t in enumerate(transcripts) if i in keep_idx]
    return kept, dropped, fallback


def read_id_list(path: str) -> set[str]:
    """One transcript id per line; blank lines and ``#`` comments ignored;
    a tab-separated file contributes its first column."""
    ids: set[str] = set()
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            ids.add(ln.split("\t")[0].split()[0])
    return ids


# --------------------------------------------------------------- npy/npz

def npy_bytes(data: bytes, dtype: str, shape: tuple[int, ...]) -> bytes:
    header = "{'descr': '%s', 'fortran_order': False, 'shape': %s, }" % (
        dtype, "(%s)" % "".join(f"{n}, " for n in shape) if shape else "()")
    pad = 64 - ((10 + len(header) + 1) % 64)
    header = header + " " * pad + "\n"
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header.encode("latin1") + data


def write_npz(path: str, arrays: dict[str, tuple[bytes, str, tuple[int, ...]]]) -> None:
    """Zip members carry a fixed timestamp, so the same inputs give the same
    bytes and a checksum of an example means something."""
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name, (data, dtype, shape) in arrays.items():
            info = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, npy_bytes(data, dtype, shape))


def f32(values) -> bytes:
    return struct.pack("<%df" % len(values), *values)


# ------------------------------------------------------------------ Newick

def parse_newick(text: str):
    """Return parallel lists (children, parent, length, name) over the nodes."""
    text = text.strip().rstrip(";")
    children: list[list[int]] = [[]]
    parent = [-1]
    length = [0.0]
    name = [""]

    def new(p: int) -> int:
        children.append([]); parent.append(p); length.append(0.0); name.append("")
        children[p].append(len(children) - 1)
        return len(children) - 1

    def flush(node: int, tok: str) -> None:
        tok = tok.strip()
        if not tok:
            return
        nm, _, ln = tok.partition(":")
        if ln:
            try:
                length[node] = float(ln)
            except ValueError:
                pass
        nm = nm.strip().strip("'")
        if nm and not name[node]:
            name[node] = nm

    cur, buf = 0, ""
    for ch in text:
        if ch == "(":
            cur = new(cur)
        elif ch == ",":
            flush(cur, buf); buf = ""; cur = new(parent[cur])
        elif ch == ")":
            flush(cur, buf); buf = ""; cur = parent[cur]
        elif not ch.isspace():
            buf += ch
    flush(cur, buf)
    return children, parent, length, name



def maf_source(src: str) -> str:
    """Assembly part of a MAF ``s`` row name, spelled as the track's tree
    spells its leaves.  UCSC names rows ``<db>.<chrom>``; GenArk assemblies
    carry a version dot in the db (``GCF_003668045.3.NC_048596.1`` in the
    mm39 35-way) and their tree leaf is ``GCF_003668045v3``."""
    m = re.match(r"^(GC[AF]_\d+)\.(\d+)\.", src)
    if m:
        return f"{m.group(1)}v{m.group(2)}"
    return src.split(".")[0]


def canonical_name(label: str, names=()) -> str:
    """Map a tree label to a MAF row name.  UCSC leaves are the row names
    themselves.  Ensembl leaves are ``species_region_start_end[strand]``: the
    longest known row name that prefixes the label wins, with a regex as the
    fallback when the row was not seen."""
    if label in names:
        return label
    best = ""
    for n in names:
        if label.startswith(n + "_") and len(n) > len(best):
            best = n
    if best:
        return best
    m = re.match(r"^(.+?)_[^_]+_\d+_\d+\[[+-]\]$", label)
    return m.group(1) if m else maf_source(label)


def newick_leaves_and_distances(text: str, ref: str, names=()) -> tuple[list[str], dict[str, float]]:
    """Depth-first leaf order and patristic distance from ``ref`` to every
    named node, leaves and named ancestors alike."""
    children, parent, length, name = parse_newick(text)
    canon = [canonical_name(n, names) if n else "" for n in name]
    leaves = [canon[i] for i in range(len(name)) if not children[i] and name[i]]
    idx: dict[str, int] = {}
    for i, n in enumerate(canon):
        if n and n not in idx:
            idx[n] = i
    key = ref if ref in idx else next((n for n in idx if n.split(".")[0] == ref), None)
    dist: dict[str, float] = {}
    if key is not None:
        up: dict[int, float] = {}
        node, d = idx[key], 0.0
        while node >= 0:
            up[node] = d; d += length[node]; node = parent[node]

        def down(node: int, d: float) -> None:
            for c in children[node]:
                if c in up:
                    continue
                dc = d + length[c]
                if canon[c]:
                    dist.setdefault(canon[c], dc)
                down(c, dc)

        for node, d in up.items():
            if canon[node] and node != idx[key]:
                dist.setdefault(canon[node], d)
            down(node, d)
    return leaves, dist


def newick_clades(text: str, names=()) -> dict[str, set[str]]:
    """Named internal nodes -> the set of leaf names under each.  On Ensembl
    block trees the internal names are the ancestral rows' names."""
    children, parent, length, name = parse_newick(text)
    canon = [canonical_name(n, names) if n else "" for n in name]

    def under(i: int) -> set[str]:
        if not children[i]:
            return {canon[i]} if canon[i] else set()
        out: set[str] = set()
        for c in children[i]:
            out |= under(c)
        return out

    return {canon[i]: under(i) for i in range(len(name)) if children[i] and canon[i]}


# -------------------------------------------------------------------- MAF

def maf_blocks(text: str) -> list[list[list[str]]]:
    blocks, cur = [], []
    for ln in text.splitlines():
        if ln.startswith("a"):
            if cur:
                blocks.append(cur)
            cur = []
        elif ln.startswith("s "):
            cur.append(ln.split())
    if cur:
        blocks.append(cur)
    return blocks


def flip_block(block: list[list[str]]) -> list[list[str]]:
    """Reverse-complement every row of a MAF block whose reference row is on
    the - strand, so the reference reads + with its forward start
    ``srcSize - start - size`` (UCSC MAF: a '-' row's start counts from the
    end of the reverse-complemented source).  Informant coordinates are not
    used downstream, so only their sequence and strand symbol are updated."""
    out = []
    for r in block:
        r = list(r)
        src_size, start, size = int(r[5]), int(r[2]), int(r[3])
        r[2] = str(src_size - start - size)
        r[4] = "-" if r[4] == "+" else "+"
        r[6] = r[6].translate(COMP)[::-1]
        out.append(r)
    return out


def paint_informants(blocks, ref: str, start: int, end: int, rows: list[str],
                     duplicate_rows: str = "identity"):
    """Return inf [K][n], ins [K][n] bytearrays and a dict describing block
    selection: overlapping blocks (the first block in file order to cover a
    reference position wins; Cactus can overlap on the reference), the window
    positions covered more than once, the aligned informant bases the losing
    blocks carried at those positions, how many of them the winning block did
    not have (unaligned or gap, so first-wins loses that base), the same
    count per informant, the blocks whose reference row was on the -
    strand and were flipped to + before painting, and the rows discarded
    because one species had several rows in one block (Cactus exports every
    copy of a duplicated region): ``duplicate_rows`` ``identity`` keeps the
    copy with the most bases identical to the reference in that block, ties
    to the first; ``first`` keeps the first row in block order."""
    n = end - start
    K = len(rows)
    row_of = {r: i for i, r in enumerate(rows)}
    inf = [bytearray([INF_UNALIGNED]) * n for _ in range(K)]
    ins = [bytearray(n) for _ in range(K)]
    seen = bytearray(n)
    twice = bytearray(n)
    overlaps = 0
    minus = 0
    discarded = 0
    lost = 0
    lost_by: dict[str, int] = {}
    dup_blocks = 0
    dup_rows = 0
    dup_bases = 0
    dup_max = 0
    # per informant, each in its own unit: MAF rows discarded, blocks in which
    # the informant had several rows, the most rows it had in one block, and
    # aligned bases inside the window that the discarded rows carried
    dup_by: dict[str, dict[str, int]] = {}
    ref_dups = 0
    for block in blocks:
        rref = next((r for r in block if maf_source(r[1]) == ref), None)
        if rref is None:
            continue
        if rref[4] != "+":
            minus += 1
            block = flip_block(block)
            rref = next(r for r in block if maf_source(r[1]) == ref)
        rseq, p = rref[6], int(rref[2])
        col_ref = []
        for ch in rseq:
            if ch == "-":
                col_ref.append(-1)
            else:
                col_ref.append(p); p += 1
        ref_dups += sum(1 for r in block if maf_source(r[1]) == ref) - 1
        # first block to cover a position wins (Cactus can overlap on the reference)
        cover = [q for q in col_ref if start <= q < end]
        if cover and any(seen[q - start] for q in cover):
            overlaps += 1
        # one row per informant per block
        by_src: dict[str, list[list[str]]] = {}
        for r in block:
            if r is rref:
                continue
            src = maf_source(r[1])
            if src in row_of:
                by_src.setdefault(src, []).append(r)
        chosen: list[list[str]] = []
        block_has_dup = False
        for src, cands in by_src.items():
            if len(cands) > 1:
                block_has_dup = True
                if duplicate_rows == "identity":
                    def ident(r):
                        return sum(1 for a, b in zip(r[6].upper(), rseq.upper()) if a == b and a != "-")
                    best = max(range(len(cands)), key=lambda i: (ident(cands[i]), -i))
                else:
                    best = 0
                d = dup_by.setdefault(src, {"rows": 0, "blocks": 0, "max_copies": 0, "bases": 0})
                d["blocks"] += 1
                d["max_copies"] = max(d["max_copies"], len(cands))
                dup_max = max(dup_max, len(cands))
                for i, r in enumerate(cands):
                    if i != best:
                        dup_rows += 1
                        d["rows"] += 1
                        b = sum(1 for col, ch in enumerate(r[6])
                                if ch != "-" and start <= col_ref[col] < end)
                        d["bases"] += b
                        dup_bases += b
                chosen.append(cands[best])
            else:
                chosen.append(cands[0])
        if block_has_dup:
            dup_blocks += 1
        for r in chosen:
            k = row_of[maf_source(r[1])]
            seq = r[6]
            last = -1
            for col, q in enumerate(col_ref):
                if q < 0:
                    if last >= 0 and seq[col] != "-" and not seen[last - start]:
                        if ins[k][last - start] < 255:
                            ins[k][last - start] += 1
                    continue
                if q < start or q >= end:
                    last = q if start <= q < end else -1
                    continue
                last = q
                ch = seq[col].upper()
                code = BASE.get(ch, INF_GAP if ch == "-" else INF_OTHER)
                if seen[q - start]:
                    if code != INF_GAP:
                        discarded += 1
                        if inf[k][q - start] in (INF_UNALIGNED, INF_GAP):
                            lost += 1
                            lost_by[rows[k]] = lost_by.get(rows[k], 0) + 1
                    continue
                inf[k][q - start] = code
        for q in cover:
            if seen[q - start]:
                twice[q - start] = 1
            seen[q - start] = 1
    kept = sum(1 for k in range(K) for v in inf[k] if v not in (INF_UNALIGNED, INF_GAP))
    stats = {
        "units": {
            "overlapping_blocks": "MAF blocks",
            "overlap_positions": "reference positions in the window",
            "overlap_informant_bases_discarded": "aligned informant bases in the window (one per row and position)",
            "overlap_informant_bases_lost": "aligned informant bases in the window",
            "overlap_lost_by_informant": "aligned informant bases in the window",
            "reference_minus_strand_blocks_flipped": "MAF blocks",
            "reference_duplicate_rows": "MAF rows",
            "duplicate_row_policy": "policy name (identity or first), not a count",
            "duplicate_row_blocks": "MAF blocks with several rows for at least one informant",
            "duplicate_informants": "informants with at least one discarded row",
            "duplicate_max_copies": "rows of one informant in one block",
            "duplicate_rows_discarded": "MAF rows",
            "duplicate_rows_discarded_bases_in_window": "aligned informant bases in the window",
            "kept_informant_bases_in_window": "aligned informant bases in the window (cells of the example that are a base)",
            "duplicates_by_informant": "one record per informant with a discarded row; fields below",
            "duplicates_by_informant.rows": "MAF rows discarded",
            "duplicates_by_informant.blocks": "MAF blocks in which the informant had several rows",
            "duplicates_by_informant.max_copies": "rows of the informant in one block",
            "duplicates_by_informant.bases": "aligned bases in the window carried by the discarded rows",
        },
        "overlapping_blocks": overlaps,
        "overlap_positions": sum(twice),
        "overlap_informant_bases_discarded": discarded,
        "overlap_informant_bases_lost": lost,
        "overlap_lost_by_informant": dict(sorted(lost_by.items(), key=lambda kv: (-kv[1], kv[0]))),
        "reference_minus_strand_blocks_flipped": minus,
        "reference_duplicate_rows": ref_dups,
        "duplicate_row_policy": duplicate_rows,
        "duplicate_row_blocks": dup_blocks,
        "duplicate_informants": len(dup_by),
        "duplicate_max_copies": dup_max,
        "duplicate_rows_discarded": dup_rows,
        "duplicate_rows_discarded_bases_in_window": dup_bases,
        "kept_informant_bases_in_window": kept,
        "duplicates_by_informant": dict(sorted(dup_by.items(), key=lambda kv: (-kv[1]["rows"], kv[0]))),
    }
    return inf, ins, stats


# ------------------------------------------------------------- annotation

def _declared_status(t: dict, which: str) -> str | None:
    """'complete', 'incomplete' or None from what the fetcher recorded."""
    v = t.get(f"cds_{which}_status")
    if v in ("complete", "incomplete"):
        return v
    if v in ("cmpl", "incmpl"):
        return "complete" if v == "cmpl" else "incomplete"
    if which == "start":
        fr = t.get("cds_start_frame")
        if isinstance(fr, int) and fr in (1, 2):
            return "incomplete"
    return None


def cds_end_status(transcripts: list[dict], seq: str, start: int, end: int,
                   stops: set[str], mode: str) -> dict:
    """Decide for every coding transcript whether its CDS is complete at the
    5' and the 3' end, set ``_partial5`` / ``_partial3`` / ``_frame0`` on
    the transcript dicts for the painter, and return the ``cds_ends``
    sidecar record.  ``seq`` is the window's reference sequence (+ strand,
    starting at ``start``).  The declaration wins where the annotation
    makes one; where it does not, the sequence decides (mode ``both``);
    ``declared`` ignores the sequence, ``sequence`` ignores the
    declaration, ``none`` marks nothing as incomplete."""
    if mode not in PARTIAL_END_MODES:
        raise ValueError(f"--partial-ends must be one of {PARTIAL_END_MODES}")
    declared = {k: 0 for k in ("start_complete", "start_incomplete", "start_undeclared",
                               "stop_complete", "stop_incomplete", "stop_undeclared")}
    sequence = {k: 0 for k in ("start_atg", "start_near_cognate", "start_other", "start_outside_window",
                               "stop_in_cds", "stop_after_cds", "stop_none", "stop_outside_window",
                               "cds_length_not_multiple_of_3")}
    partial5, partial3, disagreements, frames = [], [], [], {}

    def window(a: int, b: int) -> str | None:
        if a < start or b > end:
            return None
        return seq[a - start:b - start].upper()

    for t in transcripts:
        cds = sorted(t.get("cds", []))
        t["_partial5"] = t["_partial3"] = False
        fr = t.get("cds_start_frame")
        t["_frame0"] = fr if isinstance(fr, int) and fr in (0, 1, 2) else 0
        if not cds:
            continue
        minus = t.get("strand", "+") != "+"
        length = sum(b - a for a, b in cds)
        if length % 3:
            sequence["cds_length_not_multiple_of_3"] += 1
        pieces = [window(a, b) for a, b in cds]
        cds_seq = None if any(p is None for p in pieces) else "".join(pieces)
        if cds_seq is not None and minus:
            cds_seq = cds_seq.translate(COMP)[::-1]
        # the codon after the CDS, for the stop-excluded convention
        after = window(cds[0][0] - 3, cds[0][0]) if minus else window(cds[-1][1], cds[-1][1] + 3)
        if after is not None and minus:
            after = after.translate(COMP)[::-1]
        d5, d3 = _declared_status(t, "start"), _declared_status(t, "end")
        declared["start_" + (d5 or "undeclared")] += 1
        declared["stop_" + (d3 or "undeclared")] += 1
        if cds_seq is None or len(cds_seq) < 3:
            s5 = s3 = "outside_window"
            sequence["start_outside_window"] += 1
            sequence["stop_outside_window"] += 1
        else:
            first, last = cds_seq[:3], cds_seq[-3:]
            s5 = "atg" if first == "ATG" else ("near_cognate" if first in NEAR_COGNATE_STARTS else "other")
            s3 = "in_cds" if last in stops else ("after_cds" if after in stops else "none")
            sequence["start_" + s5] += 1
            sequence["stop_" + s3] += 1
        seq5 = None if s5 == "outside_window" else ("complete" if s5 == "atg" else "incomplete")
        seq3 = None if s3 == "outside_window" else ("complete" if s3 in ("in_cds", "after_cds") else "incomplete")
        for which, dec, sq, obs in (("start", d5, seq5, s5), ("stop", d3, seq3, s3)):
            if dec and sq and dec != sq:
                disagreements.append({"id": t.get("id"), "end": which, "declared": dec,
                                      "sequence": obs, "decided": "declared" if mode != "sequence" else "sequence"})
            if mode == "none":
                verdict = "complete"
            elif mode == "declared":
                verdict = dec or "complete"
            elif mode == "sequence":
                verdict = sq or "complete"
            else:
                verdict = dec or sq or "complete"
            if verdict == "incomplete":
                if which == "start":
                    t["_partial5"] = True
                    partial5.append(t.get("id"))
                else:
                    t["_partial3"] = True
                    partial3.append(t.get("id"))
        if mode == "none":
            t["_frame0"] = 0
        elif t["_frame0"]:
            frames[str(t.get("id"))] = t["_frame0"]
    return {"policy": mode, "declared": declared, "sequence": sequence,
            "partial_5prime": partial5, "partial_3prime": partial3,
            "frame_offsets": frames, "disagreements": disagreements}


def paint_labels(start: int, end: int, transcripts: list[dict], min_intron: int = MIN_INTRON):
    """Per-base label, CDS frame, strand and boundary marks.

    A base belongs to the transcript that gives it the highest label class
    (CDS > UTR > short_gap > intron, ``LABEL_RANK``); an exon gap shorter
    than ``min_intron`` is a short_gap, not an intron, and gets no donor or
    acceptor mark (``MIN_INTRON``).  Among transcripts tied at that class the one with
    the shortest span owns it, then annotation order.  Strand and frame are
    the owner's.  The tie rule is what makes a gene nested in another gene's
    intron (fly Adh inside the minus-strand gene that spans it, and the
    intronic genes common in large vertebrate loci) carry its own strand
    across its exons and introns instead of the enclosing gene's, so that
    frame and boundary marks, which are always counted on the owner's
    strand, agree with the strand channel.  Boundary marks are the union
    over transcripts; a position marked by two transcripts keeps the first."""
    n = end - start
    RANK = {LABEL[k]: LABEL_RANK[k] for k in LABEL}
    label = bytearray(n)
    frame = bytearray(b"\xff" * n)      # -1 as int8
    strand = bytearray(n)                # 0, 1 (+), 255 (-) as int8
    bound = bytearray(n)

    def put(arr, pos, val, force=False):
        if start <= pos < end and (force or arr[pos - start] == 0):
            arr[pos - start] = val

    # longest span first and, within a span, last in annotation order first:
    # a later painter overwrites on ties, so the nested transcript wins, and
    # among equal spans the first in annotation order does
    order = sorted(range(len(transcripts)),
                   key=lambda i: (-(transcripts[i]["end"] - transcripts[i]["start"]), -i))
    for ti in order:
        t = transcripts[ti]
        s = 1 if t.get("strand", "+") == "+" else 255
        lo, hi = max(t["start"], start), min(t["end"], end)
        if lo >= hi:
            continue
        cls = bytearray(n)
        fr = bytearray(b"\xff" * n)
        for i in range(lo, hi):
            cls[i - start] = LABEL["intron"]
        ex = sorted(t.get("exons", []))
        for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
            if 0 < a2 - b1 < min_intron:
                for i in range(max(b1, start), min(a2, end)):
                    cls[i - start] = LABEL["short_gap"]
        for a, b in t.get("exons", []):
            for i in range(max(a, start), min(b, end)):
                cls[i - start] = LABEL["utr"]
        cds = sorted(t.get("cds", []))
        for a, b in cds:
            for i in range(max(a, start), min(b, end)):
                cls[i - start] = LABEL["cds"]
        if cds:
            k = int(t.get("_frame0") or 0)  # codon position of the first CDS base (0 for a complete 5' end)
            cds_order = cds if s == 1 else [(b, a) for a, b in reversed(cds)]
            for a, b in cds_order:
                rng = range(a, b) if s == 1 else range(a - 1, b - 1, -1)
                for i in rng:
                    if start <= i < end:
                        fr[i - start] = k % 3
                    k += 1
        for i in range(lo - start, hi - start):
            c = cls[i]
            if c and RANK[c] >= RANK[label[i]]:
                label[i] = c
                strand[i] = s
                frame[i] = fr[i]
        if cds:
            # an incomplete end (cds_end_status) has no codon to mark
            mark_start, mark_stop = not t.get("_partial5"), not t.get("_partial3")
            if s == 1:
                if mark_start:
                    put(bound, cds[0][0], 1)
                if mark_stop:
                    put(bound, cds[-1][1] - 1, 2)
            else:
                if mark_start:
                    put(bound, cds[-1][1] - 1, 1)
                if mark_stop:
                    put(bound, cds[0][0], 2)
        for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
            if b1 >= a2 or a2 - b1 < min_intron:
                continue
            if s == 1:
                put(bound, b1, 3); put(bound, a2 - 1, 4)
            else:
                put(bound, a2 - 1, 3); put(bound, b1, 4)
    return label, frame, strand, bound


SITE_KINDS = ("start", "stop", "donor", "acceptor", "cds_donor", "cds_acceptor", "short_gap")


def distinct_sites(transcripts: list[dict], min_intron: int = MIN_INTRON) -> set[tuple[str, int, int]]:
    """Every distinct (kind, strand, position) the transcripts state, with
    the positions and kinds of the ``bound`` channel: start is the first
    base of the start codon, stop the last base of the stop codon, donor
    the first intron base and acceptor the last, all at + strand
    coordinates.  A donor or acceptor whose intron lies between two CDS
    segments is counted again as ``cds_donor`` / ``cds_acceptor``, which
    is the count a transcript-level CDS score sees.  Two isoforms sharing a
    site contribute it once, so the difference between the set over all
    isoforms and the set over the painted ones is what an isoform policy
    discards.  An exon gap shorter than ``min_intron`` contributes one
    ``short_gap`` at its first base and no donor or acceptor."""
    sites: set[tuple[str, int, int]] = set()
    for t in transcripts:
        s = 1 if t.get("strand", "+") == "+" else -1
        cds = sorted(t.get("cds", []))
        if cds:
            st = cds[0][0] if s == 1 else cds[-1][1] - 1
            sp = cds[-1][1] - 1 if s == 1 else cds[0][0]
            if not t.get("_partial5"):
                sites.add(("start", s, st))
            if not t.get("_partial3"):
                sites.add(("stop", s, sp))
        cds_ends = {b for _, b in cds}
        cds_starts = {a for a, _ in cds}
        ex = sorted(t.get("exons", []))
        for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
            if b1 >= a2:
                continue
            if a2 - b1 < min_intron:
                sites.add(("short_gap", s, b1))
                continue
            donor, acceptor = (b1, a2 - 1) if s == 1 else (a2 - 1, b1)
            sites.add(("donor", s, donor)); sites.add(("acceptor", s, acceptor))
            if b1 in cds_ends and a2 in cds_starts:
                sites.add(("cds_donor", s, donor)); sites.add(("cds_acceptor", s, acceptor))
    return sites


def site_counts(sites: set[tuple[str, int, int]], lo: int, hi: int) -> dict[str, int]:
    return {k: sum(1 for kind, _, pos in sites if kind == k and lo <= pos < hi) for k in SITE_KINDS}


def feature_lengths(transcripts: list[dict], lo: int, hi: int, min_intron: int = MIN_INTRON) -> dict:
    """Lengths of the introns, CDSs and transcript spans the transcripts
    state, with counts under ``LENGTH_FLOORS``.  Introns are distinct
    (strand, start, end) intervals between consecutive exons and are
    measured only when both ends lie inside [lo, hi); the ones cut by the
    window edge are counted as ``clipped`` and not measured.  A CDS length
    is a transcript's coding bases summed over its whole annotation
    record, which the fetcher carries complete even where the window
    cuts it (``outside_window`` counts the transcripts whose CDS reaches
    past the edge); a span is txStart to txEnd the same way.  Two isoforms
    sharing an intron contribute it once; CDS and span are per transcript.
    ``introns`` counts every exon gap the annotation states, the ones under
    ``min_intron`` included (``below_min_intron``), so the distribution is
    the annotation's and the floor is applied on top of it."""
    introns: set[tuple[int, int, int]] = set()
    clipped = 0
    cds_len, spans, cds_out, span_out = [], [], 0, 0
    for t in transcripts:
        s = 1 if t.get("strand", "+") == "+" else -1
        ex = sorted(t.get("exons", []))
        for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
            if b1 >= a2:
                continue
            if lo <= b1 and a2 <= hi:
                introns.add((s, b1, a2))
            else:
                clipped += 1
        cds = sorted(t.get("cds", []))
        if cds:
            cds_len.append(sum(b - a for a, b in cds))
            cds_out += int(cds[0][0] < lo or cds[-1][1] > hi)
        if ex:
            spans.append(ex[-1][1] - ex[0][0])
            span_out += int(ex[0][0] < lo or ex[-1][1] > hi)

    def summary(vals: list[int], kind: str) -> dict:
        vals = sorted(vals)
        n = len(vals)
        med = None if not n else (vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2)
        out = {"n": n, "min": vals[0] if n else None, "median": med, "max": vals[-1] if n else None}
        for f in LENGTH_FLOORS[kind]:
            out[f"below_{f}"] = sum(1 for v in vals if v < f)
        return out

    ilen = [b - a for _, a, b in introns]
    return {
        "units": {"lengths": "bases", "n": "distinct introns / transcripts", "below_F": "count with length < F",
                  "below_min_intron": "exon gaps shorter than min_intron, painted short_gap not intron"},
        "min_intron": min_intron,
        "introns": {**summary(ilen, "intron"), "clipped": clipped,
                    "below_min_intron": sum(1 for v in ilen if v < min_intron)},
        "cds": {**summary(cds_len, "cds"), "outside_window": cds_out},
        "span": {**summary(spans, "span"), "outside_window": span_out},
    }


DONORS = ("GT", "GC")


MOTIF_CLASSES = ("exact", "borrows_exon_base", "donor_only", "acceptor_only", "neither", "ambiguous")


def motif_windows(left: str, gap: str, right: str) -> dict:
    """engels' overlapping-window test (relay note 20260909T212826Z-engels-0021):
    with one exonic base on each side, ``s = left[-1] + gap + right[0]``, a
    GT/GC donor window ``s[1:3]`` and an AG acceptor window ``s[L-1:L+1]``
    can both be satisfied by a 1-base gap that borrows a base from each
    exon (A|G|T reads GT and AG), which is how bricks2marble's default
    motif mask admits an ``EI -> I -> IE`` path of one intron base.

    The two windows are reported apart, because their conjunction alone
    cannot say which side failed (engels, note 20260909T222452Z-engels-0022):
    ``donor`` and ``acceptor`` are True, False, or None when the window
    holds a base outside A/C/G/T (an N in the assembly), and ``class`` is
    one of ``MOTIF_CLASSES``: ``exact`` (both pass on the gap's own bases,
    needs 4), ``borrows_exon_base`` (both pass, a window reached into an
    exon, which is the 1-base case), ``donor_only``, ``acceptor_only``,
    ``neither`` (both decided and both fail), ``ambiguous`` (a window is
    unresolved, whatever the other window says: an unresolved window next
    to a failed one is ``ambiguous``, not ``neither``, which is reserved
    for two decided failures; engels, note 20260909T232517Z-engels-0023).
    ``both`` is the conjunction with an unresolved window counted as not
    passing, which is what ``motif_window`` returned before 0.10.

    The windows are two bases each, so an N inside a gap of three bases or
    more can sit outside both and leave every window resolved (A|GTNAG|C is
    ``exact``); ``motif_fields`` therefore also counts the gap's own bases
    outside A/C/G/T as ``gap_unresolved_bases``, and a record's
    ``gaps_with_unresolved_bases`` is what certifies the gap interior,
    which ``motif_unresolved_windows`` does not."""
    L = len(gap)
    if L < 1 or not left or not right:
        return {"donor": None, "acceptor": None, "both": False, "borrows": False, "class": "ambiguous"}
    s = (left[-1] + gap + right[0]).upper()
    dw, aw = s[1:3], s[L - 1:L + 1]
    donor = (dw in DONORS) if set(dw) <= set("ACGT") else None
    acceptor = (aw == "AG") if set(aw) <= set("ACGT") else None
    both = donor is True and acceptor is True
    borrows = both and L < 2
    if both:
        cls = "borrows_exon_base" if borrows else "exact"
    elif donor is False and acceptor is False:
        cls = "neither"
    elif donor is True and acceptor is False:
        cls = "donor_only"
    elif donor is False and acceptor is True:
        cls = "acceptor_only"
    else:
        cls = "ambiguous"
    return {"donor": donor, "acceptor": acceptor, "both": both, "borrows": borrows, "class": cls}


def motif_window(left: str, gap: str, right: str) -> tuple[bool, bool]:
    """(both windows satisfied, a window borrowed an exon base): the pre-0.10
    conjunction, kept for callers that only need the combined verdict."""
    m = motif_windows(left, gap, right)
    return m["both"], m["borrows"]


def motif_fields(left: str, gap: str, right: str) -> dict:
    """The per-gap motif fields of a ``short_gaps`` entry from its flanks."""
    m = motif_windows(left, gap, right)
    g = gap.upper()
    return {"motif_exact": len(g) >= 4 and g[:2] in DONORS and g[-2:] == "AG",
            "motif_window": m["both"], "motif_borrows_exon_base": m["borrows"],
            "motif_donor": m["donor"], "motif_acceptor": m["acceptor"], "motif_class": m["class"],
            "gap_unresolved_bases": sum(1 for c in g if c not in "ACGT")}


MOTIF_FIELDS_NONE = {"motif_exact": None, "motif_window": None, "motif_borrows_exon_base": None,
                     "motif_donor": None, "motif_acceptor": None, "motif_class": None,
                     "gap_unresolved_bases": None}


def motif_summary(gaps: list[dict]) -> dict:
    """Counts over ``short_gaps`` entries: the three pre-0.10 counters, the
    class table (every class listed, zero when absent), the unresolved
    windows, and the gaps holding a base outside A/C/G/T anywhere (0.11)."""
    return {
        "motif_exact": sum(1 for g in gaps if g["motif_exact"]),
        "motif_window": sum(1 for g in gaps if g["motif_window"]),
        "motif_borrows_exon_base": sum(1 for g in gaps if g["motif_borrows_exon_base"]),
        "motif_by_class": {c: sum(1 for g in gaps if g["motif_class"] == c) for c in MOTIF_CLASSES},
        "motif_unresolved_windows": sum((g["motif_donor"] is None) + (g["motif_acceptor"] is None)
                                        for g in gaps if g["motif_class"] is not None),
        "gaps_with_unresolved_bases": sum(1 for g in gaps if g["gap_unresolved_bases"]),
    }


def short_gaps(transcripts: list[dict], seq: str, start: int, end: int,
               min_intron: int = MIN_INTRON, flank: int = 2) -> dict:
    """Every distinct exon gap shorter than ``min_intron`` the transcripts
    state, with what a decoder or a scorer would need to classify it: its
    length and length mod 3, whether both flanking exon ends are CDS ends
    (``in_cds``: the shape of a RefSeq programmed frameshift, two CDS
    blocks and a gap, as against a short intron on the UTR side), the
    transcripts stating it, and, where the gap lies inside the window,
    ``flank`` exon bases on each side and the gap itself in transcript
    orientation, whether the gap's own ends read GT/GC..AG (``motif_exact``,
    needs 4 bases) and what the overlapping windows of ``motif_windows`` say,
    donor and acceptor apart, with ``motif_class`` naming the failing side
    or an unresolved base.  Coordinates are + strand, 0-based half-open.  ``exception`` is
    carried when the fetcher recorded one (NCBI GFF3 ``exception=ribosomal
    slippage``, stalin's note 20260909T214050Z-stalin-0021); the UCSC
    genePred sources carry none, so ``in_cds`` and the flanks are what the
    record can say."""
    gaps: dict[tuple[int, int, int], dict] = {}
    for t in transcripts:
        s = 1 if t.get("strand", "+") == "+" else -1
        ex = sorted(t.get("exons", []))
        cds_ends = {b for _, b in t.get("cds", [])}
        cds_starts = {a for a, _ in t.get("cds", [])}
        for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
            if not 0 < a2 - b1 < min_intron:
                continue
            g = gaps.setdefault((s, b1, a2), {
                "strand": "+" if s == 1 else "-", "start": b1, "end": a2, "length": a2 - b1,
                "length_mod_3": (a2 - b1) % 3, "in_cds": False, "transcripts": [], "exception": None})
            g["in_cds"] = g["in_cds"] or (b1 in cds_ends and a2 in cds_starts)
            if t.get("id") not in g["transcripts"]:
                g["transcripts"].append(t.get("id"))
            if t.get("exception") and not g["exception"]:
                g["exception"] = t["exception"]
    out = []
    for key in sorted(gaps):
        g = gaps[key]
        b1, a2 = g["start"], g["end"]
        if start <= b1 - flank and a2 + flank <= end:
            left, gap, right = seq[b1 - flank - start:b1 - start], seq[b1 - start:a2 - start], seq[a2 - start:a2 + flank - start]
            if g["strand"] == "-":
                left, gap, right = right.translate(COMP)[::-1], gap.translate(COMP)[::-1], left.translate(COMP)[::-1]
            g.update({"left": left, "gap": gap, "right": right, **motif_fields(left, gap, right)})
        else:
            g.update({"left": None, "gap": None, "right": None, **MOTIF_FIELDS_NONE})
        out.append(g)
    lengths: dict[str, int] = {}
    for g in out:
        lengths[str(g["length"])] = lengths.get(str(g["length"]), 0) + 1
    return {
        "units": {"n": "distinct (strand, start, end) exon gaps shorter than min_intron",
                  "flanks": f"{flank} exon bases each side, transcript orientation",
                  "motif_window": "engels' overlapping-window test, one exon base borrowed each side; "
                                  "the conjunction of motif_donor and motif_acceptor, an unresolved window not passing",
                  "motif_class": "exact | borrows_exon_base | donor_only | acceptor_only | neither | ambiguous "
                                 "(a window holds a base outside ACGT, whatever the other window says; "
                                 "neither means two decided failures)",
                  "gap_unresolved_bases": "bases of the gap itself outside ACGT, inside or outside the two-base windows; "
                                          "gaps_with_unresolved_bases counts the gaps where it is nonzero"},
        "min_intron": min_intron, "n": len(out),
        "in_cds": sum(1 for g in out if g["in_cds"]),
        "by_length": lengths,
        **motif_summary(out),
        "outside_window": sum(1 for g in out if g["gap"] is None),
        "gaps": out,
    }


# ---------------------------------------------------------------- cutting

def revcomp_example(ex: dict, L: int, K: int) -> dict:
    out = dict(ex)
    out["ref"] = bytes(4 if b == 4 else 3 - b for b in reversed(ex["ref"]))
    out["mask"] = ex["mask"][::-1]
    out["label"] = ex["label"][::-1]
    out["frame"] = ex["frame"][::-1]
    out["strand"] = bytes((0 if v == 0 else (255 if v == 1 else 1)) for v in reversed(ex["strand"]))
    # donors and acceptors keep their meaning in transcript orientation; only
    # the coordinate flips.  Starts and stops likewise.
    out["bound"] = ex["bound"][::-1]
    out["cons"] = ex["cons"][::-1]
    inf = bytearray()
    ins = bytearray()
    for k in range(K):
        row = ex["inf"][k * L:(k + 1) * L]
        inf += bytes(b if b >= 4 else 3 - b for b in reversed(row))
        # an insertion after base i becomes an insertion after base i-1's mirror
        irow = ex["ins"][k * L:(k + 1) * L]
        shifted = bytearray(L)
        for i in range(L):
            if irow[i] and i + 1 < L:
                shifted[L - 1 - (i + 1)] = irow[i]
        ins += shifted
    out["inf"], out["ins"] = bytes(inf), bytes(ins)
    return out


def cut(stem: str, out_dir: str, L: int, S: int, drop: set[str], both: bool, quiet: bool,
        transcript_types: str = "benchmark", drop_ancestors: bool = False,
        reference_anchored: bool = False, isoforms: str = "union",
        representatives: set[str] | None = None, partial_ends: str = "both",
        genetic_code: int = 1, stop_codons: set[str] | None = None,
        duplicate_rows: str = "identity", min_intron: int = MIN_INTRON) -> list[str]:
    stops = set(stop_codons) if stop_codons else GENETIC_CODE_STOPS.get(genetic_code)
    if stops is None:
        raise SystemExit(f"genetic code {genetic_code} is not in the table; pass --stop-codons")
    with open(stem + ".manifest.json") as fh:
        manifest = json.load(fh)
    with open(stem + ".manifest.json", "rb") as fh:
        manifest_sha = hashlib.sha256(fh.read()).hexdigest()
    ref = manifest.get("assembly") or manifest.get("species")
    if manifest.get("source") == "ensembl":
        ref = manifest.get("species", ref)
    win = manifest["window"]
    chrom, start, end = win["chrom"], int(win["start"]), int(win["end"])
    with open(stem + ".fa") as fh:
        seq = "".join(ln.strip() for ln in fh if not ln.startswith(">")).upper()
    if len(seq) != end - start:
        raise SystemExit(f"{stem}: FASTA length {len(seq)} != window {end - start}")
    with open(stem + ".maf") as fh:
        blocks = maf_blocks(fh.read())
    trees = []
    if os.path.exists(stem + ".nh"):
        with open(stem + ".nh") as fh:
            trees = [t.strip() + ";" for t in fh.read().split(";") if t.strip()]
    ref_prefix = ref.split(".")[0] if ref else None
    if ref_prefix is None:
        raise SystemExit("cannot determine the reference name from the manifest")
    sources: list[str] = []
    for b in blocks:
        for r in b:
            nm = maf_source(r[1])
            if nm not in sources:
                sources.append(nm)
    # ancestral rows: what the fetcher recorded, else by name
    declared = set(manifest.get("ancestral_sequences") or [])

    def is_ancestor(nm: str) -> bool:
        return nm in declared or ANCESTOR_RE.match(nm) is not None or nm.lower().startswith("ancestor")

    names = set(sources) | {ref_prefix}
    dist: dict[str, float] = {}
    clades: dict[str, set[str]] = {}
    tree_leaves: list[list[str]] = []
    for t in trees:
        lv, d = newick_leaves_and_distances(t, ref_prefix, names)
        tree_leaves.append(lv)
        for k, v in d.items():
            dist.setdefault(k, v)
        for k, v in newick_clades(t, names).items():
            clades.setdefault(k, set()).update(v)
    # rows: species-set members (Ensembl), else the single track tree's leaf
    # order, else the observed species sorted; then out-of-set species in
    # order of appearance; then ancestors sorted by name (see the docstring)
    members = manifest.get("species_set_members")
    if members:
        leaves, row_order = sorted(members), "species-set"
    elif len(trees) == 1:
        leaves, row_order = tree_leaves[0], "tree"
    else:
        leaves, row_order = sorted(s for s in sources if not is_ancestor(s)), "observed-sorted"
    leaves = [n for n in leaves if n != ref_prefix and n.split(".")[0] != ref_prefix and not is_ancestor(n)]
    ancestors = sorted(s for s in sources if is_ancestor(s) and s != ref_prefix)
    others = [s for s in sources if s not in leaves and s not in ancestors and s != ref_prefix]
    rows = leaves + others + ancestors
    dropped = [r for r in rows if r in drop]
    dropped_anc = [a for a in ancestors if a not in drop
                   and (drop_ancestors or bool(clades.get(a, set()) & drop))]
    rows = [r for r in rows if r not in drop and r not in dropped_anc]
    K = len(rows)
    aln_class, aln_how = classify_alignment(manifest, reference_anchored)
    inf, ins, block_stats = paint_informants(blocks, ref_prefix, start, end, rows, duplicate_rows)
    with open(stem + ".annotation.json") as fh:
        ann = json.load(fh)
    transcripts = ann.get("transcripts", ann if isinstance(ann, list) else [])
    transcripts, excluded = select_transcripts(transcripts, transcript_types)
    cds_ends = cds_end_status(transcripts, seq, start, end, stops, partial_ends)
    sites_all = distinct_sites(transcripts, min_intron)
    lengths_all = feature_lengths(transcripts, start, end, min_intron)
    gaps_all = short_gaps(transcripts, seq, start, end, min_intron)
    transcripts, dropped_iso, fallback_loci = select_isoforms(transcripts, isoforms, representatives)
    sites_kept = distinct_sites(transcripts, min_intron)
    lengths_kept = feature_lengths(transcripts, start, end, min_intron)
    gaps_kept = short_gaps(transcripts, seq, start, end, min_intron)
    label, frame, strand, bound = paint_labels(start, end, transcripts, min_intron)
    cons = [math.nan] * (end - start)
    if os.path.exists(stem + ".conservation.json"):
        with open(stem + ".conservation.json") as fh:
            c = json.load(fh)
        vals = c.get("values") if isinstance(c, dict) else None
        if isinstance(vals, list) and vals and isinstance(vals[0], (int, float)) and len(vals) == end - start:
            cons = [float(v) if v is not None else math.nan for v in vals]
        elif isinstance(vals, list):
            for item in vals:  # sparse [pos, value] or {start,end,value}
                if isinstance(item, dict):
                    for q in range(max(int(item["start"]), start), min(int(item["end"]), end)):
                        cons[q - start] = float(item.get("value", math.nan))
                elif isinstance(item, (list, tuple)) and len(item) == 3:  # [start, end, value]
                    for q in range(max(int(item[0]), start), min(int(item[1]), end)):
                        cons[q - start] = float(item[2])
                elif isinstance(item, (list, tuple)) and len(item) == 2:  # [pos, value]
                    q = int(item[0])
                    if start <= q < end:
                        cons[q - start] = float(item[1])
    n = end - start
    if L <= 0:
        L, S = n, n
    os.makedirs(out_dir, exist_ok=True)
    written = []
    base = os.path.basename(stem)
    offsets = list(range(0, max(n - L, 0) + 1, S)) or [0]
    if offsets[-1] + L < n:
        offsets.append(n - L if n >= L else 0)
    for k, off in enumerate(offsets):
        a, b = off, min(off + L, n)
        pad = L - (b - a)
        ex = {
            "ref": bytes(BASE.get(ch, 4) for ch in seq[a:b]) + b"\x04" * pad,
            "mask": b"\x01" * (b - a) + b"\x00" * pad,
            "label": bytes(label[a:b]) + b"\x00" * pad,
            "frame": bytes(frame[a:b]) + b"\xff" * pad,
            "strand": bytes(strand[a:b]) + b"\x00" * pad,
            "bound": bytes(bound[a:b]) + b"\x00" * pad,
            "inf": b"".join(bytes(inf[j][a:b]) + bytes([INF_UNALIGNED]) * pad for j in range(K)),
            "ins": b"".join(bytes(ins[j][a:b]) + b"\x00" * pad for j in range(K)),
            "cons": cons[a:b] + [math.nan] * pad,
        }
        variants = [("", ex)]
        if both:
            variants.append((".rc", revcomp_example(ex, L, K)))
        for suffix, e in variants:
            arrays = {
                "ref": (e["ref"], "|u1", (L,)),
                "mask": (e["mask"], "|u1", (L,)),
                "label": (e["label"], "|u1", (L,)),
                "frame": (e["frame"], "|i1", (L,)),
                "strand": (e["strand"], "|i1", (L,)),
                "bound": (e["bound"], "|u1", (L,)),
                "inf": (e["inf"], "|u1", (K, L)),
                "ins": (e["ins"], "|u1", (K, L)),
                "dist": (f32([dist.get(r, math.nan) for r in rows]), "<f4", (K,)),
                "cons": (f32(e["cons"]), "<f4", (L,)),
            }
            name = f"{base}.w{k}{suffix}"
            write_npz(os.path.join(out_dir, name + ".npz"), arrays)
            side = {
                "tool": "scripts/data/cut_windows.py", "version": TOOL_VERSION,
                "source_manifest_sha256": manifest_sha, "assembly": ref, "chrom": chrom,
                "start": start + a, "end": start + b, "length": L, "padding": pad,
                "orientation": "reverse_complement" if suffix else "reference",
                "alignment_track": (manifest.get("alignment_track") or manifest.get("species_set")
                                    or manifest.get("alignment_access")),
                "alignment_class": aln_class,
                "alignment_class_source": aln_how,
                "informants": rows,
                "row_order": row_order,
                "ancestral": [r for r in rows if is_ancestor(r)],
                "ancestral_clades": {r: sorted(clades.get(r, ())) for r in rows if is_ancestor(r)},
                "dropped_species": dropped,
                "dropped_ancestors": dropped_anc,
                "dropped_rows_only": bool(dropped or dropped_anc) and aln_class != "reference_anchored",
                "block_selection": block_stats,
                "transcripts": [t.get("id") for t in transcripts],
                "transcript_types": transcript_types,
                "transcripts_excluded": excluded,
                "isoform_policy": isoforms,
                "transcripts_dropped_isoforms": dropped_iso,
                "representative_fallback_loci": fallback_loci,
                "distinct_sites": {
                    "all_isoforms": site_counts(sites_all, start + a, start + b),
                    "painted": site_counts(sites_kept, start + a, start + b),
                    "dropped": site_counts(sites_all - sites_kept, start + a, start + b),
                },
                "feature_lengths": {
                    "floors": {k: list(v) for k, v in LENGTH_FLOORS.items()},
                    "all_isoforms": lengths_all,
                    "painted": lengths_kept,
                },
                "min_intron": min_intron,
                "short_gaps": {"all_isoforms": gaps_all, "painted": gaps_kept},
                "genetic_code": genetic_code if not stop_codons else None,
                "stop_codons": sorted(stops),
                "cds_ends": cds_ends,
                "label_counts": {c: sum(1 for v, m in zip(e["label"], e["mask"]) if m and v == i)
                                 for c, i in LABEL.items()},
            }
            with open(os.path.join(out_dir, name + ".json"), "w") as fh:
                json.dump(side, fh, indent=1)
            written.append(name)
            if not quiet:
                aligned = sum(1 for v in e["inf"] if v < INF_UNALIGNED) / max(1, K * (b - a))
                print(f"{name}: {chrom}:{start + a}-{start + b} K={K} pad={pad} "
                      f"labels={side['label_counts']} informant-aligned={aligned:.3f}")
    return written


# -------------------------------------------------------------- self test

def self_test() -> int:
    import tempfile
    d = tempfile.mkdtemp()
    # Checks 1-45 test painting geometry on toy sequences whose CDSs open
    # and close on arbitrary codons, so they run with --partial-ends none;
    # checks 46 onwards test the incomplete-end rules on a fixture with
    # real codons and pass the mode explicitly.
    cut_full = globals()["cut"]

    def cut(*a, **k):
        k.setdefault("partial_ends", "none")
        k.setdefault("min_intron", 0)  # the toy introns are 4 bases; checks 64+ pass the floor
        return cut_full(*a, **k)

    stem = os.path.join(d, "toy")
    # reference toy.chr1 20 bp; one + transcript with two exons, CDS 4..8 and 12..17
    ref = "ACGTACGTACGTACGTACGT"
    with open(stem + ".fa", "w") as fh:
        fh.write(">toy.chr1:0-20\n" + ref + "\n")
    with open(stem + ".maf", "w") as fh:
        fh.write("##maf version=1\n"
                 "a score=0\n"
                 "s toy.chr1 0 10 + 20 ACGTACG-TAC\n"
                 "s spA.c   0 11 + 50 ACGTTCGGTAC\n"
                 "s spB.c   0  8 + 50 AC--ACG-TAC\n\n"
                 "a score=0\n"
                 "s toy.chr1 14 6 + 20 GTACGT\n"
                 "s spA.c   20 6 + 50 GTNCGT\n\n")
    with open(stem + ".nh", "w") as fh:
        fh.write("((toy:0.1,spA:0.2):0.3,spB:0.5,spC:0.9);\n")
    with open(stem + ".annotation.json", "w") as fh:
        json.dump({"transcripts": [{"id": "t1", "strand": "+", "start": 2, "end": 19,
                                    "exons": [[2, 8], [12, 19]], "cds": [[4, 8], [12, 17]]},
                                   {"id": "t2", "strand": "-", "start": 0, "end": 2,
                                    "exons": [[0, 2]], "cds": []},
                                   {"id": "t3", "strand": "+", "start": 17, "end": 20, "type": "pseudogene",
                                    "exons": [[17, 20]], "cds": [[17, 20]]}]}, fh)
    with open(stem + ".manifest.json", "w") as fh:
        json.dump({"assembly": "toy", "source": "ucsc", "alignment_track": "multizToy",
                   "window": {"chrom": "chr1", "start": 0, "end": 20}}, fh)
    names = cut(stem, os.path.join(d, "out"), 0, 0, {"spC"}, True, True)
    assert names == ["toy.w0", "toy.w0.rc"], names
    z = zipfile.ZipFile(os.path.join(d, "out", "toy.w0.npz"))

    def arr(name: str) -> bytes:
        data = z.read(name + ".npy")
        hl = struct.unpack("<H", data[8:10])[0]
        return data[10 + hl:]

    side = json.load(open(os.path.join(d, "out", "toy.w0.json")))
    assert side["informants"] == ["spA", "spB"], side["informants"]
    assert side["dropped_species"] == ["spC"]
    label = list(arr("label"))
    assert label == [2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 3, 2, 2, 0], label  # t2's non-coding exon is UTR-class
    frame = [v - 256 if v > 127 else v for v in arr("frame")]
    assert frame[4:8] == [0, 1, 2, 0] and frame[12:17] == [1, 2, 0, 1, 2], frame
    bound = list(arr("bound"))
    assert bound[4] == 1 and bound[16] == 2 and bound[8] == 3 and bound[11] == 4, bound
    strand = [v - 256 if v > 127 else v for v in arr("strand")]
    assert strand[0] == -1 and strand[5] == 1 and strand[19] == 0, strand
    inf = list(arr("inf"))
    A, B = inf[:20], inf[20:]
    # spA: first 10 bases aligned (T mismatch at 4), 10..13 unaligned, 14..19 aligned with N at 16
    assert A[:10] == [0, 1, 2, 3, 3, 1, 2, 3, 0, 1] and A[10:14] == [5] * 4 and A[16] == 6, A
    assert B[2] == 4 and B[3] == 4 and B[14:] == [5] * 6, B
    ins = list(arr("ins"))
    assert ins[6] == 1 and ins[20 + 6] == 0, ins  # spA has a G inserted after ref base 6
    dist = struct.unpack("<2f", arr("dist"))
    assert abs(dist[0] - 0.3) < 1e-6 and abs(dist[1] - 0.9) < 1e-6, dist
    rc = zipfile.ZipFile(os.path.join(d, "out", "toy.w0.rc.npz"))
    data = rc.read("label.npy"); hl = struct.unpack("<H", data[8:10])[0]
    assert list(data[10 + hl:]) == label[::-1]
    data = rc.read("ref.npy"); hl = struct.unpack("<H", data[8:10])[0]
    assert list(data[10 + hl:]) == [3 - BASE[c] for c in reversed(ref)]
    checks = 12

    def npz(path: str) -> dict:
        zz = zipfile.ZipFile(path)
        out = {}
        for nm in zz.namelist():
            data = zz.read(nm); hl = struct.unpack("<H", data[8:10])[0]
            out[nm[:-4]] = data[10 + hl:]
        return out

    # 13-14: the pseudogene t3 is excluded by default and painted with --transcript-types all
    assert side["transcripts_excluded"] == [{"id": "t3", "type": "pseudogene"}] and label[19] == 0, side["transcripts_excluded"]
    cut(stem, os.path.join(d, "all"), 0, 0, {"spC"}, False, True, transcript_types="all")
    assert list(npz(os.path.join(d, "all", "toy.w0.npz"))["label"])[19] == 3
    checks += 2
    # 15-17: alignment class: the toy track is unknown to the table, so dropping is strict
    assert side["alignment_class"] == "unknown" and side["dropped_rows_only"] is True, side
    cut(stem, os.path.join(d, "ra"), 0, 0, {"spC"}, False, True, reference_anchored=True)
    ra = json.load(open(os.path.join(d, "ra", "toy.w0.json")))
    assert ra["alignment_class"] == "reference_anchored" and ra["dropped_rows_only"] is False, ra
    table = {"multiz124way": "reference_anchored", "multiz470way": "reference_anchored",
             "chainPanTro6": "reference_anchored", "cactus241wayBM": "jointly_inferred",
             "cactus447way": "jointly_inferred", "hprc90way": "jointly_inferred", "someNewTrack": "unknown"}
    got = {t: classify_alignment({"source": "ucsc", "alignment_track": t})[0] for t in table}
    assert got == table, got
    assert classify_alignment({"source": "ensembl", "method": "EPO"})[0] == "jointly_inferred"
    checks += 3
    # 18-20: padded example (L > window): label_counts must agree between the forward and .rc sidecars
    cut(stem, os.path.join(d, "pad"), 32, 32, {"spC"}, True, True)
    fwd = json.load(open(os.path.join(d, "pad", "toy.w0.json")))
    rcs = json.load(open(os.path.join(d, "pad", "toy.w0.rc.json")))
    assert fwd["padding"] == 12 and fwd["label_counts"] == rcs["label_counts"], (fwd["label_counts"], rcs["label_counts"])
    assert fwd["label_counts"] == {"intergenic": 1, "intron": 4, "utr": 6, "cds": 9, "short_gap": 0}, fwd["label_counts"]
    pf, pr = npz(os.path.join(d, "pad", "toy.w0.npz")), npz(os.path.join(d, "pad", "toy.w0.rc.npz"))
    assert pr["mask"] == pf["mask"][::-1] and pr["label"] == pf["label"][::-1] and pf["mask"][20:] == b"\x00" * 12
    checks += 3
    # 21-26: an EPO-shaped window: one tree per block, ancestral rows named after their descendants
    e = os.path.join(d, "epo")
    with open(e + ".fa", "w") as fh:
        fh.write(">1:0-20\n" + ref + "\n")
    with open(e + ".maf", "w") as fh:
        fh.write("##maf version=1\n"
                 "a score=0\n"
                 "s gallus_gallus.1       0 10 + 20 ACGTACGTAC\n"
                 "s meleagris_gallopavo.5 0 10 + 30 ACGTTCGTAC\n"
                 "s Ggal-Mgal[2].0        0 10 + 0  ACGTACGTAC\n\n"
                 "a score=0\n"
                 "s gallus_gallus.1       10 10 + 20 GTACGTACGT\n"
                 "s taeniopygia_guttata.7  0 10 + 30 GTACGTACGT\n"
                 "s Ggal-Tgut[3].0         0 10 + 0  GTACGTACGT\n\n")
    with open(e + ".nh", "w") as fh:
        fh.write("(gallus_gallus_1_0_10[+]:0.1,meleagris_gallopavo_5_0_10[+]:0.1)Ggal-Mgal[2]:0.05;\n"
                 "(gallus_gallus_1_10_10[+]:0.2,taeniopygia_guttata_7_0_10[+]:0.3)Ggal-Tgut[3]:0.1;\n")
    with open(e + ".annotation.json", "w") as fh:
        json.dump({"transcripts": []}, fh)
    members = ["taeniopygia_guttata", "gallus_gallus", "meleagris_gallopavo", "anas_platyrhynchos"]
    man = {"assembly": "gallus_gallus", "source": "ensembl", "species_set": "sauropsids", "method": "EPO",
           "ancestral_sequences": ["Ggal-Mgal[2]", "Ggal-Tgut[3]"], "species_set_members": members,
           "window": {"chrom": "1", "start": 0, "end": 20}}
    with open(e + ".manifest.json", "w") as fh:
        json.dump(man, fh)
    cut(e, os.path.join(d, "epo1"), 0, 0, set(), False, True)
    s1 = json.load(open(os.path.join(d, "epo1", "epo.w0.json")))
    assert s1["informants"] == ["anas_platyrhynchos", "meleagris_gallopavo", "taeniopygia_guttata",
                                "Ggal-Mgal[2]", "Ggal-Tgut[3]"], s1["informants"]
    assert s1["row_order"] == "species-set" and s1["ancestral"] == ["Ggal-Mgal[2]", "Ggal-Tgut[3]"], s1
    assert s1["ancestral_clades"]["Ggal-Tgut[3]"] == ["gallus_gallus", "taeniopygia_guttata"], s1["ancestral_clades"]
    z1 = npz(os.path.join(d, "epo1", "epo.w0.npz"))
    d1 = struct.unpack("<5f", z1["dist"])
    assert math.isnan(d1[0]) and abs(d1[1] - 0.2) < 1e-6 and abs(d1[2] - 0.5) < 1e-6 \
        and abs(d1[3] - 0.1) < 1e-6 and abs(d1[4] - 0.2) < 1e-6, d1
    assert set(z1["inf"][:20]) == {INF_UNALIGNED} and z1["inf"][20 + 4] == 3, "anas unaligned, meleagris T at 4"
    # dropping a held-out leaf also drops the ancestor of its clade, and the class is jointly inferred
    cut(e, os.path.join(d, "epo2"), 0, 0, {"taeniopygia_guttata"}, False, True)
    s2 = json.load(open(os.path.join(d, "epo2", "epo.w0.json")))
    assert s2["informants"] == ["anas_platyrhynchos", "meleagris_gallopavo", "Ggal-Mgal[2]"], s2["informants"]
    assert s2["dropped_species"] == ["taeniopygia_guttata"] and s2["dropped_ancestors"] == ["Ggal-Tgut[3]"] \
        and s2["alignment_class"] == "jointly_inferred" and s2["dropped_rows_only"] is True, s2
    # a manifest from before species_set_members / ancestral_sequences: names still classify, order is observed-sorted
    with open(e + ".manifest.json", "w") as fh:
        json.dump({k: v for k, v in man.items() if k not in ("species_set_members", "ancestral_sequences")}, fh)
    cut(e, os.path.join(d, "epo3"), 0, 0, set(), False, True, drop_ancestors=True)
    s3 = json.load(open(os.path.join(d, "epo3", "epo.w0.json")))
    assert s3["informants"] == ["meleagris_gallopavo", "taeniopygia_guttata"] and s3["row_order"] == "observed-sorted" \
        and s3["dropped_ancestors"] == ["Ggal-Mgal[2]", "Ggal-Tgut[3]"], s3
    checks += 7
    # GenArk source names carry a version dot and their tree leaf spells it with a v (mm39 35-way)
    assert maf_source("GCF_003668045.3.NC_048596.1") == "GCF_003668045v3" and maf_source("mm39.chr19") == "mm39" \
        and maf_source("C_sp38_MB_2015.chrI") == "C_sp38_MB_2015" and canonical_name("GCF_003668045.3.NC_048596.1") == "GCF_003668045v3"
    checks += 1
    # 29-33: a gene nested in another gene's intron on the opposite strand owns its
    # bases (strand, frame) even though the enclosing transcript is listed first
    nst = os.path.join(d, "nest")
    for ext in (".fa", ".maf", ".nh", ".manifest.json"):
        shutil.copyfile(stem + ext, nst + ext)
    with open(nst + ".annotation.json", "w") as fh:
        json.dump({"transcripts": [
            {"id": "outer", "strand": "-", "start": 0, "end": 20, "exons": [[0, 3], [17, 20]], "cds": [[0, 3], [17, 20]]},
            {"id": "inner", "strand": "+", "start": 5, "end": 15, "exons": [[5, 8], [11, 15]], "cds": [[6, 8], [11, 14]]},
        ]}, fh)
    cut(nst, os.path.join(d, "nest1"), 0, 0, set(), True, True)
    z2 = npz(os.path.join(d, "nest1", "nest.w0.npz"))
    lab2 = list(z2["label"])
    st2 = [v - 256 if v > 127 else v for v in z2["strand"]]
    fr2 = [v - 256 if v > 127 else v for v in z2["frame"]]
    bd2 = list(z2["bound"])
    assert lab2 == [3, 3, 3, 1, 1, 2, 3, 3, 1, 1, 1, 3, 3, 3, 2, 1, 1, 3, 3, 3], lab2
    assert st2[:5] == [-1] * 5 and st2[5:15] == [1] * 10 and st2[15:] == [-1] * 5, st2
    assert fr2[:3] == [2, 1, 0] and fr2[17:] == [2, 1, 0] and fr2[6:8] == [0, 1] and fr2[11:14] == [2, 0, 1] \
        and fr2[3:6] == [-1] * 3 and fr2[8:11] == [-1] * 3, fr2
    assert bd2[19] == 1 and bd2[0] == 2 and bd2[16] == 3 and bd2[3] == 4, bd2  # outer, read on -
    assert bd2[6] == 1 and bd2[13] == 2 and bd2[8] == 3 and bd2[10] == 4, bd2  # inner, read on +
    checks += 5
    # 34: the reverse complement of the nested case flips strand signs and reverses everything else
    r2 = npz(os.path.join(d, "nest1", "nest.w0.rc.npz"))
    assert [v - 256 if v > 127 else v for v in r2["strand"]] == [-v for v in st2][::-1] \
        and list(r2["label"]) == lab2[::-1] and list(r2["frame"]) == list(z2["frame"])[::-1]
    checks += 1
    # 35: the same inputs give byte-identical archives (fixed zip timestamps), so checksums mean something
    cut(nst, os.path.join(d, "nest2"), 0, 0, set(), True, True)
    for nm in ("nest.w0.npz", "nest.w0.rc.npz"):
        with open(os.path.join(d, "nest1", nm), "rb") as f1, open(os.path.join(d, "nest2", nm), "rb") as f2:
            assert f1.read() == f2.read(), nm
    checks += 1
    # 36-41: isoform policy.  Gene G has two + isoforms: ta (9 CDS bases, exons
    # 2-8 and 12-19) and tb (8 CDS bases, exons 2-8 and 10-19, so bases 10-11
    # are CDS in tb and intronic in ta, and 14-16 are UTR in tb, CDS in ta).
    # Two unnamed - transcripts overlap at 0-3 and cluster into one locus by span.
    iso = os.path.join(d, "iso")
    for ext in (".fa", ".maf", ".nh", ".manifest.json"):
        shutil.copyfile(stem + ext, iso + ext)
    with open(iso + ".annotation.json", "w") as fh:
        json.dump({"transcripts": [
            {"id": "ta.1", "gene": "G", "strand": "+", "start": 2, "end": 19, "exons": [[2, 8], [12, 19]], "cds": [[4, 8], [12, 17]]},
            {"id": "tb.2", "gene": "G", "strand": "+", "start": 2, "end": 19, "exons": [[2, 8], [10, 19]], "cds": [[4, 8], [10, 14]]},
            {"id": "tc", "strand": "-", "start": 0, "end": 2, "exons": [[0, 2]], "cds": []},
            {"id": "td", "strand": "-", "start": 0, "end": 3, "exons": [[0, 3]], "cds": [[0, 3]]},
        ]}, fh)
    cut(iso, os.path.join(d, "iso_u"), 0, 0, set(), False, True)
    su = json.load(open(os.path.join(d, "iso_u", "iso.w0.json")))
    lu = list(npz(os.path.join(d, "iso_u", "iso.w0.npz"))["label"])
    assert su["isoform_policy"] == "union" and su["transcripts_dropped_isoforms"] == [] \
        and su["transcripts"] == ["ta.1", "tb.2", "tc", "td"], su
    assert lu[10:12] == [3, 3] and lu[14:17] == [3, 3, 3] and lu[0:3] == [3, 3, 3], lu  # union: tb's CDS and ta's CDS both paint
    checks += 1
    # distinct sites over all isoforms: one shared start (4) and donor (8), two
    # acceptors (11 and 9) and two + stops (16, 13), both introns CDS-flanked;
    # td adds a - start at 2 and stop at 0; tc has no CDS and no intron
    all_sites = {"start": 2, "stop": 3, "donor": 1, "acceptor": 2, "cds_donor": 1, "cds_acceptor": 2, "short_gap": 0}
    assert su["distinct_sites"] == {"all_isoforms": all_sites, "painted": all_sites,
                                    "dropped": dict.fromkeys(SITE_KINDS, 0)}, su["distinct_sites"]
    checks += 1
    cut(iso, os.path.join(d, "iso_l"), 0, 0, set(), False, True, isoforms="longest-cds")
    sl = json.load(open(os.path.join(d, "iso_l", "iso.w0.json")))
    ll = list(npz(os.path.join(d, "iso_l", "iso.w0.npz"))["label"])
    assert sl["transcripts"] == ["ta.1", "td"] and sl["representative_fallback_loci"] == [], sl["transcripts"]
    assert sl["transcripts_dropped_isoforms"] == [{"id": "tb.2", "locus": "G", "reason": "shorter-cds"},
                                                  {"id": "tc", "locus": "-0-3", "reason": "shorter-cds"}], sl["transcripts_dropped_isoforms"]
    assert ll[10:12] == [1, 1] and ll[14:17] == [3, 3, 3] and ll[0:3] == [3, 3, 3], ll
    checks += 2
    # longest-cds keeps ta: tb's acceptor at 9 and stop at 13 are the loss, the shared start and donor are not
    assert sl["distinct_sites"]["all_isoforms"] == all_sites and sl["distinct_sites"]["painted"] == \
        {"start": 2, "stop": 2, "donor": 1, "acceptor": 1, "cds_donor": 1, "cds_acceptor": 1, "short_gap": 0} and \
        sl["distinct_sites"]["dropped"] == {"start": 0, "stop": 1, "donor": 0, "acceptor": 1,
                                            "cds_donor": 0, "cds_acceptor": 1, "short_gap": 0}, sl["distinct_sites"]
    checks += 1
    # a UTR-only alternative intron is a donor/acceptor but not a cds_donor/cds_acceptor
    assert site_counts(distinct_sites([{"strand": "+", "start": 0, "end": 20, "exons": [[0, 4], [8, 20]], "cds": [[10, 16]]}], min_intron=0), 0, 20) == \
        {"start": 1, "stop": 1, "donor": 1, "acceptor": 1, "cds_donor": 0, "cds_acceptor": 0, "short_gap": 0}
    checks += 1
    # a representative list names tb without its version; the unnamed - locus falls back to longest-cds
    cut(iso, os.path.join(d, "iso_r"), 0, 0, set(), False, True, isoforms="representative", representatives={"tb"})
    sr = json.load(open(os.path.join(d, "iso_r", "iso.w0.json")))
    lr = list(npz(os.path.join(d, "iso_r", "iso.w0.npz"))["label"])
    assert sr["transcripts"] == ["tb.2", "td"] and sr["representative_fallback_loci"] == ["-0-3"], sr
    assert sr["transcripts_dropped_isoforms"] == [{"id": "ta.1", "locus": "G", "reason": "not-representative"},
                                                  {"id": "tc", "locus": "-0-3", "reason": "shorter-cds"}], sr["transcripts_dropped_isoforms"]
    assert lr[10:12] == [3, 3] and lr[14:17] == [2, 2, 2], lr
    checks += 2
    assert sr["distinct_sites"]["dropped"] == {"start": 0, "stop": 1, "donor": 0, "acceptor": 1,
                                               "cds_donor": 0, "cds_acceptor": 1, "short_gap": 0}, sr["distinct_sites"]
    checks += 1
    # the id file reader: comments, blank lines, first column of a TSV, versioned ids
    with open(os.path.join(d, "reps.tsv"), "w") as fh:
        fh.write("# MANE-like list\n\nta.1\tG\tMANE Select\n")
    assert read_id_list(os.path.join(d, "reps.tsv")) == {"ta.1"}
    cut(iso, os.path.join(d, "iso_r2"), 0, 0, set(), False, True, isoforms="representative",
        representatives=read_id_list(os.path.join(d, "reps.tsv")))
    assert json.load(open(os.path.join(d, "iso_r2", "iso.w0.json")))["transcripts"] == ["ta.1", "td"]
    checks += 1
    # 46-55: incomplete CDS ends.  48 bp reference with four single-exon
    # CDSs: tA (+, 2-14, ATG AAA CCC TAG, declared complete both ends), tB
    # (+, 16-25, CCC AAA TGA, nothing declared), tC (-, 27-36, reads CCC AAA
    # TAA on its strand, declared 5' incomplete with frame 1), tD (+, 38-47,
    # GGG CCC AAA, declared complete at both ends against the sequence).
    pe = os.path.join(d, "pe")
    pref = "GG" + "ATGAAACCCTAG" + "GG" + "CCCAAATGA" + "GG" + "TTAAAAGGG" + "GG" + "GGGCCCAAA" + "G"
    assert len(pref) == 48
    with open(pe + ".fa", "w") as fh:
        fh.write(">pe.chr1:0-48\n" + pref + "\n")
    with open(pe + ".maf", "w") as fh:
        fh.write("##maf version=1\na score=0\n"
                 f"s pe.chr1 0 48 + 48 {pref}\n"
                 f"s spA.c   0 48 + 48 {pref}\n\n")
    with open(pe + ".nh", "w") as fh:
        fh.write("(pe:0.1,spA:0.2);\n")
    with open(pe + ".manifest.json", "w") as fh:
        json.dump({"assembly": "pe", "source": "ucsc", "alignment_track": "multizToy",
                   "window": {"chrom": "chr1", "start": 0, "end": 48}}, fh)
    with open(pe + ".annotation.json", "w") as fh:
        json.dump({"transcripts": [
            {"id": "tA", "strand": "+", "start": 2, "end": 14, "exons": [[2, 14]], "cds": [[2, 14]],
             "cds_start_status": "complete", "cds_end_status": "complete"},
            {"id": "tB", "strand": "+", "start": 16, "end": 25, "exons": [[16, 25]], "cds": [[16, 25]]},
            {"id": "tC", "strand": "-", "start": 27, "end": 36, "exons": [[27, 36]], "cds": [[27, 36]],
             "cds_start_status": "incomplete", "cds_start_frame": 1},
            {"id": "tD", "strand": "+", "start": 38, "end": 47, "exons": [[38, 47]], "cds": [[38, 47]],
             "cds_start_status": "complete", "cds_end_status": "complete"},
        ]}, fh)
    # 46: default (both): declaration where stated, sequence where not
    cut_full(pe, os.path.join(d, "pe_b"), 0, 0, set(), True, True)
    sb = json.load(open(os.path.join(d, "pe_b", "pe.w0.json")))
    zb = npz(os.path.join(d, "pe_b", "pe.w0.npz"))
    ce = sb["cds_ends"]
    assert ce["policy"] == "both" and ce["partial_5prime"] == ["tB", "tC"] and ce["partial_3prime"] == [], ce
    checks += 1
    # 47: marks follow the verdicts: tB and tC lose their start mark, everything else is marked
    bb = list(zb["bound"])
    assert bb[2] == 1 and bb[13] == 2 and bb[16] == 0 and bb[24] == 2 and bb[35] == 0 and bb[27] == 2 \
        and bb[38] == 1 and bb[46] == 2, bb
    assert sb["distinct_sites"]["all_isoforms"]["start"] == 2 and sb["distinct_sites"]["all_isoforms"]["stop"] == 4
    checks += 1
    # 48: the declared frame offsets tC's frame channel (read on -, first CDS base at 35 is codon position 1)
    fb = [v - 256 if v > 127 else v for v in zb["frame"]]
    assert fb[35] == 1 and fb[34] == 2 and fb[33] == 0 and fb[27] == 0 and ce["frame_offsets"] == {"tC": 1}, fb
    assert fb[2:5] == [0, 1, 2] and fb[16] == 0, fb  # complete and undeclared 5' ends start at 0
    checks += 1
    # 49: evidence counts and the two tD disagreements are on record, decided by the declaration
    assert ce["declared"] == {"start_complete": 2, "start_incomplete": 1, "start_undeclared": 1,
                              "stop_complete": 2, "stop_incomplete": 0, "stop_undeclared": 2}, ce["declared"]
    assert ce["sequence"] == {"start_atg": 1, "start_near_cognate": 0, "start_other": 3, "start_outside_window": 0,
                              "stop_in_cds": 3, "stop_after_cds": 0, "stop_none": 1, "stop_outside_window": 0,
                              "cds_length_not_multiple_of_3": 0}, ce["sequence"]
    assert [(x["id"], x["end"], x["declared"], x["sequence"], x["decided"]) for x in ce["disagreements"]] == \
        [("tD", "start", "complete", "other", "declared"), ("tD", "stop", "complete", "none", "declared")], ce
    checks += 1
    # 50: the reverse-complement example carries the same record and mirrored marks
    sr = json.load(open(os.path.join(d, "pe_b", "pe.w0.rc.json")))
    assert sr["cds_ends"] == ce and list(npz(os.path.join(d, "pe_b", "pe.w0.rc.npz"))["bound"]) == bb[::-1]
    checks += 1
    # 51: none paints every end, as before 0.5, and applies no frame offset
    cut_full(pe, os.path.join(d, "pe_n"), 0, 0, set(), False, True, partial_ends="none")
    sn = json.load(open(os.path.join(d, "pe_n", "pe.w0.json")))
    zn = npz(os.path.join(d, "pe_n", "pe.w0.npz"))
    assert sn["cds_ends"]["partial_5prime"] == [] and list(zn["bound"])[16] == 1 and list(zn["bound"])[35] == 1 \
        and sn["distinct_sites"]["all_isoforms"]["start"] == 4 and [v - 256 if v > 127 else v for v in zn["frame"]][35] == 0
    assert sn["cds_ends"]["disagreements"] == ce["disagreements"]  # still reported
    checks += 1
    # 52: sequence only: tD loses both marks, tC its start (CCC), tB its start; tA keeps both
    cut_full(pe, os.path.join(d, "pe_s"), 0, 0, set(), False, True, partial_ends="sequence")
    ss = json.load(open(os.path.join(d, "pe_s", "pe.w0.json")))
    assert ss["cds_ends"]["partial_5prime"] == ["tB", "tC", "tD"] and ss["cds_ends"]["partial_3prime"] == ["tD"], ss["cds_ends"]
    assert list(npz(os.path.join(d, "pe_s", "pe.w0.npz"))["bound"])[38] == 0
    assert ss["cds_ends"]["disagreements"][0]["decided"] == "sequence"
    checks += 1
    # 53: declared only: tB is complete (nothing declared), tC's start incomplete, tD complete
    cut_full(pe, os.path.join(d, "pe_d"), 0, 0, set(), False, True, partial_ends="declared")
    sd = json.load(open(os.path.join(d, "pe_d", "pe.w0.json")))
    assert sd["cds_ends"]["partial_5prime"] == ["tC"] and sd["cds_ends"]["partial_3prime"] == []
    checks += 1
    # 54: genetic code 6 (TAA and TAG read Gln): under the sequence rule tA and tC end on a non-stop
    cut_full(pe, os.path.join(d, "pe_6"), 0, 0, set(), False, True, partial_ends="sequence", genetic_code=6)
    s6 = json.load(open(os.path.join(d, "pe_6", "pe.w0.json")))
    assert s6["stop_codons"] == ["TGA"] and s6["cds_ends"]["partial_3prime"] == ["tA", "tC", "tD"], s6["cds_ends"]
    assert s6["cds_ends"]["sequence"]["stop_in_cds"] == 1 and s6["genetic_code"] == 6
    checks += 1
    # 55: a stop-excluded convention (CDS ends just before its stop) reads as complete via the next codon,
    #     a CDS reaching past the window is outside_window and left complete, and --stop-codons overrides
    with open(pe + ".annotation.json", "w") as fh:
        json.dump({"transcripts": [
            {"id": "tA", "strand": "+", "start": 2, "end": 11, "exons": [[2, 11]], "cds": [[2, 11]]},
            {"id": "tE", "strand": "+", "start": 44, "end": 60, "exons": [[44, 60]], "cds": [[44, 60]]},
        ]}, fh)
    cut_full(pe, os.path.join(d, "pe_x"), 0, 0, set(), False, True, partial_ends="both", stop_codons={"TAG"})
    sx = json.load(open(os.path.join(d, "pe_x", "pe.w0.json")))
    assert sx["cds_ends"]["sequence"]["stop_after_cds"] == 1 and sx["cds_ends"]["sequence"]["stop_outside_window"] == 1 \
        and sx["cds_ends"]["partial_3prime"] == [] and sx["cds_ends"]["partial_5prime"] == [] \
        and sx["stop_codons"] == ["TAG"] and sx["genetic_code"] is None, sx["cds_ends"]
    assert list(npz(os.path.join(d, "pe_x", "pe.w0.npz"))["bound"])[10] == 2  # last CDS base still marks the end
    checks += 1
    # 56: overlapping blocks: the first block in file order to cover a position wins,
    #     and the sidecar counts what the losing block carried there.  Reference
    #     ACGTTGCA; block 1 covers 1-4 (spA), block 2 covers 3-6 (spA with a gap, spB).
    ov = os.path.join(d, "ov")
    with open(ov + ".fa", "w") as fh:
        fh.write(">ov.chr1:0-8\nACGTTGCA\n")
    with open(ov + ".nh", "w") as fh:
        fh.write("((ov:0.1,spA:0.2):0.3,spB:0.5,spC:0.9);\n")
    with open(ov + ".annotation.json", "w") as fh:
        json.dump({"transcripts": [{"id": "t1", "strand": "+", "start": 1, "end": 7,
                                    "exons": [[1, 7]], "cds": [[1, 7]]}]}, fh)
    with open(ov + ".manifest.json", "w") as fh:
        json.dump({"assembly": "ov", "source": "ucsc", "alignment_track": "cactusToy",
                   "window": {"chrom": "chr1", "start": 0, "end": 8}}, fh)
    with open(ov + ".maf", "w") as fh:
        fh.write("##maf version=1\n"
                 "a score=0\n"
                 "s ov.chr1 0 2 - 8 TG\n"      # reference on the - strand: forward start 8-0-2 = 6, CA
                 "s spC.c   0 2 - 8 TG\n\n"
                 "a score=0\n"
                 "s ov.chr1 1 3 + 8 CGT\n"
                 "s spA.c   0 3 + 8 CGT\n\n"
                 "a score=0\n"
                 "s ov.chr1 3 3 + 8 TTG\n"
                 "s spA.c   0 2 + 8 T-G\n"
                 "s spB.c   0 3 + 8 TTG\n\n")
    cut(ov, os.path.join(d, "ov_out"), 0, 0, set(), False, True)
    so = json.load(open(os.path.join(d, "ov_out", "ov.w0.json")))
    zo = npz(os.path.join(d, "ov_out", "ov.w0.npz"))
    assert so["informants"] == ["spA", "spB", "spC"], so["informants"]
    bs = so["block_selection"]
    assert bs["overlapping_blocks"] == 1 and bs["overlap_positions"] == 1, bs
    # at position 3 the losing block carried spA=T (winner had T: kept) and spB=T (winner had no spB: lost)
    assert bs["overlap_informant_bases_discarded"] == 2 and bs["overlap_informant_bases_lost"] == 1 \
        and bs["overlap_lost_by_informant"] == {"spB": 1}, bs
    inf = list(zo["inf"])
    U, G = INF_UNALIGNED, INF_GAP
    assert inf[0:8] == [U, 1, 2, 3, G, 2, U, U], inf[0:8]      # spA: block 1 wins at 3, block 2 gives gap at 4 and G at 5
    assert inf[8:16] == [U, U, U, U, 3, 2, U, U], inf[8:16]    # spB: position 3 lost to block 1, 4-5 from block 2
    checks += 1
    # 57: a block whose reference row is on the - strand is reverse-complemented into
    #     forward coordinates (srcSize - start - size) rather than skipped: spC aligns CA at 6-8
    assert bs["reference_minus_strand_blocks_flipped"] == 1, bs
    assert inf[16:24] == [U, U, U, U, U, U, 1, 0], inf[16:24]  # spC: C at 6, A at 7
    checks += 1
    # 58: duplicate rows: two spA copies in one block; identity keeps the copy matching the
    #     reference (the second), first keeps the first, and the sidecar counts the loser
    with open(ov + ".maf", "w") as fh:
        fh.write("##maf version=1\n"
                 "a score=0\n"
                 "s ov.chr1 1 4 + 8 CGTT\n"
                 "s spA.c   0 4 + 8 CAAT\n"
                 "s spA.d   0 3 + 8 CGT-\n"
                 "s spB.c   0 4 + 8 CGTT\n\n")
    cut(ov, os.path.join(d, "dup_i"), 0, 0, set(), False, True)
    si = json.load(open(os.path.join(d, "dup_i", "ov.w0.json")))
    bi = si["block_selection"]
    assert bi["duplicate_row_policy"] == "identity" and bi["duplicate_row_blocks"] == 1 \
        and bi["duplicate_rows_discarded"] == 1 and bi["duplicate_rows_discarded_bases_in_window"] == 4 \
        and bi["duplicates_by_informant"] == {"spA": {"rows": 1, "blocks": 1, "max_copies": 2, "bases": 4}} \
        and bi["duplicate_informants"] == 1 and bi["duplicate_max_copies"] == 2 \
        and bi["reference_duplicate_rows"] == 0, bi
    # kept bases: spA keeps CGT (3), spB keeps CGTT (4): the unit the discarded count compares against
    assert bi["kept_informant_bases_in_window"] == 7, bi
    assert set(bi["units"]) >= {k for k in bi if k != "units"}, "every counter has a stated unit"
    assert list(npz(os.path.join(d, "dup_i", "ov.w0.npz"))["inf"])[0:8] == [U, 1, 2, 3, G, U, U, U]
    cut(ov, os.path.join(d, "dup_f"), 0, 0, set(), False, True, duplicate_rows="first")
    sf = json.load(open(os.path.join(d, "dup_f", "ov.w0.json")))
    assert sf["block_selection"]["duplicate_rows_discarded_bases_in_window"] == 3
    assert sf["block_selection"]["duplicates_by_informant"]["spA"]["bases"] == 3 \
        and sf["block_selection"]["duplicates_by_informant"]["spA"]["rows"] == 1
    assert sf["block_selection"]["kept_informant_bases_in_window"] == 8
    checks += 1
    # 59: the same informant duplicated in two blocks: rows and blocks are counted separately,
    #     max_copies is per block, and a third row in one block counts two discarded rows
    with open(ov + ".maf", "w") as fh:
        fh.write("##maf version=1\n"
                 "a score=0\n"
                 "s ov.chr1 0 3 + 8 ACG\n"
                 "s spA.c   0 3 + 8 ACG\n"
                 "s spA.d   0 3 + 8 AAA\n"
                 "s spA.e   0 3 + 8 CCC\n\n"
                 "a score=0\n"
                 "s ov.chr1 3 3 + 8 TTG\n"
                 "s spA.c   0 3 + 8 TTG\n"
                 "s spA.d   0 2 + 8 T-G\n\n")
    cut(ov, os.path.join(d, "dup_2"), 0, 0, set(), False, True)
    b2 = json.load(open(os.path.join(d, "dup_2", "ov.w0.json")))["block_selection"]
    assert b2["duplicates_by_informant"] == {"spA": {"rows": 3, "blocks": 2, "max_copies": 3, "bases": 8}}, b2
    assert b2["duplicate_row_blocks"] == 2 and b2["duplicate_rows_discarded"] == 3 \
        and b2["duplicate_max_copies"] == 3 and b2["duplicate_informants"] == 1 \
        and b2["kept_informant_bases_in_window"] == 6, b2
    assert list(npz(os.path.join(d, "dup_f", "ov.w0.npz"))["inf"])[0:8] == [U, 1, 0, 0, 3, U, U, U]
    checks += 1
    # 60-62: feature_lengths.  Two isoforms share the 20-base intron; the
    # 45-base one is the second isoform's alone; the 100-base intron reaches
    # past the window edge and is clipped, as is that transcript's span.
    fl = feature_lengths([
        {"strand": "+", "exons": [[0, 10], [30, 40], [85, 100]], "cds": [[5, 10], [30, 40], [85, 90]]},
        {"strand": "+", "exons": [[0, 10], [30, 40], [85, 100], [200, 210]], "cds": [[5, 10], [30, 40]]},
        {"strand": "-", "exons": [[50, 60]], "cds": []},
    ], 0, 150)
    assert fl["introns"] == {"n": 2, "min": 20, "median": 32.5, "max": 45, "below_30": 1, "below_50": 2,
                             "clipped": 1, "below_min_intron": 0} and fl["min_intron"] == 20, fl["introns"]
    checks += 1
    assert fl["cds"] == {"n": 2, "min": 15, "median": 17.5, "max": 20, "below_60": 2, "outside_window": 0}, fl["cds"]
    checks += 1
    assert fl["span"] == {"n": 3, "min": 10, "median": 100, "max": 210, "below_81": 1, "outside_window": 1} \
        and set(fl["units"]) == {"lengths", "n", "below_F", "below_min_intron"}, fl["span"]
    checks += 1
    # 63: the toy window's sidecar carries the record for all and painted isoforms
    st = json.load(open(os.path.join(d, "out", "toy.w0.json")))["feature_lengths"]
    assert st["floors"] == {"intron": [30, 50], "cds": [60], "span": [81]} \
        and st["all_isoforms"]["introns"]["n"] == 1 and st["painted"]["introns"]["below_30"] == 1 \
        and st["painted"]["span"]["below_81"] == 2, st
    checks += 1
    # 64-65: the intron floor.  stalin's length sweep (relay note
    # 20260909T214050Z-stalin-0021): two exons with gaps of 1, 15, 16, 19, 20
    # and 30 bases on both strands; under MIN_INTRON 20 the first four are
    # short_gap with no donor or acceptor and one short_gap site, the last
    # two are introns with both marks.
    for strand in ("+", "-"):
        for g in (1, 15, 16, 19, 20, 30):
            tx = [{"id": "t", "strand": strand, "start": 100, "end": 400 + g,
                   "exons": [[100, 200], [200 + g, 400 + g]], "cds": [[100, 200], [200 + g, 400 + g]]}]
            lab, _, _, bnd = paint_labels(0, 500, tx)
            kinds = {k for k, _, _ in distinct_sites(tx)}
            gap = set(lab[200:200 + g])
            marks = {v for v in bnd[200:200 + g] if v}
            if g < 20:
                assert gap == {LABEL["short_gap"]} and not marks and kinds == {"start", "stop", "short_gap"}, (strand, g, gap, marks, kinds)
            else:
                assert gap == {LABEL["intron"]} and marks == {3, 4} \
                    and kinds == {"start", "stop", "donor", "acceptor", "cds_donor", "cds_acceptor"}, (strand, g, gap, marks, kinds)
            fl = feature_lengths(tx, 0, 500)
            assert fl["introns"]["n"] == 1 and fl["introns"]["below_min_intron"] == int(g < 20), (g, fl["introns"])
    checks += 1
    # 65: the floor is a policy: --min-intron 15 makes a 15-base gap an intron
    #     (Stentor's 15- and 16-base introns) and 0 turns the floor off for a 1-base gap
    tx = [{"id": "t", "strand": "+", "start": 100, "end": 415, "exons": [[100, 200], [215, 415]], "cds": [[100, 200], [215, 415]]}]
    lab, _, _, bnd = paint_labels(0, 500, tx, min_intron=15)
    assert set(lab[200:215]) == {LABEL["intron"]} and bnd[200] == 3 and bnd[214] == 4
    assert {k for k, _, _ in distinct_sites(tx, min_intron=15)} >= {"donor", "acceptor"}
    tx1 = [{"id": "t", "strand": "+", "start": 100, "end": 401, "exons": [[100, 200], [201, 401]], "cds": [[100, 200], [201, 401]]}]
    lab, _, _, bnd = paint_labels(0, 500, tx1, min_intron=0)
    assert lab[200] == LABEL["intron"] and bnd[200] == 3, (lab[200], bnd[200])  # a 1-base intron: donor and acceptor at one base, the donor is put first and wins
    checks += 1
    # 66: engels' overlapping motif windows (relay note 20260909T212826Z-engels-0021):
    #     exon ending in A, gap G, exon starting T reads GT and AG by borrowing a base
    #     from each exon; GTAG satisfies both windows on its own bases; GT, GTA and
    #     ATAC do not; the same on the minus strand through the reverse complement.
    assert motif_window("CA", "G", "TC") == (True, True)
    assert motif_window("CA", "G", "CC") == (True, True)   # A|G|C: GC donor, AG acceptor
    assert motif_window("CC", "G", "TC") == (False, False)  # no A before the gap
    assert motif_window("CA", "GTAG", "TC") == (True, False)
    assert motif_window("CA", "GT", "TC") == (False, False)
    assert motif_window("CA", "GTA", "TC") == (False, False)
    assert motif_window("CA", "ATAC", "TC") == (False, False)
    # 0.10: the two windows apart, and the class names which side failed (engels,
    #       note 20260909T222452Z-engels-0022): a GT donor with no AG, an AG with no
    #       donor, neither, and an N in a window leaves it unresolved rather than failed
    mw = lambda l, g, r: (motif_windows(l, g, r)["donor"], motif_windows(l, g, r)["acceptor"], motif_windows(l, g, r)["class"])
    assert mw("CA", "G", "TC") == (True, True, "borrows_exon_base")
    assert mw("CA", "GTAG", "TC") == (True, True, "exact")
    assert mw("CA", "GTAA", "CC") == (True, False, "donor_only")
    assert mw("CA", "AAAG", "CC") == (False, True, "acceptor_only")
    assert mw("CA", "ATAC", "CC") == (False, False, "neither")
    assert mw("CN", "G", "TC") == (True, None, "ambiguous")      # N|G|T: donor reads GT, the acceptor window is N,G
    assert mw("CA", "N", "TC") == (None, None, "ambiguous")
    assert mw("CA", "GTAN", "CC") == (True, None, "ambiguous")
    assert mw("CA", "ATAN", "CC") == (False, None, "ambiguous")  # one side undecided, the other failed: still ambiguous, not neither
    assert mw("CN", "GTAA", "CC") == (True, False, "donor_only")  # the N sits outside both windows
    assert motif_window("CN", "G", "TC") == (False, False)       # the conjunction counts unresolved as not passing
    assert mw("ca", "g", "tc") == (True, True, "borrows_exon_base")
    seq = "CCCCCCCCAGTCCCCCCCCCCCCC"  # + strand: exon ..CA | G | TC.. with the gap at 9
    tx = [{"id": "p", "strand": "+", "start": 0, "end": 24, "exons": [[0, 9], [10, 24]], "cds": [[0, 9], [10, 24]]}]
    sg = short_gaps(tx, seq, 0, 24)
    assert sg["n"] == 1 and sg["in_cds"] == 1 and sg["by_length"] == {"1": 1} and sg["motif_window"] == 1 \
        and sg["motif_borrows_exon_base"] == 1 and sg["motif_exact"] == 0, sg
    assert sg["motif_by_class"] == {"exact": 0, "borrows_exon_base": 1, "donor_only": 0, "acceptor_only": 0,
                                    "neither": 0, "ambiguous": 0} and sg["motif_unresolved_windows"] == 0, sg
    assert (sg["gaps"][0]["motif_donor"], sg["gaps"][0]["motif_acceptor"], sg["gaps"][0]["motif_class"]) == (True, True, "borrows_exon_base")
    g = sg["gaps"][0]
    assert (g["left"], g["gap"], g["right"], g["length_mod_3"], g["transcripts"]) == ("CA", "G", "TC", 1, ["p"]), g
    rc = seq.translate(COMP)[::-1]  # minus strand: the same gap at 24-10=14
    txm = [{"id": "m", "strand": "-", "start": 0, "end": 24, "exons": [[0, 14], [15, 24]], "cds": [[0, 14], [15, 24]]}]
    gm = short_gaps(txm, rc, 0, 24)["gaps"][0]
    assert (gm["left"], gm["gap"], gm["right"], gm["motif_window"], gm["start"]) == ("CA", "G", "TC", True, 14), gm
    # a gap on the UTR side of the CDS is not in_cds; a 4-base GTAG gap is motif_exact
    seq4 = "CCCCCCCCAGTAGTCCCCCCCCCCCC"
    tx4 = [{"id": "u", "strand": "+", "start": 0, "end": 26, "exons": [[0, 9], [13, 26]], "cds": [[15, 26]]}]
    s4 = short_gaps(tx4, seq4, 0, 26)
    assert s4["in_cds"] == 0 and s4["motif_exact"] == 1 and s4["motif_window"] == 1 \
        and s4["motif_borrows_exon_base"] == 0 and s4["gaps"][0]["gap"] == "GTAG", s4
    assert s4["gaps"][0]["motif_class"] == "exact" and s4["motif_by_class"]["exact"] == 1
    # a gap whose flanks reach past the window edge is listed without sequence
    s5 = short_gaps(tx, seq, 8, 24)
    assert s5["outside_window"] == 1 and s5["gaps"][0]["motif_window"] is None and s5["gaps"][0]["motif_class"] is None, s5
    assert s5["motif_by_class"] == {c: 0 for c in MOTIF_CLASSES} and s5["motif_unresolved_windows"] == 0, s5
    checks += 1
    # 69: the classes on real-shaped gaps through short_gaps, both strands: a
    #     donor-only 4-base gap, an acceptor-only one, the Ty1 site TT|A|GG (neither),
    #     and a minus-strand gap over an N, whose class is ambiguous and whose
    #     unresolved window is counted once
    seq9 = "CCCCCCCCAGTAACCCCCCCCAAAAGCCCCCCTTAGGCCCCCCCNGTCCCC"
    tx9 = [{"id": "q", "strand": "+", "start": 0, "end": 51,
            "exons": [[0, 9], [13, 22], [26, 34], [35, 44], [45, 51]],
            "cds": [[0, 9], [13, 22], [26, 34], [35, 44], [45, 51]]}]
    s9 = short_gaps(tx9, seq9, 0, 51)
    assert [g["gap"] for g in s9["gaps"]] == ["GTAA", "AAAG", "A", "N"], s9["gaps"]
    assert [g["motif_class"] for g in s9["gaps"]] == ["donor_only", "acceptor_only", "neither", "ambiguous"], s9["gaps"]
    assert s9["gaps"][2]["left"] == "TT" and s9["gaps"][2]["right"] == "GG"
    assert s9["motif_by_class"] == {"exact": 0, "borrows_exon_base": 0, "donor_only": 1, "acceptor_only": 1,
                                    "neither": 1, "ambiguous": 1} and s9["motif_window"] == 0, s9
    assert s9["motif_unresolved_windows"] == 2, s9["motif_unresolved_windows"]  # the N sits in both windows of a 1-base gap
    rc9 = seq9.translate(COMP)[::-1]
    txm9 = [{"id": "qm", "strand": "-", "start": 0, "end": 51,
             "exons": [[51 - b, 51 - a] for a, b in reversed(tx9[0]["exons"])],
             "cds": [[51 - b, 51 - a] for a, b in reversed(tx9[0]["cds"])]}]
    sm9 = short_gaps(txm9, rc9, 0, 51)
    assert sorted(g["motif_class"] for g in sm9["gaps"]) == sorted(g["motif_class"] for g in s9["gaps"]), sm9["gaps"]
    assert sm9["motif_by_class"] == s9["motif_by_class"] and sm9["motif_unresolved_windows"] == 2
    assert {g["gap"] for g in sm9["gaps"]} == {"GTAA", "AAAG", "A", "N"}, sm9["gaps"]
    checks += 1
    # 70: engels' three fixtures (note 20260909T232517Z-engels-0023): an N
    #     inside the gap but outside both two-base windows leaves the class
    #     exact and the windows resolved, so only gap_unresolved_bases sees
    #     it; an unresolved window beside a failed one is ambiguous on
    #     either side; the counters agree through short_gaps on both strands
    assert mw("A", "GTNAG", "C") == (True, True, "exact")
    assert mw("A", "AANA", "C") == (False, None, "ambiguous")
    assert mw("A", "NTAA", "C") == (None, False, "ambiguous")
    assert motif_fields("A", "GTNAG", "C")["gap_unresolved_bases"] == 1
    assert motif_fields("A", "GTAG", "C")["gap_unresolved_bases"] == 0
    assert motif_fields("CA", "N", "TC")["gap_unresolved_bases"] == 1
    assert motif_fields("a", "gtnag", "c")["gap_unresolved_bases"] == 1
    seq10 = "CCCCCCCCA" + "GTNAG" + "CCCCCCCCA" + "AANA" + "CCCCCCCA" + "NTAA" + "CCCCCCC"  # 46 bases
    tx10 = [{"id": "r", "strand": "+", "start": 0, "end": 46,
             "exons": [[0, 9], [14, 23], [27, 35], [39, 46]],
             "cds": [[0, 9], [14, 23], [27, 35], [39, 46]]}]
    s10 = short_gaps(tx10, seq10, 0, 46)
    assert [g["gap"] for g in s10["gaps"]] == ["GTNAG", "AANA", "NTAA"], s10["gaps"]
    assert [g["motif_class"] for g in s10["gaps"]] == ["exact", "ambiguous", "ambiguous"], s10["gaps"]
    assert [g["gap_unresolved_bases"] for g in s10["gaps"]] == [1, 1, 1], s10["gaps"]
    assert s10["motif_unresolved_windows"] == 2 and s10["gaps_with_unresolved_bases"] == 3, s10
    assert s10["motif_by_class"]["neither"] == 0 and s10["motif_exact"] == 1, s10
    rc10 = seq10.translate(COMP)[::-1]
    txm10 = [{"id": "rm", "strand": "-", "start": 0, "end": 46,
              "exons": [[46 - b, 46 - a] for a, b in reversed(tx10[0]["exons"])],
              "cds": [[46 - b, 46 - a] for a, b in reversed(tx10[0]["cds"])]}]
    sm10 = short_gaps(txm10, rc10, 0, 46)
    assert {g["gap"] for g in sm10["gaps"]} == {"GTNAG", "AANA", "NTAA"}, sm10["gaps"]
    assert sm10["motif_unresolved_windows"] == 2 and sm10["gaps_with_unresolved_bases"] == 3, sm10
    assert s9["gaps_with_unresolved_bases"] == 1 and sg["gaps_with_unresolved_bases"] == 0
    assert s5["gaps_with_unresolved_bases"] == 0, s5
    checks += 1
    # 67: precedence: a CDS in another isoform over a short gap paints CDS; a short gap
    #     over another transcript's intron paints short_gap
    two = [{"id": "a", "strand": "+", "start": 100, "end": 401, "exons": [[100, 200], [201, 401]], "cds": [[100, 200], [201, 401]]},
           {"id": "b", "strand": "+", "start": 100, "end": 401, "exons": [[100, 401]], "cds": [[100, 401]]}]
    lab, _, _, _ = paint_labels(0, 500, two)
    assert lab[200] == LABEL["cds"], lab[200]
    two[1] = {"id": "c", "strand": "-", "start": 0, "end": 500, "exons": [[0, 10], [490, 500]], "cds": []}
    lab, _, _, _ = paint_labels(0, 500, two)
    assert lab[200] == LABEL["short_gap"] and lab[50] == LABEL["intron"], (lab[200], lab[50])
    checks += 1
    # 68: the toy window under the default floor: its 4-base gap is short_gap in the
    #     sidecar counts and in both orientations, with no donor or acceptor mark and
    #     the short_gaps record carrying the flanks; under min_intron 0 nothing changes
    cut(stem, os.path.join(d, "floor"), 0, 0, {"spC"}, True, True, min_intron=20)
    sf = json.load(open(os.path.join(d, "floor", "toy.w0.json")))
    zf = npz(os.path.join(d, "floor", "toy.w0.npz"))
    assert sf["min_intron"] == 20 and sf["label_counts"]["short_gap"] == 4 and sf["label_counts"]["intron"] == 0, sf["label_counts"]
    assert list(zf["label"])[8:12] == [4, 4, 4, 4] and 3 not in zf["bound"] and 4 not in zf["bound"], list(zf["label"])
    assert sf["distinct_sites"]["painted"]["short_gap"] == 1 and sf["distinct_sites"]["painted"]["donor"] == 0
    assert sf["feature_lengths"]["painted"]["introns"]["below_min_intron"] == 1
    sgf = sf["short_gaps"]["painted"]
    assert sgf["n"] == 1 and sgf["gaps"][0]["gap"] == ref[8:12] and sgf["gaps"][0]["in_cds"] and sgf["gaps"][0]["length"] == 4, sgf
    zr = npz(os.path.join(d, "floor", "toy.w0.rc.npz"))
    assert list(zr["label"])[8:12] == [4, 4, 4, 4], list(zr["label"])
    assert side["label_counts"].get("short_gap", 0) == 0 and side["min_intron"] == 0
    checks += 1
    print(f"self-test passed ({checks} checks)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Cut fetch_window.py output into fixed-length training examples.")
    ap.add_argument("--stem", action="append", default=[], help="fetch_window output stem (path without extension); repeatable")
    ap.add_argument("--out", default=None, help="output directory")
    ap.add_argument("--length", type=int, default=4096, help="example length L; 0 = whole window")
    ap.add_argument("--stride", type=int, default=2048)
    ap.add_argument("--drop-species", default="", help="comma-separated informant names to remove (leakage rule)")
    ap.add_argument("--both-strands", action="store_true", help="also write the reverse-complement example")
    ap.add_argument("--transcript-types", default="benchmark",
                    help="which transcripts paint labels: 'benchmark' (default; drops pseudogenes and Ig/TCR "
                         "segments as docs/benchmark.md section 4 does), 'all', or a comma-separated list of types")
    ap.add_argument("--isoforms", default="union", choices=["union", "longest-cds", "representative"],
                    help="which isoforms of a locus paint labels: 'union' (default; every transcript), "
                         "'longest-cds' (one per locus, most CDS bases), or 'representative' (the ids in "
                         "--representatives, e.g. MANE Select; a locus with none listed falls back to longest-cds "
                         "and the sidecar names it)")
    ap.add_argument("--representatives", default=None,
                    help="file of transcript ids, one per line (first column of a TSV), for --isoforms representative")
    ap.add_argument("--drop-ancestors", action="store_true", help="remove every inferred ancestral row (Ensembl EPO)")
    ap.add_argument("--reference-anchored", action="store_true",
                    help="declare a track the built-in table does not know to be reference-anchored, so dropped "
                         "rows count as exact removal; recorded in the sidecar as the operator's claim")
    ap.add_argument("--duplicate-rows", default="identity", choices=["identity", "first"],
                    help="when one species has several rows in one block (Cactus exports every copy of a "
                         "duplicated region): keep the copy with most bases identical to the reference, "
                         "or the first row; the sidecar counts what was discarded (default identity)")
    ap.add_argument("--partial-ends", default="both", choices=list(PARTIAL_END_MODES),
                    help="how an incomplete CDS end is recognised: 'both' (default; the annotation's "
                         "declaration where it makes one, else the reference sequence), 'declared', "
                         "'sequence', or 'none' (mark every end as a codon, the pre-0.5 behaviour)")
    ap.add_argument("--genetic-code", type=int, default=1,
                    help="NCBI genetic code table for the stop-codon check (1 standard, 6 ciliate, ...)")
    ap.add_argument("--stop-codons", default=None,
                    help="comma-separated stop codons, overriding --genetic-code")
    ap.add_argument("--min-intron", type=int, default=MIN_INTRON,
                    help="an exon gap shorter than this is painted short_gap, not intron, and gets no donor "
                         "or acceptor mark (default 20, the benchmark scorer's MIN_INTRON; 0 turns the floor off)")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if not args.stem or not args.out:
        ap.error("--stem and --out are required")
    drop = {s.strip() for s in args.drop_species.split(",") if s.strip()}
    stop_codons = {c.strip().upper() for c in args.stop_codons.split(",") if c.strip()} if args.stop_codons else None
    reps = read_id_list(args.representatives) if args.representatives else None
    if args.isoforms == "representative" and not reps:
        ap.error("--isoforms representative needs --representatives FILE with at least one id")
    total = 0
    for stem in args.stem:
        total += len(cut(stem, args.out, args.length, args.stride, drop, args.both_strands, args.quiet,
                         args.transcript_types, args.drop_ancestors, args.reference_anchored,
                         args.isoforms, reps, args.partial_ends, args.genetic_code, stop_codons,
                         args.duplicate_rows, args.min_intron))
    print(f"wrote {total} examples to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
