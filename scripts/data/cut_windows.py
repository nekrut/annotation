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
                           (union over transcripts, CDS > UTR > intron)
  ``frame``  int8  [L]     codon position 0,1,2 of a CDS base counted from the
                           start codon on the transcript's strand; -1 elsewhere
  ``strand`` int8  [L]     +1 / -1 strand of the transcript that labelled the
                           base, 0 for intergenic
  ``bound``  uint8 [L]     0 none, 1 first base of start codon, 2 last base of
                           stop codon, 3 donor (first intron base),
                           4 acceptor (last intron base); all placed at their
                           + strand coordinate, strand given by ``strand``
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
excluded by type, the source manifest's SHA-256 and the tool version.

Conventions worth stating once:

* Examples are cut on the fetched window with length ``L`` and stride ``S``
  (defaults 4096 / 2048); the last example is padded, never dropped, and
  ``mask`` says where the padding starts.  ``--length 0`` emits the whole
  fetched window as one unpadded example.
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
import struct
import sys
import zipfile

TOOL_VERSION = "0.1"
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
INF_GAP, INF_UNALIGNED, INF_OTHER = 4, 5, 6
LABEL = {"intergenic": 0, "intron": 1, "utr": 2, "cds": 3}
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


# --------------------------------------------------------------- npy/npz

def npy_bytes(data: bytes, dtype: str, shape: tuple[int, ...]) -> bytes:
    header = "{'descr': '%s', 'fortran_order': False, 'shape': %s, }" % (
        dtype, "(%s)" % "".join(f"{n}, " for n in shape) if shape else "()")
    pad = 64 - ((10 + len(header) + 1) % 64)
    header = header + " " * pad + "\n"
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header.encode("latin1") + data


def write_npz(path: str, arrays: dict[str, tuple[bytes, str, tuple[int, ...]]]) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name, (data, dtype, shape) in arrays.items():
            z.writestr(name + ".npy", npy_bytes(data, dtype, shape))


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
    return m.group(1) if m else label.split(".")[0]


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


def paint_informants(blocks, ref: str, start: int, end: int, rows: list[str]):
    """Return inf [K][n], ins [K][n] bytearrays, the count of overlapping blocks
    and the count of blocks skipped because their reference row is on the - strand."""
    n = end - start
    K = len(rows)
    row_of = {r: i for i, r in enumerate(rows)}
    inf = [bytearray([INF_UNALIGNED]) * n for _ in range(K)]
    ins = [bytearray(n) for _ in range(K)]
    seen = bytearray(n)
    overlaps = 0
    minus = 0
    for block in blocks:
        rref = next((r for r in block if r[1].split(".")[0] == ref), None)
        if rref is None:
            continue
        if rref[4] != "+":
            minus += 1
            continue
        rseq, p = rref[6], int(rref[2])
        col_ref = []
        for ch in rseq:
            if ch == "-":
                col_ref.append(-1)
            else:
                col_ref.append(p); p += 1
        # first block to cover a position wins (Cactus can overlap on the reference)
        cover = [q for q in col_ref if start <= q < end]
        if cover and any(seen[q - start] for q in cover):
            overlaps += 1
        for r in block:
            if r is rref:
                continue
            k = row_of.get(r[1].split(".")[0])
            if k is None:
                continue
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
                if seen[q - start]:
                    continue
                ch = seq[col].upper()
                inf[k][q - start] = BASE.get(ch, INF_GAP if ch == "-" else INF_OTHER)
        for q in cover:
            seen[q - start] = 1
    return inf, ins, overlaps, minus


# ------------------------------------------------------------- annotation

def paint_labels(start: int, end: int, transcripts: list[dict]):
    n = end - start
    label = bytearray(n)
    frame = bytearray(b"\xff" * n)      # -1 as int8
    strand = bytearray(n)                # 0, 1 (+), 255 (-) as int8
    bound = bytearray(n)

    def put(arr, pos, val, force=False):
        if start <= pos < end and (force or arr[pos - start] == 0):
            arr[pos - start] = val

    for t in transcripts:
        s = 1 if t.get("strand", "+") == "+" else 255
        for i in range(max(t["start"], start), min(t["end"], end)):
            if label[i - start] < LABEL["intron"]:
                label[i - start] = LABEL["intron"]
            if strand[i - start] == 0:
                strand[i - start] = s
        for a, b in t.get("exons", []):
            for i in range(max(a, start), min(b, end)):
                if label[i - start] < LABEL["utr"]:
                    label[i - start] = LABEL["utr"]
        cds = sorted(t.get("cds", []))
        for a, b in cds:
            for i in range(max(a, start), min(b, end)):
                label[i - start] = LABEL["cds"]
        if cds:
            k = 0
            order = cds if s == 1 else [(b, a) for a, b in reversed(cds)]
            for a, b in order:
                rng = range(a, b) if s == 1 else range(a - 1, b - 1, -1)
                for i in rng:
                    if start <= i < end and frame[i - start] == 255:
                        frame[i - start] = k % 3
                    k += 1
            if s == 1:
                put(bound, cds[0][0], 1); put(bound, cds[-1][1] - 1, 2)
            else:
                put(bound, cds[-1][1] - 1, 1); put(bound, cds[0][0], 2)
        ex = sorted(t.get("exons", []))
        for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
            if b1 >= a2:
                continue
            if s == 1:
                put(bound, b1, 3); put(bound, a2 - 1, 4)
            else:
                put(bound, a2 - 1, 3); put(bound, b1, 4)
    return label, frame, strand, bound


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
        reference_anchored: bool = False) -> list[str]:
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
            nm = r[1].split(".")[0]
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
    inf, ins, overlaps, minus_blocks = paint_informants(blocks, ref_prefix, start, end, rows)
    with open(stem + ".annotation.json") as fh:
        ann = json.load(fh)
    transcripts = ann.get("transcripts", ann if isinstance(ann, list) else [])
    transcripts, excluded = select_transcripts(transcripts, transcript_types)
    label, frame, strand, bound = paint_labels(start, end, transcripts)
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
                "overlapping_blocks": overlaps,
                "reference_minus_strand_blocks": minus_blocks,
                "transcripts": [t.get("id") for t in transcripts],
                "transcript_types": transcript_types,
                "transcripts_excluded": excluded,
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
    assert fwd["label_counts"] == {"intergenic": 1, "intron": 4, "utr": 6, "cds": 9}, fwd["label_counts"]
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
    ap.add_argument("--drop-ancestors", action="store_true", help="remove every inferred ancestral row (Ensembl EPO)")
    ap.add_argument("--reference-anchored", action="store_true",
                    help="declare a track the built-in table does not know to be reference-anchored, so dropped "
                         "rows count as exact removal; recorded in the sidecar as the operator's claim")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    if not args.stem or not args.out:
        ap.error("--stem and --out are required")
    drop = {s.strip() for s in args.drop_species.split(",") if s.strip()}
    total = 0
    for stem in args.stem:
        total += len(cut(stem, args.out, args.length, args.stride, drop, args.both_strands, args.quiet,
                         args.transcript_types, args.drop_ancestors, args.reference_anchored))
    print(f"wrote {total} examples to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
