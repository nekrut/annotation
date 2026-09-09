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

The sidecar JSON carries the coordinates, the informant names in row order
(the tree's depth-first leaf order, so the same track always yields the
same rows in the same order), which rows are inferred ancestors, which
species were dropped for leakage (``--drop-species``), the transcripts, the
source manifest's SHA-256 and the tool version.

Conventions worth stating once:

* Examples are cut on the fetched window with length ``L`` and stride ``S``
  (defaults 4096 / 2048); the last example is padded, never dropped, and
  ``mask`` says where the padding starts.  ``--length 0`` emits the whole
  fetched window as one unpadded example.
* Reverse-strand genes are *not* flipped here.  ``--both-strands`` writes a
  second example ``.rc.npz`` with every array reverse-complemented (labels
  and boundaries recomputed on the flipped coordinates), so a model can be
  trained orientation-symmetric; the default keeps one orientation and
  leaves augmentation to the loader.
* Stop codons are inside the CDS, as in RefSeq genePred and GFF3 from NCBI;
  ``bound == 2`` marks the last base of that codon.
* ``--drop-species`` implements the benchmark leakage rule for training
  alignments (docs/benchmark.md 3.2, docs/data-sources.md section 7): a
  held-out species' row is removed before the tensors are built.  For multiz
  this is exact; for Cactus and EPO it removes the sequence but not its
  influence on the alignment, which the sidecar records as
  ``"dropped_rows_only": true``.
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
import struct
import sys
import zipfile

TOOL_VERSION = "0.1"
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
INF_GAP, INF_UNALIGNED, INF_OTHER = 4, 5, 6
LABEL = {"intergenic": 0, "intron": 1, "utr": 2, "cds": 3}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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

def newick_leaves_and_distances(text: str, ref: str) -> tuple[list[str], dict[str, float]]:
    """Depth-first leaf order and patristic distance from ``ref`` to each leaf."""
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
    leaves = [name[i] for i in range(len(name)) if not children[i] and name[i]]
    idx = {n: i for i, n in enumerate(name) if n}
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
                if not children[c] and name[c]:
                    dist[name[c]] = dc
                down(c, dc)

        for node, d in up.items():
            down(node, d)
    return leaves, dist


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
    """Return inf [K][n], ins [K][n] bytearrays and the count of overlapping blocks."""
    n = end - start
    K = len(rows)
    row_of = {r: i for i, r in enumerate(rows)}
    inf = [bytearray([INF_UNALIGNED]) * n for _ in range(K)]
    ins = [bytearray(n) for _ in range(K)]
    seen = bytearray(n)
    overlaps = 0
    for block in blocks:
        rref = next((r for r in block if r[1].split(".")[0] == ref), None)
        if rref is None or rref[4] != "+":
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
    return inf, ins, overlaps


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


def cut(stem: str, out_dir: str, L: int, S: int, drop: set[str], both: bool, quiet: bool) -> list[str]:
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
    # rows: tree leaf order when there is one tree; else names in order of appearance
    leaves, dist = ([], {})
    if len(trees) == 1:
        leaves, dist = newick_leaves_and_distances(trees[0], ref_prefix)
    sources: list[str] = []
    for b in blocks:
        for r in b:
            nm = r[1].split(".")[0]
            if nm not in sources:
                sources.append(nm)
    rows = [n for n in leaves if n != ref_prefix and n.split(".")[0] != ref_prefix]
    rows += [s for s in sources if s not in rows and s != ref_prefix]
    if len(trees) > 1:
        for t in trees:
            _, d = newick_leaves_and_distances(t, ref_prefix)
            for k, v in d.items():
                dist.setdefault(k, v)
    dropped = [r for r in rows if r in drop]
    rows = [r for r in rows if r not in drop]
    K = len(rows)
    inf, ins, overlaps = paint_informants(blocks, ref_prefix, start, end, rows)
    with open(stem + ".annotation.json") as fh:
        ann = json.load(fh)
    transcripts = ann.get("transcripts", ann if isinstance(ann, list) else [])
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
                "alignment_track": manifest.get("alignment_track") or manifest.get("species_set"),
                "informants": rows,
                "ancestral": [r for r in rows if r.lower().startswith("ancestor")],
                "dropped_species": dropped,
                "dropped_rows_only": bool(dropped) and (manifest.get("source") == "ensembl"
                                                         or "cactus" in str(manifest.get("alignment_track", "")).lower()),
                "overlapping_blocks": overlaps,
                "transcripts": [t.get("id") for t in transcripts],
                "label_counts": {c: sum(1 for v in e["label"][:b - a] if v == i) for c, i in LABEL.items()},
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
                                    "exons": [[0, 2]], "cds": []}]}, fh)
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
    print("self-test passed (12 checks)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Cut fetch_window.py output into fixed-length training examples.")
    ap.add_argument("--stem", action="append", default=[], help="fetch_window output stem (path without extension); repeatable")
    ap.add_argument("--out", default=None, help="output directory")
    ap.add_argument("--length", type=int, default=4096, help="example length L; 0 = whole window")
    ap.add_argument("--stride", type=int, default=2048)
    ap.add_argument("--drop-species", default="", help="comma-separated informant names to remove (leakage rule)")
    ap.add_argument("--both-strands", action="store_true", help="also write the reverse-complement example")
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
        total += len(cut(stem, args.out, args.length, args.stride, drop, args.both_strands, args.quiet))
    print(f"wrote {total} examples to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
