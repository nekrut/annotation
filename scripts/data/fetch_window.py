#!/usr/bin/env python3
"""Fetch one locus window with its multiple alignment, tree, annotation and
conservation scores, in a form a model can consume.

Standard library only (Python 3.11).  Nothing is cached in the repository:
the output directory is wherever you point ``--out``.

Two backends:

``--source ucsc`` (default for UCSC assembly names such as ``hg38``, ``dm6``,
and for GenArk assembly-hub accessions such as ``GCF_000002765.6``)

  * alignment: the UCSC JSON API (``api.genome.ucsc.edu``).  For *bigMaf*
    tracks (hg38 ``multiz470way``, ``cactus241wayBM``, ``cactus447way``) the
    API returns the MAF blocks directly.  For the older *wigMaf* tracks
    (hg38 ``multiz100way``, dm6 ``multiz124way``, mm39 ``multiz35way``, ...)
    the API returns only an index (file id + byte offset); the script then
    issues HTTP Range requests against the uncompressed per-chromosome MAF
    that hgdownload serves under ``/gbdb/<db>/<track>/maf/<chrom>.maf`` and
    pulls exactly the overlapping blocks.  No whole-chromosome download.
  * tree: the ``.nh`` Newick file in ``goldenPath/<db>/<track>/`` whose leaf
    labels match the MAF source names.
  * annotation: ``ncbiRefSeq`` (genePred on UCSC assemblies, bigGenePred on
    GenArk hubs).  RefSeq sequence names (``NC_000011.10``) are translated to
    UCSC names through the ``chromAlias`` table, and the alias is recorded.
  * conservation: a phyloP or phastCons track matched to the alignment
    (``multiz470way`` -> ``phyloP470way``, ``cactus241wayBM`` -> ``phyloP241wayBW``,
    ``multiz124way`` -> ``phyloP124way``), overridable with ``--conservation``.

``--source ensembl`` (default for Ensembl species names such as
``gallus_gallus``)

  * alignment: ``rest.ensembl.org/alignment/region`` for one of the Compara
    multiple-alignment species sets (``mammals``, ``primates``, ``sauropsids``,
    ``fish``, ``amniotes``, ``murinae``) and method (EPO, EPO_EXTENDED, PECAN,
    CACTUS_DB).  Each block carries its own Newick tree with branch lengths,
    and EPO blocks include inferred *ancestral* sequences, which are kept and
    flagged.
  * annotation: ``rest.ensembl.org/overlap/region`` gene/transcript/exon/cds.
  * conservation: GERP constrained elements from the same species set.

Outputs (``<out>/<name>.*``):

  ``.maf``               alignment blocks overlapping the window (whole blocks,
                         not clipped; the reference row is first in each block)
  ``.fa``                the reference sequence of the window, + strand
  ``.nh``                Newick tree(s); one line per block for Ensembl
  ``.annotation.json``   transcripts with exon and CDS intervals, 0-based
                         half-open, on the reference sequence, and what the
                         record declares about CDS completeness
                         (``cds_start_status``, ``cds_end_status``,
                         ``cds_start_frame``; UCSC tracks only)
  ``.conservation.json`` per-base scores (UCSC) or constrained elements (Ensembl)
  ``.manifest.json``     every URL fetched, byte counts, SHA-256 of each output,
                         the sequence-name aliases used, and a per-species
                         alignment coverage table split by annotation class
                         (cds, utr, intron, intergenic)

Examples::

    # human beta-globin, 470-way multiz through the API (bigMaf)
    python3 scripts/data/fetch_window.py --assembly hg38 \
        --locus chr11:5225464-5229395 --flank 500 --track multiz470way \
        --out /tmp/win/HBB

    # fruit fly Adh, 124-way multiz through Range reads on the MAF (wigMaf)
    python3 scripts/data/fetch_window.py --assembly dm6 \
        --locus chr2L:14615552-14618902 --flank 500 --track multiz124way \
        --out /tmp/win/Adh

    # chicken GAPDH, Ensembl Compara sauropsids EPO with ancestral sequences
    python3 scripts/data/fetch_window.py --source ensembl \
        --assembly gallus_gallus --locus 1:76900598-76906236 --flank 500 \
        --species-set sauropsids --out /tmp/win/GAPDH

    # RefSeq sequence names are accepted on the UCSC side
    python3 scripts/data/fetch_window.py --assembly hg38 \
        --locus NC_000011.10:5225464-5229395 --track cactus241wayBM --out /tmp/win/HBB241

Exit status is non-zero if the alignment could not be fetched.  Missing
annotation or conservation is a warning, not an error, because GenArk hubs
have RefSeq but no conservation and some UCSC assemblies have alignments but
no RefSeq.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UCSC_API = "https://api.genome.ucsc.edu"
UCSC_DL = "https://hgdownload.soe.ucsc.edu"
# hgdownload mirrors, tried in order when a host resets or times out; the
# primary reset every connection for several minutes on 2026-09-09 while
# hgdownload2 served the same files.  Override the primary with --download-host.
UCSC_DL_MIRRORS = ["https://hgdownload.soe.ucsc.edu", "https://hgdownload2.soe.ucsc.edu",
                   "https://hgdownload2.gi.ucsc.edu"]
ENSEMBL = "https://rest.ensembl.org"
UA = "relay-annotation-fetch-window/0.1 (https://github.com/nekrut/annotation)"

# alignment track -> conservation track, by convention of the UCSC names
UCSC_CONS_DEFAULT = {
    "multiz470way": "phyloP470wayBW",
    "cactus241wayBM": "phyloP241wayBW",
    "cactus447way": "phyloP447wayBW",
    "multiz100way": "phyloP100way",
    "multiz30way": "phyloP30way",
}


# ---------------------------------------------------------------- utilities
class Fetcher:
    """HTTP GET with a log of everything fetched, polite pacing, retries."""

    def __init__(self, timeout: float, pause: float, log: list) -> None:
        self.timeout = timeout
        self.pause = pause
        self.log = log
        self._last = 0.0

    def get(self, url: str, headers: dict | None = None, what: str = "") -> bytes:
        wait = self.pause - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        last_err: Exception | None = None
        # mirror rotation for hgdownload URLs: attempt k uses mirror k (mod n)
        host = next((m for m in UCSC_DL_MIRRORS if url.startswith(m + "/")), None)
        mirrors = ([host] + [m for m in UCSC_DL_MIRRORS if m != host]) if host else [None]
        for attempt in range(4):
            attempt_url = url
            if host and attempt:
                attempt_url = mirrors[attempt % len(mirrors)] + url[len(host):]
            req = urllib.request.Request(attempt_url, headers={"User-Agent": UA, **(headers or {})})
            try:
                t0 = time.monotonic()
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    data = r.read()
                    status = r.status
                self._last = time.monotonic()
                entry = {"url": attempt_url, "what": what, "status": status, "bytes": len(data),
                         "seconds": round(self._last - t0, 2),
                         "range": (headers or {}).get("Range")}
                if attempt_url != url:
                    entry["mirror_for"] = url
                self.log.append(entry)
                return data
            except urllib.error.HTTPError as e:
                self._last = time.monotonic()
                body = e.read()[:400].decode("utf-8", "replace")
                self.log.append({"url": attempt_url, "what": what, "status": e.code, "bytes": 0, "error": body})
                if e.code in (429, 500, 502, 503, 504) and attempt < 3:
                    time.sleep(2 ** attempt)
                    last_err = e
                    continue
                raise RuntimeError(f"HTTP {e.code} for {url}: {body}") from e
            except (urllib.error.URLError, TimeoutError, OSError) as e:
                last_err = e
                self.log.append({"url": attempt_url, "what": what, "status": 0, "bytes": 0, "error": str(e)[:200]})
                time.sleep(2 ** attempt)
        raise RuntimeError(f"failed after retries: {url}: {last_err}")

    def json(self, url: str, what: str = "") -> dict:
        return json.loads(self.get(url, what=what).decode("utf-8"))



def maf_source(src: str) -> str:
    """Assembly part of a MAF ``s`` row name, spelled as the track's tree
    spells its leaves.  UCSC names rows ``<db>.<chrom>``; GenArk assemblies
    carry a version dot in the db (``GCF_003668045.3.NC_048596.1`` in the
    mm39 35-way) and their tree leaf is ``GCF_003668045v3``."""
    m = re.match(r"^(GC[AF]_\d+)\.(\d+)\.", src)
    if m:
        return f"{m.group(1)}v{m.group(2)}"
    return src.split(".")[0]


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_locus(text: str) -> tuple[str, int, int]:
    """'chr11:5225464-5229395' -> ('chr11', 5225464, 5229395); commas allowed.

    Interpreted as 0-based half-open like the UCSC API and BED.  If you copy
    1-based coordinates from a browser, the window is off by one base at the
    start, which the flank makes irrelevant.
    """
    m = re.fullmatch(r"([^:]+):([\d,]+)-([\d,]+)", text.strip())
    if not m:
        sys.exit(f"cannot parse locus {text!r}; expected chrom:start-end")
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if end <= start:
        sys.exit("locus end must be greater than start")
    return m.group(1), start, end


def newick_leaves(newick: str) -> list[str]:
    return re.findall(r"[(,]\s*([^(),:;\s]+)\s*:", newick)


# ------------------------------------------------------------ UCSC backend
def ucsc_track_url(genome: str, track: str, chrom: str, start: int, end: int, max_items: int) -> str:
    return (f"{UCSC_API}/getData/track?genome={genome};track={track};chrom={chrom};"
            f"start={start};end={end};maxItemsOutput={max_items}")


def ucsc_resolve_chrom(fx: Fetcher, genome: str, chrom: str, manifest: dict) -> str:
    """Accept UCSC, RefSeq, GenBank or Ensembl sequence names; return the
    name the assembly uses, recording the alias table row."""
    chroms = fx.json(f"{UCSC_API}/list/chromosomes?genome={genome}", "chromosome list")
    sizes = chroms.get("chromosomes", {})
    manifest["reference_sequence_count"] = chroms.get("chromCount")
    if chrom in sizes:
        manifest["chrom_size"] = sizes[chrom]
        return chrom
    try:
        alias = fx.json(f"{UCSC_API}/getData/track?genome={genome};track=chromAlias;maxItemsOutput=1000000",
                        "chromAlias table")
    except RuntimeError as e:
        sys.exit(f"{chrom!r} is not a sequence of {genome} and chromAlias is unavailable: {e}")
    rows = alias.get("chromAlias", [])
    if isinstance(rows, dict):  # some assemblies key the table by chromosome
        rows = [r for v in rows.values() for r in v]
    for r in rows:
        if r.get("alias") == chrom:
            manifest["chrom_alias"] = r
            manifest["chrom_size"] = sizes.get(r["chrom"])
            return r["chrom"]
    sys.exit(f"{chrom!r} is neither a sequence name nor an alias in {genome}")


def maf_block_from_api(item: dict) -> str:
    """bigMaf items carry the block as one string with ';' between lines."""
    return item["mafBlock"].replace(";", "\n").rstrip("\n") + "\n"


def maf_blocks_by_range(fx: Fetcher, genome: str, track: str, chrom: str, offsets: list[int],
                        chunk: int = 1 << 20, limit: int = 1 << 28) -> list[str]:
    """Read the MAF blocks that start at the given byte offsets of the
    uncompressed per-chromosome MAF backing a wigMaf track.

    The index the API returns lists every block overlapping the window; in a
    position-sorted MAF those blocks are contiguous in the file, so one
    contiguous byte range from the first offset to the end of the last block
    is fetched in ``chunk``-sized Range requests and then cut at blank lines.
    Blocks that lie between wanted offsets but were not in the index (there
    should be none) are dropped, so the result matches the index exactly."""
    wanted = sorted(set(offsets))
    candidates = [f"{UCSC_DL}/gbdb/{genome}/{track}/maf/{chrom}.maf",
                  f"{UCSC_DL}/gbdb/{genome}/{track}/{chrom}.maf"]
    err = None
    for url in candidates:
        first = wanted[0]
        buf = b""
        pos = first
        done = False
        try:
            while len(buf) < limit and not done:
                data = fx.get(url, headers={"Range": f"bytes={pos}-{pos + chunk - 1}"}, what="maf range")
                if not data:
                    break
                buf += data
                pos += len(data)
                # finished once the block starting at the last wanted offset has ended
                tail = wanted[-1] - first
                if len(buf) > tail and buf.find(b"\n\n", tail) >= 0:
                    done = True
                if len(data) < chunk:
                    break
        except RuntimeError as e:  # 404: try the other layout
            err = e
            continue
        if not buf.startswith(b"a "):
            raise RuntimeError(f"byte {first} of {url} is not the start of a MAF block: {buf[:60]!r}")
        blocks: list[str] = []
        want = set(w - first for w in wanted)
        i = 0
        while i < len(buf):
            j = buf.find(b"\n\n", i)
            if j < 0:
                j = len(buf)
            if i in want:
                blocks.append(buf[i:j + 1].decode("utf-8", "replace").rstrip("\n") + "\n")
            i = j + 2
            # skip any further blank lines
            while i < len(buf) and buf[i:i + 1] == b"\n":
                i += 1
        if len(blocks) != len(wanted):
            raise RuntimeError(f"expected {len(wanted)} blocks from {url}, cut {len(blocks)}; "
                               "the index offsets and the file on hgdownload may be out of step")
        return blocks
    raise RuntimeError(f"no uncompressed MAF for {genome}/{track}/{chrom} on hgdownload: {err}")


def ucsc_alignment(fx: Fetcher, genome: str, track: str, chrom: str, start: int, end: int,
                   max_items: int, manifest: dict) -> list[str]:
    d = fx.json(ucsc_track_url(genome, track, chrom, start, end, max_items), "alignment index/blocks")
    items = d.get(track)
    manifest["alignment_track_type"] = d.get("trackType")
    manifest["alignment_data_time"] = d.get("dataTime")
    if not isinstance(items, list):
        raise RuntimeError(f"track {track} on {genome}: {d.get('error') or d.get('statusMessage') or 'no items'}")
    blocks: list[str] = []
    if items and "mafBlock" in items[0]:
        manifest["alignment_access"] = "api bigMaf"
        blocks = [maf_block_from_api(it) for it in items]
    elif items and "offset" in items[0]:
        manifest["alignment_access"] = "api wigMaf index + HTTP Range on hgdownload /gbdb MAF"
        blocks = maf_blocks_by_range(fx, genome, track, chrom, [int(it["offset"]) for it in items])
    manifest["alignment_blocks"] = len(blocks)
    return blocks


def ucsc_tree(fx: Fetcher, genome: str, track: str, manifest: dict) -> str | None:
    """Pick the Newick file in goldenPath/<db>/<track>/ whose leaf labels are
    the database names used as MAF sources (not common or scientific names)."""
    # the bigMaf track cactus241wayBM downloads from goldenPath/hg38/cactus241way/
    dirs = [track] + ([track[:-2]] if track.endswith(("BM", "BW")) else [])
    names: list[str] = []
    base = ""
    for d in dirs:
        base = f"{UCSC_DL}/goldenPath/{genome}/{d}/"
        try:
            listing = fx.get(base, what="track download directory").decode("utf-8", "replace")
        except RuntimeError:
            continue
        names = sorted(set(re.findall(r'href="([^"/]+\.nh(?:\.txt)?)"', listing)))
        if names:
            break
    if not names:
        return None

    def rank(n: str) -> int:
        n = n.lower()
        if "taxid" in n:
            return 4
        if "common" in n:
            return 3
        if "scientific" in n:
            return 2
        return 0  # e.g. hg38.470way.nh, dm6.124way.sequenceNames.nh

    names.sort(key=rank)
    manifest["tree_candidates"] = names
    tree = fx.get(base + names[0], what="newick tree").decode("utf-8")
    manifest["tree_url"] = base + names[0]
    return tree


def ucsc_annotation(fx: Fetcher, genome: str, chrom: str, start: int, end: int, track: str,
                    manifest: dict) -> list[dict]:
    d = fx.json(ucsc_track_url(genome, track, chrom, start, end, 100000), "annotation")
    items = d.get(track)
    manifest["annotation_track_type"] = d.get("trackType")
    manifest["annotation_data_time"] = d.get("dataTime")
    if not isinstance(items, list):
        raise RuntimeError(f"annotation track {track}: {d.get('error') or d.get('statusMessage') or 'no items'}")
    out = []
    for it in items:
        if "exonStarts" in it:  # genePred
            es = [int(x) for x in it["exonStarts"].rstrip(",").split(",") if x]
            ee = [int(x) for x in it["exonEnds"].rstrip(",").split(",") if x]
            cs, ce = int(it["cdsStart"]), int(it["cdsEnd"])
            tx = {"id": it["name"], "gene": it.get("name2"), "strand": it["strand"],
                  "start": int(it["txStart"]), "end": int(it["txEnd"]),
                  "exons": list(zip(es, ee)), "source": track}
        elif "blockCount" in it:  # bigGenePred (GenArk hubs)
            sizes = [int(x) for x in str(it["blockSizes"]).rstrip(",").split(",") if x]
            starts = [int(x) for x in str(it["chromStarts"]).rstrip(",").split(",") if x]
            es = [int(it["chromStart"]) + s for s in starts]
            ee = [a + b for a, b in zip(es, sizes)]
            cs, ce = int(it["thickStart"]), int(it["thickEnd"])
            tx = {"id": it["name"], "gene": it.get("geneName") or it.get("name2"), "strand": it["strand"],
                  "start": int(it["chromStart"]), "end": int(it["chromEnd"]),
                  "exons": list(zip(es, ee)), "source": track, "type": it.get("geneType")}
        else:
            continue
        tx["cds"] = [(max(a, cs), min(b, ce)) for a, b in tx["exons"] if min(b, ce) > max(a, cs)] if ce > cs else []
        if tx["cds"]:
            tx.update(cds_completeness(it, tx["strand"]))
        out.append(tx)
    return out


def cds_completeness(it: dict, strand: str) -> dict:
    """What a UCSC genePred / bigGenePred record declares about its CDS
    ends.  ``cdsStartStat`` / ``cdsEndStat`` (``cmpl`` / ``incmpl``; ``unk``
    and ``none`` declare nothing, and UCSC's GENCODE tracks carry ``none``
    throughout) come first; GENCODE's ``cds_start_NF`` / ``cds_end_NF``
    tags in bigGenePred ``tag`` next; and the frame of the first coding
    exon (``exonFrames``, transcript orientation) is recorded as
    ``cds_start_frame``, a non-zero value being a third statement of a 5'
    truncation.  The cutter reads all three (cut_windows.py
    ``--partial-ends``).  Keys are absent when nothing is declared."""
    out: dict = {}
    src = {}
    for which, col in (("start", "cdsStartStat"), ("end", "cdsEndStat")):
        v = str(it.get(col) or "").lower()
        if v in ("cmpl", "incmpl"):
            out[f"cds_{which}_status"] = "complete" if v == "cmpl" else "incomplete"
            src[which] = col
    tags = {t.strip() for t in str(it.get("tag") or "").split(",") if t.strip()}
    for which, tag in (("start", "cds_start_NF"), ("end", "cds_end_NF")):
        if tag in tags and f"cds_{which}_status" not in out:
            out[f"cds_{which}_status"] = "incomplete"
            src[which] = "tag"
    frames = [int(x) for x in str(it.get("exonFrames") or "").rstrip(",").split(",") if x.strip().lstrip("-").isdigit()]
    coding = [f for f in frames if f in (0, 1, 2)]
    if coding:
        out["cds_start_frame"] = coding[0] if strand == "+" else coding[-1]
    if src:
        out["cds_status_source"] = src
    return out


def ucsc_conservation(fx: Fetcher, genome: str, track: str, chrom: str, start: int, end: int,
                      manifest: dict) -> dict:
    d = fx.json(ucsc_track_url(genome, track, chrom, start, end, 1000000), "conservation")
    vals = d.get(track)
    if not isinstance(vals, list):
        raise RuntimeError(f"conservation track {track}: {d.get('error') or d.get('statusMessage') or 'no items'}")
    manifest["conservation_track_type"] = d.get("trackType")
    return {"track": track, "kind": "per-base", "values": [[v["start"], v["end"], v["value"]] for v in vals]}


# --------------------------------------------------------- Ensembl backend
def ensembl_alignment(fx: Fetcher, species: str, region: str, species_set: str, method: str,
                      manifest: dict) -> list[dict]:
    url = (f"{ENSEMBL}/alignment/region/{species}/{region}?species_set_group={species_set};"
           f"method={method};content-type=application/json")
    d = fx.json(url, "alignment blocks")
    if isinstance(d, dict):
        raise RuntimeError(f"Ensembl alignment/region: {d.get('error')}")
    manifest["alignment_access"] = f"rest.ensembl.org alignment/region {method} {species_set}"
    manifest["alignment_blocks"] = len(d)
    return d


def ensembl_srclen(fx: Fetcher, cache: dict, species: str, seq_region: str) -> int:
    key = (species, seq_region)
    if key not in cache:
        try:
            info = fx.json(f"{ENSEMBL}/info/assembly/{species}/{seq_region}?content-type=application/json",
                           "sequence length")
            cache[key] = int(info.get("length", 0))
        except RuntimeError:
            cache[key] = 0
    return cache[key]


def ensembl_blocks_to_maf(fx: Fetcher, blocks: list[dict], species: str, manifest: dict) -> tuple[list[str], list[str]]:
    """Ensembl JSON -> MAF text.  Coordinates become 0-based half-open on the
    + strand as MAF requires; sequences on the - strand are reported by
    Ensembl already reverse-complemented, and MAF's start for a '-' row is
    measured from the end of the source sequence."""
    cache: dict = {}
    mafs, trees = [], []
    ancestral = set()
    for b in blocks:
        rows = sorted(b["alignments"], key=lambda a: 0 if a["species"] == species else 1)
        lines = ["a score=0.0"]
        for a in rows:
            sp, reg = a["species"], str(a["seq_region"])
            seq = a["seq"]
            size = len(seq) - seq.count("-")
            is_anc = a.get("description", "") == "" and re.fullmatch(r"[A-Za-z]+(-[A-Za-z]+)*\[\d+\]", sp) is not None
            if is_anc:
                ancestral.add(sp)
                src_len = 0
            else:
                src_len = ensembl_srclen(fx, cache, sp, reg)
            start1, end1, strand = int(a["start"]), int(a["end"]), int(a["strand"])
            if strand == 1 or src_len == 0:
                start0 = start1 - 1
                sym = "+" if strand == 1 else "-"
            else:
                start0 = src_len - end1
                sym = "-"
            lines.append(f"s {sp}.{reg} {start0} {size} {sym} {src_len} {seq}")
        mafs.append("\n".join(lines) + "\n")
        trees.append(b.get("tree", ""))
    manifest["ancestral_sequences"] = sorted(ancestral)
    return mafs, trees


def ensembl_species_set_members(fx: Fetcher, species_set: str, method: str, manifest: dict) -> None:
    """Record the species set's full member list (sorted) in the manifest, so
    ``cut_windows.py`` can give every window of a set the same rows in the
    same order whatever subset happens to align there.  Failure is a warning,
    not an error: the cutter falls back to the observed rows."""
    url = f"{ENSEMBL}/info/compara/species_sets/{method}?content-type=application/json"
    try:
        sets = fx.json(url, "species_sets")
    except (RuntimeError, ValueError) as err:
        manifest["warnings"].append(f"species_sets: {err}")
        return
    for s in sets if isinstance(sets, list) else []:
        if s.get("species_set_group") == species_set and s.get("method", method) == method:
            manifest["species_set_name"] = s.get("name")
            manifest["species_set_members"] = sorted(s.get("species_set", []))
            return
    manifest["warnings"].append(f"species_sets: group {species_set!r} not listed for method {method}")


def ensembl_annotation(fx: Fetcher, species: str, region: str) -> list[dict]:
    url = (f"{ENSEMBL}/overlap/region/{species}/{region}?feature=gene;feature=transcript;"
           f"feature=exon;feature=cds;content-type=application/json")
    feats = fx.json(url, "annotation")
    if isinstance(feats, dict):
        raise RuntimeError(f"Ensembl overlap/region: {feats.get('error')}")
    txs: dict[str, dict] = {}
    genes = {f["id"]: f for f in feats if f["feature_type"] == "gene"}
    for f in feats:
        if f["feature_type"] == "transcript":
            g = genes.get(f.get("Parent"), {})
            txs[f["id"]] = {"id": f["id"], "gene": g.get("external_name") or f.get("Parent"),
                            "strand": "+" if f["strand"] == 1 else "-",
                            "start": f["start"] - 1, "end": f["end"], "exons": [], "cds": [],
                            "source": "ensembl", "type": f.get("biotype")}
    for f in feats:
        p = f.get("Parent")
        if f["feature_type"] == "exon" and p in txs:
            txs[p]["exons"].append((f["start"] - 1, f["end"]))
        elif f["feature_type"] == "cds" and p in txs:
            txs[p]["cds"].append((f["start"] - 1, f["end"]))
    for t in txs.values():
        t["exons"].sort()
        t["cds"].sort()
    return list(txs.values())


def ensembl_constrained(fx: Fetcher, species: str, region: str, species_set: str) -> dict:
    url = (f"{ENSEMBL}/overlap/region/{species}/{region}?feature=constrained;"
           f"species_set_group={species_set};content-type=application/json")
    feats = fx.json(url, "constrained elements")
    if isinstance(feats, dict):
        raise RuntimeError(f"Ensembl constrained elements: {feats.get('error')}")
    return {"track": f"GERP constrained elements ({species_set})", "kind": "elements",
            "values": [[f["start"] - 1, f["end"], f.get("score")] for f in feats]}


def ucsc_sequence(fx: Fetcher, genome: str, chrom: str, start: int, end: int) -> str:
    d = fx.json(f"{UCSC_API}/getData/sequence?genome={genome};chrom={chrom};start={start};end={end}", "reference sequence")
    return d.get("dna", "")


def ensembl_sequence(fx: Fetcher, species: str, region: str) -> str:
    d = fx.json(f"{ENSEMBL}/sequence/region/{species}/{region}?content-type=application/json", "reference sequence")
    return d.get("seq", "")


# --------------------------------------------------- coverage by annotation
def classify_window(start: int, end: int, transcripts: list[dict]) -> list[str]:
    """Per reference base: cds > utr > intron > intergenic, union over transcripts."""
    n = end - start
    cls = ["intergenic"] * n
    rank = {"intergenic": 0, "intron": 1, "utr": 2, "cds": 3}

    def paint(a: int, b: int, label: str) -> None:
        for i in range(max(a, start), min(b, end)):
            if rank[label] > rank[cls[i - start]]:
                cls[i - start] = label

    for t in transcripts:
        paint(t["start"], t["end"], "intron")
        for a, b in t["exons"]:
            paint(a, b, "utr")
        for a, b in t["cds"]:
            paint(a, b, "cds")
    return cls


def coverage_by_species(maf_blocks: list[str], ref_prefix: str, start: int, end: int,
                        classes: list[str]) -> dict:
    """Fraction of window reference bases, per annotation class, at which each
    other source has an aligned (non-gap) base."""
    n = end - start
    covered: dict[str, bytearray] = {}
    for block in maf_blocks:
        rows = [ln.split() for ln in block.splitlines() if ln.startswith("s ")]
        if not rows:
            continue
        ref = next((r for r in rows if maf_source(r[1]) == ref_prefix), None)
        if ref is None:
            continue
        ref_seq = ref[6]
        ref_pos = int(ref[2])
        if ref[4] == "-":  # reference rows in these tracks are always +, keep it simple
            continue
        # map alignment columns to reference coordinates
        col_to_ref = []
        p = ref_pos
        for ch in ref_seq:
            col_to_ref.append(p if ch != "-" else -1)
            if ch != "-":
                p += 1
        for r in rows:
            if r is ref:
                continue
            sp = maf_source(r[1])
            arr = covered.setdefault(sp, bytearray(n))
            seq = r[6]
            for col, refp in enumerate(col_to_ref):
                if refp < start or refp >= end or refp < 0:
                    continue
                if seq[col] not in "-":
                    arr[refp - start] = 1
    totals = {c: classes.count(c) for c in ("cds", "utr", "intron", "intergenic")}
    table = {}
    for sp, arr in covered.items():
        counts = {c: 0 for c in totals}
        for i, hit in enumerate(arr):
            if hit:
                counts[classes[i]] += 1
        table[sp] = {c: (round(counts[c] / totals[c], 3) if totals[c] else None) for c in totals}
        table[sp]["all"] = round(sum(counts.values()) / n, 3)
    return {"window_bases_by_class": totals, "species": table}


# ------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", required=True,
                    help="UCSC db (hg38, dm6), GenArk accession (GCF_000002765.6) or Ensembl species (gallus_gallus)")
    ap.add_argument("--locus", required=True, help="chrom:start-end; RefSeq/GenBank names accepted on the UCSC side")
    ap.add_argument("--flank", type=int, default=0, help="bases added on each side")
    ap.add_argument("--source", choices=["auto", "ucsc", "ensembl"], default="auto")
    ap.add_argument("--track", default=None,
                    help="UCSC alignment track (default: multiz470way on hg38, else the first *way track); "
                         "'none' skips the alignment, for GenArk assemblies that have RefSeq but no alignment")
    ap.add_argument("--annotation", default="ncbiRefSeq", help="UCSC annotation track")
    ap.add_argument("--conservation", default=None, help="UCSC conservation track; 'none' to skip")
    ap.add_argument("--species-set", default="mammals", help="Ensembl Compara species_set_group")
    ap.add_argument("--method", default="EPO", help="Ensembl Compara method: EPO, EPO_EXTENDED, PECAN, CACTUS_DB")
    ap.add_argument("--out", required=True, help="output directory (outside the repository)")
    ap.add_argument("--name", default=None, help="output file stem (default: assembly_chrom_start_end)")
    ap.add_argument("--max-items", type=int, default=100000, help="UCSC maxItemsOutput for alignment blocks (API maximum 1000000)")
    ap.add_argument("--timeout", type=float, default=120.0)
    ap.add_argument("--pause", type=float, default=0.34, help="minimum seconds between requests (NCBI/UCSC politeness)")
    ap.add_argument("--download-host", default=None,
                    help="primary hgdownload host, e.g. https://hgdownload2.soe.ucsc.edu (mirrors are tried on connection errors)")
    ap.add_argument("--dry-run", action="store_true", help="print what would be fetched and exit")
    args = ap.parse_args()
    if args.download_host:
        global UCSC_DL
        UCSC_DL = args.download_host.rstrip("/")
        if UCSC_DL in UCSC_DL_MIRRORS:
            UCSC_DL_MIRRORS.remove(UCSC_DL)
        UCSC_DL_MIRRORS.insert(0, UCSC_DL)

    chrom, s, e = parse_locus(args.locus)
    start, end = max(0, s - args.flank), e + args.flank
    source = args.source
    if source == "auto":
        source = "ensembl" if re.fullmatch(r"[a-z]+_[a-z_]+", args.assembly) else "ucsc"
    name = args.name or f"{args.assembly}_{chrom}_{start}_{end}".replace("/", "_")
    manifest: dict = {"tool": "scripts/data/fetch_window.py", "version": "0.1", "source": source,
                      "assembly": args.assembly, "locus_requested": args.locus, "flank": args.flank,
                      "window": {"chrom": chrom, "start": start, "end": end, "coordinates": "0-based half-open"},
                      "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "requests": [],
                      "warnings": []}
    fx = Fetcher(args.timeout, args.pause, manifest["requests"])
    if args.dry_run:
        print(json.dumps({k: v for k, v in manifest.items() if k != "requests"}, indent=2))
        return 0
    os.makedirs(args.out, exist_ok=True)
    stem = os.path.join(args.out, name)

    maf_blocks: list[str] = []
    trees: list[str] = []
    sequence = ""
    transcripts: list[dict] = []
    cons: dict | None = None

    if source == "ucsc":
        genome = args.assembly
        chrom = ucsc_resolve_chrom(fx, genome, chrom, manifest)
        manifest["window"]["chrom"] = chrom
        track = args.track
        if track is None:
            tracks = fx.json(f"{UCSC_API}/list/tracks?genome={genome};trackLeavesOnly=1", "track list")
            leaves = next((v for k, v in tracks.items() if k == genome or k.startswith("GC")), {})
            way = [k for k, v in leaves.items() if isinstance(v, dict)
                   and str(v.get("type", "")).split()[0] in ("bigMaf", "wigMaf")]
            if not way:
                print(f"error: {genome} has no multiple-alignment track", file=sys.stderr)
                return 2
            track = "multiz470way" if "multiz470way" in way else sorted(way)[0]
            manifest["alignment_tracks_available"] = sorted(way)
        manifest["alignment_track"] = track
        if track == "none":
            manifest["alignment_access"] = "skipped (--track none)"
        else:
            try:
                maf_blocks = ucsc_alignment(fx, genome, track, chrom, start, end, args.max_items, manifest)
            except RuntimeError as err:
                print(f"error: {err}", file=sys.stderr)
                return 2
            tree = ucsc_tree(fx, genome, track, manifest)
            if tree:
                trees = [tree.strip()]
            else:
                manifest["warnings"].append("no Newick tree found in the track download directory")
        try:
            sequence = ucsc_sequence(fx, genome, chrom, start, end)
        except RuntimeError as err:
            manifest["warnings"].append(f"sequence: {err}")
        try:
            transcripts = ucsc_annotation(fx, genome, chrom, start, end, args.annotation, manifest)
            manifest["annotation_track"] = args.annotation
        except RuntimeError as err:
            manifest["warnings"].append(f"annotation: {err}")
        cons_track = args.conservation
        if cons_track is None:
            m = re.search(r"(\d+)way", track)
            cons_track = UCSC_CONS_DEFAULT.get(track, f"phyloP{m.group(1)}way" if m else "none")
        if cons_track != "none" and track != "none":
            try:
                cons = ucsc_conservation(fx, genome, cons_track, chrom, start, end, manifest)
                manifest["conservation_track"] = cons_track
            except RuntimeError as err:
                manifest["warnings"].append(f"conservation: {err}")
        ref_prefix = genome
    else:
        species = args.assembly
        region = f"{chrom}:{start + 1}-{end}"  # Ensembl is 1-based closed
        manifest["ensembl_region"] = region
        manifest["species_set"] = args.species_set
        manifest["method"] = args.method
        ensembl_species_set_members(fx, args.species_set, args.method, manifest)
        try:
            blocks = ensembl_alignment(fx, species, region, args.species_set, args.method, manifest)
        except RuntimeError as err:
            print(f"error: {err}", file=sys.stderr)
            return 2
        maf_blocks, trees = ensembl_blocks_to_maf(fx, blocks, species, manifest)
        try:
            sequence = ensembl_sequence(fx, species, region)
        except RuntimeError as err:
            manifest["warnings"].append(f"sequence: {err}")
        try:
            transcripts = ensembl_annotation(fx, species, region)
            manifest["annotation_track"] = "ensembl overlap/region"
        except RuntimeError as err:
            manifest["warnings"].append(f"annotation: {err}")
        if args.conservation != "none":
            try:
                cons = ensembl_constrained(fx, species, region, args.species_set)
            except RuntimeError as err:
                manifest["warnings"].append(f"constrained elements: {err}")
        ref_prefix = species

    # ---- write outputs
    with open(stem + ".maf", "w") as f:
        f.write(f"##maf version=1 scoring=none\n# {source} {args.assembly} {chrom}:{start}-{end} "
                f"track={manifest.get('alignment_track', manifest.get('alignment_access'))}\n\n")
        for b in maf_blocks:
            f.write(b.rstrip("\n") + "\n\n")
    with open(stem + ".fa", "w") as f:
        f.write(f">{args.assembly}.{chrom}:{start}-{end} 0-based half-open, + strand\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i + 80] + "\n")
    with open(stem + ".nh", "w") as f:
        for t in trees:
            f.write(t.rstrip("\n") + "\n")
    with open(stem + ".annotation.json", "w") as f:
        json.dump({"chrom": chrom, "window": [start, end], "coordinates": "0-based half-open",
                   "transcripts": transcripts}, f, indent=1)
    with open(stem + ".conservation.json", "w") as f:
        json.dump(cons or {"track": None, "values": []}, f)

    # ---- summary statistics
    sources = sorted({maf_source(ln.split()[1]) for b in maf_blocks for ln in b.splitlines() if ln.startswith("s ")})
    manifest["alignment_sources"] = len(sources)
    if trees and source == "ucsc":
        leaves = set(newick_leaves(trees[0]))
        manifest["tree_leaves"] = len(leaves)
        manifest["sources_not_in_tree"] = sorted(set(sources) - leaves)[:20]
    classes = classify_window(start, end, transcripts)
    manifest["coverage"] = coverage_by_species(maf_blocks, ref_prefix, start, end, classes)
    manifest["transcripts"] = len(transcripts)
    manifest["outputs"] = {}
    for suffix in (".maf", ".fa", ".nh", ".annotation.json", ".conservation.json"):
        p = stem + suffix
        manifest["outputs"][os.path.basename(p)] = {"bytes": os.path.getsize(p), "sha256": sha256_file(p)}
    manifest["bytes_downloaded"] = sum(r["bytes"] for r in manifest["requests"])
    with open(stem + ".manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)

    cov = manifest["coverage"]
    print(f"{name}: {len(maf_blocks)} blocks, {len(sources)} sources, {len(transcripts)} transcripts, "
          f"{manifest['bytes_downloaded'] / 1e6:.1f} MB downloaded in {len(manifest['requests'])} requests")
    print("window bases by class:", cov["window_bases_by_class"])
    for w in manifest["warnings"]:
        print("warning:", w, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
