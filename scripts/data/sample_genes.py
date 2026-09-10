#!/usr/bin/env python3
"""Draw a fixed-seed, non-overlapping sample of protein-coding genes from a
UCSC RefSeq track on one chromosome, for the coverage runs of
``docs/data-sources.md`` section 6.1.

Standard library only.  Reads ``ncbiRefSeqCurated`` (or another genePred
track) through the UCSC JSON API, keeps NM_ transcripts whose CDS is
complete at both ends, with at least ``--min-exons`` exons and a transcript
length inside ``--length``, takes the longest such transcript per gene,
drops genes whose transcripts overlap any other kept gene, shuffles with
``--seed`` and prints the first ``--n`` as one line per gene::

    gene<TAB>transcript<TAB>chrom:start-end (0-based half-open)<TAB>exons

The same command with the same seed reproduces the same list as long as
the track on the API host has not changed (the manifest of every fetched
window records the track's ``dataTime``).

Example (the mouse and worm samples of section 6.1)::

    python3 scripts/data/sample_genes.py --genome mm39 --chrom chr19 --length 3000-12000 --n 10
    python3 scripts/data/sample_genes.py --genome ce11 --chrom chrIII --length 2000-8000 --n 10
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import urllib.request

API = "https://api.genome.ucsc.edu"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genome", required=True, help="UCSC database, e.g. mm39")
    ap.add_argument("--chrom", required=True, help="one chromosome, e.g. chr19")
    ap.add_argument("--track", default="ncbiRefSeqCurated", help="genePred track (default ncbiRefSeqCurated)")
    ap.add_argument("--length", default="3000-12000", help="transcript length range, inclusive, e.g. 2000-8000")
    ap.add_argument("--min-exons", type=int, default=3)
    ap.add_argument("--n", type=int, default=10, help="genes to print")
    ap.add_argument("--seed", type=int, default=20260909)
    ap.add_argument("--max-items", type=int, default=200000)
    a = ap.parse_args()
    lo, hi = (int(x) for x in a.length.split("-"))

    url = f"{API}/getData/track?genome={a.genome};track={a.track};chrom={a.chrom};maxItemsOutput={a.max_items}"
    with urllib.request.urlopen(url, timeout=180) as r:
        d = json.load(r)
    items = d.get(a.track)
    if not isinstance(items, list):
        print(f"no items for {a.track} on {a.genome} {a.chrom}: {d.get('error') or d.get('statusMessage')}", file=sys.stderr)
        return 1
    if d.get("itemsReturned") == a.max_items:
        print(f"warning: itemsReturned == max-items ({a.max_items}); raise --max-items", file=sys.stderr)
    print(f"# {a.genome} {a.chrom} {a.track} dataTime={d.get('dataTime')} items={len(items)} "
          f"length={lo}-{hi} min_exons={a.min_exons} seed={a.seed}")

    # one transcript per gene: longest NM_ with complete CDS and enough exons
    best: dict[str, dict] = {}
    for it in items:
        name = it.get("name", "")
        if not name.startswith("NM_"):
            continue
        if it.get("cdsStartStat") != "cmpl" or it.get("cdsEndStat") != "cmpl":
            continue
        if int(it.get("exonCount", 0)) < a.min_exons:
            continue
        length = int(it["txEnd"]) - int(it["txStart"])
        if not lo <= length <= hi:
            continue
        gene = it.get("name2") or name
        if gene not in best or length > int(best[gene]["txEnd"]) - int(best[gene]["txStart"]):
            best[gene] = it

    # drop genes overlapping any other transcript of the track on either strand
    spans = sorted((int(it["txStart"]), int(it["txEnd"]), it.get("name2") or it["name"]) for it in items)
    kept = []
    for gene, it in best.items():
        s, e = int(it["txStart"]), int(it["txEnd"])
        clash = any(os < e and oe > s and og != gene for os, oe, og in spans)
        if not clash:
            kept.append((gene, it))
    kept.sort(key=lambda x: (int(x[1]["txStart"]), x[0]))
    print(f"# candidates after filters: {len(best)}; non-overlapping: {len(kept)}")
    random.Random(a.seed).shuffle(kept)
    for gene, it in kept[: a.n]:
        print(f"{gene}\t{it['name']}\t{a.chrom}:{it['txStart']}-{it['txEnd']}\t{it['exonCount']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
