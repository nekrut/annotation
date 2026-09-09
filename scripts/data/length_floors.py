#!/usr/bin/env python3
"""Count how many of a chromosome's annotated introns, CDSs and transcript
spans lie under the length floors a decoder cannot represent
(``LENGTH_FLOORS`` in cut_windows.py: introns shorter than 30 and 50
bases, CDSs shorter than 60, spans of at most 80), for
``docs/data-sources.md`` section 6.3.  Standard library only.

Reads one UCSC genePred or bigGenePred track for one chromosome through
the JSON API, converts every record with ``fetch_window.ucsc_item_to_transcript``
(so the CDS bounds are the ones the record states, not the exon bounds),
keeps the transcripts with a CDS (``--all`` keeps every transcript) and
tabulates ``cut_windows.feature_lengths`` over the whole chromosome, where
nothing is clipped.  Introns are distinct (strand, start, end) intervals,
so an intron shared by isoforms counts once; CDS and span are per
transcript.  ``--markdown`` prints one table row per chromosome::

    python3 scripts/data/length_floors.py --markdown \\
        hg38:chr21:ncbiRefSeqCurated dm6:chr2L:ncbiRefSeqCurated \\
        ce11:chrIII:ncbiRefSeqCurated sacCer3:chrIV:ncbiRefSeq \\
        GCF_000002765.6:NC_037283.1:ncbiRefSeq

Each argument is ``genome:chrom[:track]`` (track default ``ncbiRefSeqCurated``).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cut_windows import LENGTH_FLOORS, feature_lengths  # noqa: E402
from fetch_window import ucsc_item_to_transcript  # noqa: E402

API = "https://api.genome.ucsc.edu"


def fetch_track(genome: str, chrom: str, track: str, max_items: int, timeout: float) -> tuple[list[dict], dict]:
    url = f"{API}/getData/track?genome={genome};track={track};chrom={chrom};maxItemsOutput={max_items}"
    t0 = time.time()
    with urllib.request.urlopen(url, timeout=timeout) as r:
        raw = r.read()
    d = json.loads(raw)
    items = d.get(track)
    if not isinstance(items, list):
        raise RuntimeError(f"{genome} {chrom} {track}: {d.get('error') or d.get('statusMessage') or 'no items'}")
    if d.get("itemsReturned") == max_items:
        raise RuntimeError(f"{genome} {chrom} {track}: itemsReturned == max-items; raise --max-items")
    meta = {"url": url, "bytes": len(raw), "seconds": round(time.time() - t0, 1),
            "dataTime": d.get("dataTime"), "trackType": d.get("trackType"), "items": len(items)}
    return items, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("targets", nargs="+", help="genome:chrom[:track]")
    ap.add_argument("--all", action="store_true", help="keep transcripts without a CDS too")
    ap.add_argument("--max-items", type=int, default=200000)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between requests")
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--json", default=None, help="also write every record here")
    a = ap.parse_args()
    fi, fc, fs = LENGTH_FLOORS["intron"], LENGTH_FLOORS["cds"], LENGTH_FLOORS["span"]
    if a.markdown:
        head = (["Assembly", "Chromosome", "Track", "Coding transcripts", "Single-exon", "Distinct introns",
                 "Shortest intron", "Median intron"] + [f"Introns < {f}" for f in fi]
                + [f"CDS < {f}" for f in fc] + [f"Span < {f}" for f in fs])
        print("| " + " | ".join(head) + " |")
        print("|" + "---|" * len(head))
    records = []
    for k, target in enumerate(a.targets):
        parts = target.split(":")
        genome, chrom = parts[0], parts[1]
        track = parts[2] if len(parts) > 2 else "ncbiRefSeqCurated"
        if k:
            time.sleep(a.pause)
        items, meta = fetch_track(genome, chrom, track, a.max_items, a.timeout)
        txs = [t for t in (ucsc_item_to_transcript(it, track) for it in items) if t is not None]
        skipped = len(items) - len(txs)
        if not a.all:
            txs = [t for t in txs if t["cds"]]
        fl = feature_lengths(txs, 0, 1 << 40)
        single = sum(1 for t in txs if len(t["exons"]) == 1)
        rec = {"genome": genome, "chrom": chrom, "track": track, "fetch": meta, "items_skipped": skipped,
               "transcripts": len(txs), "single_exon": single, "floors": {k: list(v) for k, v in LENGTH_FLOORS.items()},
               "lengths": fl}
        records.append(rec)
        if a.markdown:
            i, c, s = fl["introns"], fl["cds"], fl["span"]
            pct = lambda n, d: f"{n:,} ({100 * n / d:.2f}%)" if d else "0"
            row = [genome, chrom, track, f"{len(txs):,}", pct(single, len(txs)), f"{i['n']:,}",
                   str(i["min"]), str(i["median"])] + [pct(i[f"below_{f}"], i["n"]) for f in fi] \
                + [pct(c[f"below_{f}"], c["n"]) for f in fc] + [pct(s[f"below_{f}"], s["n"]) for f in fs]
            print("| " + " | ".join(row) + " |")
        else:
            print(json.dumps(rec))
        print(f"# {genome} {chrom} {track}: {meta['items']} items, {meta['bytes']:,} bytes, {meta['seconds']} s, "
              f"dataTime {meta['dataTime']}, {skipped} items of another shape skipped", file=sys.stderr)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(records, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
