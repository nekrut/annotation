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

``--short-gaps`` also lists every exon gap shorter than the cutter's
``MIN_INTRON`` (20, the benchmark scorer's floor) with two flanking exon
bases fetched from the API's sequence endpoint (one request per gap), the
transcripts stating it, whether both flanks are CDS ends, and engels'
overlapping motif-window test (``cut_windows.motif_windows``): whether a
motif-masked decoder could read GT/GC..AG across the gap by borrowing an
exon base on each side, with the donor and acceptor windows reported
apart and a class naming which side failed.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cut_windows import LENGTH_FLOORS, MIN_INTRON, feature_lengths, short_gaps  # noqa: E402
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


def fetch_sequence(genome: str, chrom: str, start: int, end: int, timeout: float) -> str:
    """Reference bases of [start, end) (0-based half-open), upper case."""
    url = f"{API}/getData/sequence?genome={genome};chrom={chrom};start={start};end={end}"
    with urllib.request.urlopen(url, timeout=timeout) as r:
        d = json.loads(r.read())
    return str(d["dna"]).upper()


def short_gap_table(genome: str, chrom: str, txs: list[dict], flank: int, pause: float, timeout: float) -> dict:
    """``cut_windows.short_gaps`` over a chromosome, with the flanks of every
    gap fetched separately so no chromosome sequence is downloaded."""
    rec = short_gaps(txs, "", 0, 0, MIN_INTRON, flank)  # no sequence yet: every gap outside_window
    from cut_windows import COMP, motif_fields, motif_summary
    for g in rec["gaps"]:
        time.sleep(pause)
        seq = fetch_sequence(genome, chrom, g["start"] - flank, g["end"] + flank, timeout)
        left, gap, right = seq[:flank], seq[flank:flank + g["length"]], seq[flank + g["length"]:]
        if g["strand"] == "-":
            left, gap, right = right.translate(COMP)[::-1], gap.translate(COMP)[::-1], left.translate(COMP)[::-1]
        g.update({"left": left, "gap": gap, "right": right, **motif_fields(left, gap, right)})
    rec.update(motif_summary(rec["gaps"]))
    rec["outside_window"] = 0
    rec["sequence_requests"] = len(rec["gaps"])
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("targets", nargs="+", help="genome:chrom[:track]")
    ap.add_argument("--all", action="store_true", help="keep transcripts without a CDS too")
    ap.add_argument("--max-items", type=int, default=200000)
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between requests")
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--json", default=None, help="also write every record here")
    ap.add_argument("--short-gaps", action="store_true",
                    help=f"list every exon gap shorter than {MIN_INTRON} with its flanks (one sequence request each)")
    ap.add_argument("--flank", type=int, default=2, help="exon bases on each side of a short gap to fetch")
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
        if a.short_gaps:
            rec["short_gaps"] = sg = short_gap_table(genome, chrom, txs, a.flank, a.pause, a.timeout)
            print(f"# {genome} {chrom}: {sg['n']} exon gaps under {MIN_INTRON} bases, {sg['in_cds']} between two CDS "
                  f"ends, lengths {sg['by_length']}, motif_window {sg['motif_window']} "
                  f"(borrowing an exon base {sg['motif_borrows_exon_base']}), motif_exact {sg['motif_exact']}, "
                  f"by class {sg['motif_by_class']}, unresolved windows {sg['motif_unresolved_windows']}, "
                  f"gaps with a base outside ACGT {sg['gaps_with_unresolved_bases']}",
                  file=sys.stderr)
            if a.markdown:
                print("| Assembly | Chrom | Strand | Gap (0-based, half-open) | Length | In CDS | Transcripts | "
                      "Flank / gap / flank | Donor window | Acceptor window | Class |")
                print("|---|---|---|---|---|---|---|---|---|---|---|")
                for g in sg["gaps"]:
                    print(f"| {genome} | {chrom} | {g['strand']} | {g['start']}-{g['end']} | {g['length']} | "
                          f"{'yes' if g['in_cds'] else 'no'} | {', '.join(str(t) for t in g['transcripts'])} | "
                          f"`{g['left']} {g['gap']} {g['right']}` | "
                          f"{ {True: 'yes', False: 'no', None: 'unresolved'}[g['motif_donor']] } | "
                          f"{ {True: 'yes', False: 'no', None: 'unresolved'}[g['motif_acceptor']] } | "
                          f"{g['motif_class'].replace('_', ' ')} |")
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
