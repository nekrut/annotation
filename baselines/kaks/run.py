#!/usr/bin/env python3
"""The single command of T-human-010: run the KA/KS sliding-window
classifier over a seeded gene sample for every species pair in a table,
sweep window size and informant distance, and write the results.

Standard library only (Python 3.11).  From a fresh clone::

    python3 baselines/kaks/run.py --cache /tmp/kaks-cache

does, for every row of ``baselines/kaks/pairs.tsv``:

1. draw the gene sample with ``scripts/data/sample_genes.py`` (fixed seed,
   one chromosome, non-overlapping RefSeq genes with complete CDS);
2. fetch each gene's window with ``scripts/data/fetch_window.py`` into
   ``<cache>/windows/<pair>/<gene>.*`` (skipped when the manifest exists,
   so a rerun is offline);
3. run ``windows.run_one`` for every informant of the row and every window
   size in ``--windows``, and pool the confusion counts over the genes.

Outputs, under ``--out`` (default ``baselines/kaks/results``):

``summary.tsv``
    one row per pair x informant x window size, pooled over genes:
    strand-aware nucleotide sensitivity, precision, F1 and MCC
    (``docs/benchmark.md`` 4.1), the same restricted to aligned bases, the
    aligned fraction, window status counts, and sensitivity by CDS segment
    length.
``per_gene.tsv``
    the same per gene, so the spread behind each pooled number is visible.
``summary.md``
    the sweep as Markdown tables, one per pair.
``manifest.json``
    command line, sample headers (track ``dataTime``), SHA-256 of every
    fetched input, informants that were absent from a window, and the
    tool versions.

``pairs.tsv`` columns: ``pair``, ``reference`` (UCSC db), ``track``,
``chrom``, ``informants`` (comma-separated MAF source names; a missing
one is recorded, not fatal), ``length`` (transcript length range for the
sample), ``n`` (genes), ``flank`` (bases on each side), ``table`` (genetic
code), ``distance_note`` (free text).

``--limit-genes`` and ``--pairs-filter`` bound a run; ``--offline``
refuses to fetch and uses the cache only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts", "data"))
import windows as wm  # noqa: E402
import cut_windows as cw  # noqa: E402
import kaks  # noqa: E402

TOOL_VERSION = "0.1"
SAMPLE = os.path.join(ROOT, "scripts", "data", "sample_genes.py")
FETCH = os.path.join(ROOT, "scripts", "data", "fetch_window.py")
METRIC_KEYS = ("tp", "fp", "fn", "tn")


def read_pairs(path: str) -> list[dict]:
    rows = []
    with open(path) as fh:
        header = None
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def sha256(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def sample_genes(pair: dict, cache: str, offline: bool, seed: int) -> tuple[list[dict], str]:
    """The gene sample of one pair, cached as ``<cache>/samples/<pair>.tsv``."""
    d = os.path.join(cache, "samples")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, f"{pair['pair']}.tsv")
    if not os.path.exists(path):
        if offline:
            raise SystemExit(f"--offline and no sample at {path}")
        cmd = [sys.executable, SAMPLE, "--genome", pair["reference"], "--chrom", pair["chrom"],
               "--length", pair["length"], "--n", pair["n"], "--seed", str(seed)]
        for attempt in range(3):  # the API resets connections now and then
            out = subprocess.run(cmd, capture_output=True, text=True)
            if out.returncode == 0:
                break
            time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError(f"sample_genes failed for {pair['pair']}: {out.stderr[-400:]}")
        with open(path, "w") as fh:
            fh.write(out.stdout)
    genes, header = [], []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                header.append(ln.strip())
                continue
            g, tx, locus, exons = ln.rstrip("\n").split("\t")
            genes.append({"gene": g, "transcript": tx, "locus": locus, "exons": int(exons)})
    return genes, " | ".join(header)


def fetch_window(pair: dict, gene: dict, cache: str, offline: bool, pause: float) -> str | None:
    """Fetch one gene's window if it is not cached; return the stem or None."""
    d = os.path.join(cache, "windows", pair["pair"])
    os.makedirs(d, exist_ok=True)
    stem = os.path.join(d, gene["gene"])
    if os.path.exists(stem + ".manifest.json") and os.path.exists(stem + ".maf"):
        return stem
    if offline:
        return None
    cmd = [sys.executable, FETCH, "--assembly", pair["reference"], "--locus", gene["locus"],
           "--flank", pair.get("flank", "500"), "--track", pair["track"], "--conservation", "none",
           "--out", d, "--name", gene["gene"], "--pause", str(pause)]
    for attempt in range(2):
        out = subprocess.run(cmd, capture_output=True, text=True)
        if out.returncode == 0 and os.path.exists(stem + ".maf"):
            break
        time.sleep(5 * (attempt + 1))
    if out.returncode != 0 or not os.path.exists(stem + ".maf"):
        sys.stderr.write(f"fetch failed for {pair['pair']} {gene['gene']} {gene['locus']}: {out.stderr[-500:]}\n")
        return None
    return stem


def zero_metrics() -> dict:
    return {k: 0 for k in METRIC_KEYS}


def pooled_strata(acc: dict, strata: dict) -> None:
    for k, v in strata.items():
        a = acc.setdefault(k, {"tp": 0, "fn": 0})
        a["tp"] += v["tp"]
        a["fn"] += v["fn"]


def finish_strata(acc: dict) -> dict:
    return {k: {**v, "sensitivity": (v["tp"] / (v["tp"] + v["fn"])) if v["tp"] + v["fn"] else None}
            for k, v in acc.items()}


def fmt(v, nd=3) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def run_pair(pair: dict, genes: list[dict], stems: dict, sizes: list[int], a, manifest: dict) -> tuple[list[dict], list[dict]]:
    """Per-gene rows and pooled rows of one pair."""
    informants = [s.strip() for s in pair["informants"].split(",") if s.strip()]
    table = int(pair.get("table", "1") or 1)
    per_gene, pooled = [], {}
    absent = manifest.setdefault("informants_absent", {})
    for gene in genes:
        stem = stems.get(gene["gene"])
        if not stem:
            continue
        mf, ref, chrom, start, end, seq, blocks, rows, dist, transcripts = wm.load(stem, None, set(), a.transcript_types)
        present = [x for x in informants if x in rows]
        for x in informants:
            if x not in rows:
                absent.setdefault(f"{pair['pair']}:{x}", []).append(gene["gene"])
        for ext in (".maf", ".fa", ".annotation.json", ".manifest.json"):
            manifest["inputs"][os.path.relpath(stem + ext, a.cache)] = sha256(stem + ext)
        if not present:
            continue
        inf_all, ins_all, _ = cw.paint_informants(blocks, ref, start, end, present)
        truth = wm.truth_arrays(transcripts, start, end)
        cds_bases = sum(truth[0]) + sum(truth[1])
        for k, informant in enumerate(present):
            inf, ins = wm.informant_string(inf_all[k], ins_all[k])
            for W in sizes:
                S = a.step or max(3, W // 3)
                res, _ = wm.run_one(seq, inf, ins, start, W, S, a.alpha, a.min_codons, a.max_ref_stops, truth,
                                    table, not a.no_bonferroni)
                m, ma = res["nucleotide"], res["nucleotide_aligned_only"]
                row = {"pair": pair["pair"], "reference": ref, "informant": informant, "distance": dist.get(informant),
                       "gene": gene["gene"], "chrom": chrom, "start": start, "end": end, "window_bp": end - start,
                       "cds_bases": cds_bases, "W": W, "step": S, "aligned_fraction": res["aligned_fraction"],
                       **{f"win_{s}": res["window_status"][s] for s in wm.STATUSES},
                       **{k2: m[k2] for k2 in METRIC_KEYS},
                       "sensitivity": m["sensitivity"], "precision": m["precision"], "f1": m["f1"], "mcc": m["mcc"],
                       **{f"aligned_{k2}": ma[k2] for k2 in METRIC_KEYS},
                       "aligned_sensitivity": ma["sensitivity"], "aligned_precision": ma["precision"],
                       "aligned_f1": ma["f1"], "aligned_mcc": ma["mcc"],
                       "strata": res["sensitivity_by_cds_segment_length"]}
                per_gene.append(row)
                key = (informant, W)
                p = pooled.setdefault(key, {"pair": pair["pair"], "reference": ref, "informant": informant,
                                            "distance": dist.get(informant), "W": W, "step": S, "genes": 0,
                                            "window_bp": 0, "cds_bases": 0, "aligned_bases": 0,
                                            **{f"win_{s}": 0 for s in wm.STATUSES},
                                            "all": zero_metrics(), "aligned": zero_metrics(), "strata": {}})
                p["genes"] += 1
                p["window_bp"] += end - start
                p["cds_bases"] += cds_bases
                p["aligned_bases"] += round(res["aligned_fraction"] * (end - start))
                for s in wm.STATUSES:
                    p[f"win_{s}"] += res["window_status"][s]
                p["all"] = wm.add(p["all"], m)
                p["aligned"] = wm.add(p["aligned"], ma)
                pooled_strata(p["strata"], res["sensitivity_by_cds_segment_length"])
    rows = []
    for key in sorted(pooled, key=lambda k: ((pooled[k]["distance"] if pooled[k]["distance"] is not None else 9e9), k[1])):
        p = pooled[key]
        m, ma = wm.metrics(p["all"]), wm.metrics(p["aligned"])
        rows.append({k2: v for k2, v in p.items() if k2 not in ("all", "aligned", "strata", "aligned_bases")} |
                    {"aligned_fraction": p["aligned_bases"] / p["window_bp"] if p["window_bp"] else None}
                    | {k2: m[k2] for k2 in METRIC_KEYS}
                    | {"sensitivity": m["sensitivity"], "precision": m["precision"], "f1": m["f1"], "mcc": m["mcc"]}
                    | {f"aligned_{k2}": ma[k2] for k2 in METRIC_KEYS}
                    | {"aligned_sensitivity": ma["sensitivity"], "aligned_precision": ma["precision"],
                       "aligned_f1": ma["f1"], "aligned_mcc": ma["mcc"], "strata": finish_strata(p["strata"])})
    return per_gene, rows


STRATA_COLS = ("0-60", "60-100", "100-200", "200-400", "400-")


def write_tsv(path: str, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols + [f"sens_{s}" for s in STRATA_COLS]) + "\n")
        for r in rows:
            vals = [fmt(r.get(c)) for c in cols]
            vals += [fmt(r["strata"].get(s, {}).get("sensitivity")) for s in STRATA_COLS]
            fh.write("\t".join(vals) + "\n")


def write_md(path: str, pairs: list[dict], rows: list[dict], sizes: list[int], a) -> None:
    with open(path, "w") as fh:
        fh.write(f"# KA/KS sweep (`baselines/kaks/run.py` {TOOL_VERSION}, windows {','.join(map(str, sizes))}, "
                 f"alpha {a.alpha}{'' if a.no_bonferroni else ' / 6 per frame'}, min codons {a.min_codons})\n\n")
        fh.write("Pooled over the genes of each pair: strand-aware nucleotide metrics over whole windows "
                 "(`docs/benchmark.md` 4.1), then the same over aligned bases only. `aligned` is the fraction of "
                 "window bases with an informant base. Window counts: called / tested-not-called / identical / "
                 "saturated / too few codons. The last five columns are sensitivity by CDS segment length (bp).\n\n")
        for pair in pairs:
            prs = [r for r in rows if r["pair"] == pair["pair"]]
            if not prs:
                continue
            fh.write(f"## {pair['pair']}: {pair['reference']} {pair['track']} {pair['chrom']}, "
                     f"{prs[0]['genes']} genes ({pair.get('distance_note', '')})\n\n")
            fh.write("| informant | dist | W | genes | aligned | windows c/t/i/s/f | Sens | Prec | F1 | MCC | "
                     "aligned F1 | aligned MCC | <60 | 60-100 | 100-200 | 200-400 | 400+ |\n")
            fh.write("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
            for r in prs:
                st = "/".join(str(r[f"win_{s}"]) for s in wm.STATUSES)
                sv = [fmt(r["strata"].get(s, {}).get("sensitivity"), 2) for s in STRATA_COLS]
                fh.write(f"| {r['informant']} | {fmt(r['distance'], 2)} | {r['W']} | {r['genes']} | {fmt(r['aligned_fraction'], 2)} | {st} | "
                         f"{fmt(r['sensitivity'], 2)} | {fmt(r['precision'], 2)} | {fmt(r['f1'], 2)} | {fmt(r['mcc'], 2)} | "
                         f"{fmt(r['aligned_f1'], 2)} | {fmt(r['aligned_mcc'], 2)} | " + " | ".join(sv) + " |\n")
            fh.write("\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=os.path.join(HERE, "pairs.tsv"))
    ap.add_argument("--pairs-filter", default="", help="comma-separated pair names to run (default: all)")
    ap.add_argument("--cache", default=os.environ.get("KAKS_CACHE", "/tmp/kaks-cache"), help="samples and windows (outside the repository)")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--windows", default="90,150,300,600")
    ap.add_argument("--step", type=int, default=None, help="step in bases (default: window / 3)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--no-bonferroni", action="store_true")
    ap.add_argument("--min-codons", type=int, default=15)
    ap.add_argument("--max-ref-stops", type=int, default=None)
    ap.add_argument("--transcript-types", default="benchmark")
    ap.add_argument("--seed", type=int, default=20260909, help="gene sample seed (sample_genes.py)")
    ap.add_argument("--limit-genes", type=int, default=None, help="use only the first N genes of every sample")
    ap.add_argument("--offline", action="store_true", help="never fetch; use the cache only")
    ap.add_argument("--pause", type=float, default=0.5, help="seconds between UCSC requests inside fetch_window.py")
    a = ap.parse_args()

    pairs = read_pairs(a.pairs)
    if a.pairs_filter:
        keep = {s.strip() for s in a.pairs_filter.split(",")}
        pairs = [p for p in pairs if p["pair"] in keep]
    sizes = [int(x) for x in a.windows.split(",")]
    os.makedirs(a.out, exist_ok=True)
    manifest = {"tool": "baselines/kaks/run.py", "version": TOOL_VERSION, "kaks_version": kaks.__doc__.split("\n")[0][:80],
                "windows_version": wm.TOOL_VERSION, "argv": sys.argv[1:], "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "pairs": pairs, "samples": {}, "inputs": {}, "windows_missing": {}}
    all_gene_rows, all_rows = [], []
    t0 = time.time()
    for pair in pairs:
        try:
            genes, header = sample_genes(pair, a.cache, a.offline, a.seed)
        except (RuntimeError, SystemExit) as e:
            print(f"{pair['pair']}: skipped, {e}", file=sys.stderr)
            manifest["windows_missing"][pair["pair"]] = ["(no sample)"]
            continue
        if a.limit_genes:
            genes = genes[: a.limit_genes]
        manifest["samples"][pair["pair"]] = {"header": header, "genes": [g["gene"] for g in genes]}
        stems = {}
        for g in genes:
            stem = fetch_window(pair, g, a.cache, a.offline, a.pause)
            if stem:
                stems[g["gene"]] = stem
            else:
                manifest["windows_missing"].setdefault(pair["pair"], []).append(g["gene"])
        print(f"{pair['pair']}: {len(stems)} / {len(genes)} windows, {time.time() - t0:.0f} s", file=sys.stderr)
        gene_rows, rows = run_pair(pair, genes, stems, sizes, a, manifest)
        all_gene_rows += gene_rows
        all_rows += rows
        print(f"{pair['pair']}: {len(rows)} pooled rows, {time.time() - t0:.0f} s", file=sys.stderr)
    cols = ["pair", "reference", "informant", "distance", "W", "step", "genes", "window_bp", "cds_bases", "aligned_fraction"] + \
           [f"win_{s}" for s in wm.STATUSES] + list(METRIC_KEYS) + ["sensitivity", "precision", "f1", "mcc"] + \
           [f"aligned_{k}" for k in METRIC_KEYS] + ["aligned_sensitivity", "aligned_precision", "aligned_f1", "aligned_mcc"]
    gcols = ["pair", "reference", "informant", "distance", "gene", "chrom", "start", "end", "window_bp", "cds_bases", "W", "step",
             "aligned_fraction"] + [f"win_{s}" for s in wm.STATUSES] + list(METRIC_KEYS) + \
            ["sensitivity", "precision", "f1", "mcc"] + [f"aligned_{k}" for k in METRIC_KEYS] + \
            ["aligned_sensitivity", "aligned_precision", "aligned_f1", "aligned_mcc"]
    write_tsv(os.path.join(a.out, "summary.tsv"), all_rows, cols)
    write_tsv(os.path.join(a.out, "per_gene.tsv"), all_gene_rows, gcols)
    write_md(os.path.join(a.out, "summary.md"), pairs, all_rows, sizes, a)
    manifest["finished"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(a.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    with open(os.path.join(a.out, "summary.md")) as fh:
        sys.stdout.write(fh.read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
