#!/usr/bin/env python3
"""Fold the cluster jobs' return files into measured.tsv rows and a Markdown table.

Reads what `gagarin` copies into `relay/artifacts/T-gagarin-001/<label>/`
(decision 20260915T143249Z-human-0016): the Tiberius `summary.tsv` written by
`tiberius_panel.sbatch` (one row per species), the EGAPx `summary.tsv` written
by `summarize_trace.py`, the scorer JSON next to each (`<tool>-<species>.json`
from `benchmark/score.py`), and `jobs.tsv` / `usage.txt` for provenance.

Rewrites `docs/cost-baseline/measured.tsv` with two added columns, `tool`
(first) and `peak_gpu_mem_mb` (after `peak_rss_mb`); existing AUGUSTUS rows
keep their values with `tool=augustus` and an empty GPU cell. Rows are keyed
by (tool, species) and replaced on rerun, so the fold is idempotent. Rows
whose job exit status is non-zero are reported on stderr and left out unless
`--include-failed` is given. `--markdown` prints the same rows in the layout
of section 3.2 of docs/cost-baseline.md. Standard library only.

    python3 docs/cost-baseline/cluster/fold_returns.py            # rewrite measured.tsv
    python3 docs/cost-baseline/cluster/fold_returns.py --markdown # print table rows
    python3 docs/cost-baseline/cluster/fold_returns.py --dry-run  # show, do not write
"""
import argparse
import csv
import json
import os
import sys

COLUMNS = ["tool", "species", "genome_mb", "processes", "wall_s", "cpu_user_s",
           "cpu_s_per_mb", "peak_rss_mb", "peak_gpu_mem_mb", "genes",
           "nucleotide_f1", "locus_f1", "transcript_f1", "source"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_usage(label_dir):
    """usage.txt is one `key=value ...` line then free text; return the keys."""
    out = {}
    path = os.path.join(label_dir, "usage.txt")
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        first = fh.readline()
    for tok in first.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def read_jobs(art):
    """jobs.tsv: label -> row (job_id, submitted, cpus, mem_gb, gpus, ...)."""
    path = os.path.join(art, "jobs.tsv")
    if not os.path.exists(path):
        return {}
    return {r["label"]: r for r in read_tsv(path)}


def scores(path):
    """(nucleotide_f1, locus_f1, transcript_f1) from a score.py JSON, or blanks."""
    if not os.path.exists(path):
        return "", "", ""
    with open(path) as fh:
        d = json.load(fh)
    f = lambda k: f"{d[k]['f1']:.5f}" if k in d and "f1" in d[k] else ""
    return f("nucleotide"), f("locus"), f("transcript")


def fnum(x, digits=1):
    try:
        return f"{float(x):.{digits}f}"
    except (TypeError, ValueError):
        return ""


def provenance(label, job, usage, extra=""):
    parts = ["T-gagarin-001"]
    if job or usage.get("job_id"):
        # usage.txt names the job that actually produced the files; jobs.tsv
        # may still carry an earlier, failed submission of the same label
        parts.append(f"Slurm job {usage.get('job_id') or job.get('job_id', '?')}")
        job = job or {}
        sub = job.get("submitted", "")
        if sub:
            parts.append(f"submitted {sub[:10]}")
        parts.append(f"{job.get('cpus', '?')} CPUs x {job.get('mem_gb', '?')} GB"
                     + (f" x {job['gpus']} GPU" if job.get("gpus", "0") not in ("", "0") else ""))
    if usage.get("node"):
        parts.append(usage["node"])
    parts.append(label)
    if extra:
        parts.append(extra)
    return ", ".join(parts)


def tiberius_rows(art, label, job, include_failed):
    d = os.path.join(art, label)
    summ = os.path.join(d, "summary.tsv")
    if not os.path.exists(summ):
        return []
    usage = read_usage(d)
    rows = []
    for r in read_tsv(summ):
        sp = r["species"]
        exit_ = r.get("exit", "")
        mb = float(r["genome_mb"])
        user = float(r.get("cpu_user_s") or 0)
        if exit_ not in ("0", ""):
            print(f"fold: {label} {sp} exit={exit_}; not measured"
                  + (" (kept: --include-failed)" if include_failed else ""), file=sys.stderr)
            if not include_failed:
                continue
        n, l, t = scores(os.path.join(d, f"tiberius-{sp}.json"))
        rows.append({
            "tool": "tiberius", "species": sp, "genome_mb": f"{mb:.4f}", "processes": "1",
            "wall_s": fnum(r.get("wall_s")), "cpu_user_s": fnum(user),
            "cpu_s_per_mb": f"{user / mb:.1f}" if mb and exit_ in ("0", "") else "",
            "peak_rss_mb": r.get("peak_rss_mb", ""), "peak_gpu_mem_mb": r.get("peak_gpu_mem_mb", ""),
            "genes": r.get("genes", "") if exit_ in ("0", "") else "",
            "nucleotide_f1": n, "locus_f1": l, "transcript_f1": t,
            "source": provenance(label, job, usage,
                                 f"model {r.get('model_cfg', '?')}, batch {r.get('batch_size', '?')}"
                                 + (f", exit {exit_}" if exit_ not in ("0", "") else "")),
        })
    return rows


def egapx_rows(art, label, job, include_failed):
    d = os.path.join(art, label)
    summ = os.path.join(d, "summary.tsv")
    if not os.path.exists(summ):
        return []
    usage = read_usage(d)
    rows = []
    for r in read_tsv(summ):
        sp = r["species"]
        exit_ = r.get("exit", "")
        mb = float(r["genome_mb"])
        cpu = float(r.get("cpu_s_trace") or 0)
        if exit_ not in ("0", ""):
            print(f"fold: {label} {sp} exit={exit_}; not measured"
                  + (" (kept: --include-failed)" if include_failed else ""), file=sys.stderr)
            if not include_failed:
                continue
        n, l, t = scores(os.path.join(d, f"egapx-{sp}.json"))
        rows.append({
            "tool": "egapx", "species": sp, "genome_mb": f"{mb:.4f}",
            "processes": r.get("tasks", ""), "wall_s": fnum(r.get("wall_s")),
            "cpu_user_s": fnum(cpu),
            "cpu_s_per_mb": f"{cpu / mb:.1f}" if mb and exit_ in ("0", "") else "",
            "peak_rss_mb": r.get("peak_task_rss_mb", ""), "peak_gpu_mem_mb": "",
            "genes": "", "nucleotide_f1": n, "locus_f1": l, "transcript_f1": t,
            "source": provenance(label, job, usage,
                                 f"CPU-s from the Nextflow trace over {r.get('completed', '?')} completed"
                                 f" of {r.get('tasks', '?')} tasks (GNU time on the driver: "
                                 f"{fnum(r.get('cpu_s_gnu_time'))} s); peak RSS is the largest task"
                                 + (f", exit {exit_}" if exit_ not in ("0", "") else "")),
        })
    return rows


def markdown(rows):
    out = ["| tool | genome | size (Mb) | how run | wall clock (s) | CPU (user s) | **CPU-s / Mb** |"
           " peak RSS (MB) | peak GPU (MB) | genes | nucleotide F1 | locus F1 | transcript F1 | source |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        name = "*" + r["species"].replace("_", " ") + "*"
        how = ("1 process, whole genome, 1 GPU" if r["tool"] == "tiberius"
               else f"{r['processes']} Nextflow tasks" if r["tool"] == "egapx" else f"{r['processes']} processes")
        g = f"{int(float(r['genes'])):,}" if r["genes"] else ""
        cells = [r["tool"], name, f"{float(r['genome_mb']):.2f}", how,
                 fnum(r["wall_s"], 0), fnum(r["cpu_user_s"], 0),
                 f"**{fnum(r['cpu_s_per_mb'], 0)}**" if r["cpu_s_per_mb"] else "",
                 f"{int(float(r['peak_rss_mb'])):,}" if r["peak_rss_mb"] else "",
                 f"{int(float(r['peak_gpu_mem_mb'])):,}" if r["peak_gpu_mem_mb"] else "",
                 g, r["nucleotide_f1"], r["locus_f1"], r["transcript_f1"], r["source"]]
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(here, "..", "..", ".."))
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--artifacts", default=os.path.join(repo, "relay", "artifacts", "T-gagarin-001"))
    ap.add_argument("--measured", default=os.path.join(repo, "docs", "cost-baseline", "measured.tsv"))
    ap.add_argument("--out", default=None, help="write here instead of over --measured")
    ap.add_argument("--include-failed", action="store_true")
    ap.add_argument("--markdown", action="store_true", help="print the cluster rows as a section 3 table")
    ap.add_argument("--dry-run", action="store_true", help="print the folded TSV, write nothing")
    a = ap.parse_args()

    jobs = read_jobs(a.artifacts)
    new = []
    if os.path.isdir(a.artifacts):
        for label in sorted(os.listdir(a.artifacts)):
            job = jobs.get(label)
            if label.startswith("tiberius"):
                new += tiberius_rows(a.artifacts, label, job, a.include_failed)
            elif label.startswith("egapx"):
                new += egapx_rows(a.artifacts, label, job, a.include_failed)
    new.sort(key=lambda r: (r["tool"], float(r["genome_mb"])))

    existing = read_tsv(a.measured) if os.path.exists(a.measured) else []
    for r in existing:
        r.setdefault("tool", "augustus")
        r.setdefault("peak_gpu_mem_mb", "")
    keys = {(r["tool"], r["species"]) for r in new}
    kept = [r for r in existing if (r["tool"], r["species"]) not in keys]
    rows = kept + new
    if a.markdown:
        print(markdown(new))
        return
    if a.dry_run:
        w = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
        return
    out = a.out or a.measured
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"fold: {len(new)} cluster row(s) folded, {len(kept)} kept; wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
