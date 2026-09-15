#!/usr/bin/env python3
"""Fold a Nextflow trace and a GNU time -v file into one measured.tsv-style row.

CPU-seconds are reported two ways: the sum over trace tasks of
realtime * %cpu / 100 (what the workers used) and user+sys from GNU time on
the driver (which includes children on a single node). Peak RSS is the largest
trace peak_rss, since Nextflow tasks run as separate processes and GNU time's
maximum resident set covers one process at a time. Standard library only.
"""
import argparse, csv, re, sys

UNITS = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3, "TB": 1024**4}


def parse_bytes(s):
    m = re.match(r"\s*([\d.]+)\s*([KMGT]?B)?", s or "")
    if not m or not m.group(1):
        return 0.0
    return float(m.group(1)) * UNITS.get(m.group(2) or "B", 1)


def parse_duration(s):
    """Nextflow durations: '1h 2m 3.4s', '45.2s', '3ms', '1d 2h'."""
    total = 0.0
    for num, unit in re.findall(r"([\d.]+)\s*(ms|[dhms])", s or ""):
        total += float(num) * {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 0.001}[unit]
    return total


def gnu_time(path):
    out = {}
    try:
        for line in open(path):
            if ":" not in line:
                continue
            k, v = line.rsplit(":", 1)
            out[k.strip()] = v.strip()
    except OSError:
        pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trace", required=True)
    ap.add_argument("--time", required=True, help="GNU time -v output")
    ap.add_argument("--wall-s", type=float, required=True)
    ap.add_argument("--genome-bp", type=int, required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--exit", type=int, default=0)
    a = ap.parse_args()

    tasks = ok = failed = 0
    cpu_s = 0.0
    peak_rss = 0.0
    try:
        for row in csv.DictReader(open(a.trace), delimiter="\t"):
            tasks += 1
            if (row.get("status") or "").upper() == "COMPLETED":
                ok += 1
            elif (row.get("status") or "").upper() in ("FAILED", "ABORTED"):
                failed += 1
            pct = float((row.get("%cpu") or "0").rstrip("%") or 0)
            cpu_s += parse_duration(row.get("realtime")) * pct / 100.0
            peak_rss = max(peak_rss, parse_bytes(row.get("peak_rss")))
    except OSError:
        pass
    t = gnu_time(a.time)
    user = float(t.get("User time (seconds)", 0) or 0)
    sys_ = float(t.get("System time (seconds)", 0) or 0)
    driver_rss_kb = float(t.get("Maximum resident set size (kbytes)", 0) or 0)
    mb = a.genome_bp / 1e6
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["species", "genome_mb", "exit", "tasks", "completed", "failed", "wall_s",
                "cpu_s_trace", "cpu_s_gnu_time", "cpu_s_per_mb_trace", "cpu_h_trace",
                "peak_task_rss_mb", "driver_rss_mb"])
    w.writerow([a.species, f"{mb:.4f}", a.exit, tasks, ok, failed, f"{a.wall_s:.1f}",
                f"{cpu_s:.1f}", f"{user + sys_:.1f}", f"{cpu_s / mb:.1f}" if mb else "",
                f"{cpu_s / 3600:.2f}", f"{peak_rss / 1024**2:.0f}", f"{driver_rss_kb / 1024:.0f}"])


if __name__ == "__main__":
    main()
