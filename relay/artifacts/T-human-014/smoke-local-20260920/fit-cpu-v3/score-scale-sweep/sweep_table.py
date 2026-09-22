#!/usr/bin/env python3
"""The section 3.7 curve: one row per (scale, chromosome) from the sweep dirs.
Usage: sweep_table.py DIR [DIR...]"""
import json, sys
from pathlib import Path

NAMES = {"NC_001133.9": "chr I (yeast)", "NC_003283.11": "chr V (nematode)"}
hdr = ("scale", "chrom", "cpu_s/Mb", "exact_tx", "exon_exF1", "nt_F1", "nt_sens", "nt_prec",
       "locus_F1", "fusion", "split", "donF1", "accF1", "GTAG_tp", "GTAG_fp", "GTAG_fn", "pred_tx")
print("\t".join(hdr))
for d in map(Path, sys.argv[1:]):
    for m in sorted(d.glob("measure_*.json")):
        j = json.load(open(m))
        n = m.name[len("measure_"):-5]
        s = json.load(open(d / f"score_{n}.json"))
        sc = j.get("splice_pwm_scale")
        sp, ex, nt, lo, tr = s["splice"], s["exon"]["all"], s["nucleotide"], s["locus"], s["transcript"]
        gt = sp["by_dinucleotide"].get("GT-AG", {"tp": 0, "fp": 0, "fn": 0})
        f1 = lambda p: p["f1"]
        print("\t".join(str(v) for v in (
            "none" if sc is None else f"{sc:g}", NAMES.get(j["seqid"], j["seqid"]),
            f"{j['cpu_s_per_mb']:.2f}", tr["tp"], f"{ex['f1']:.5f}", f"{nt['f1']:.5f}",
            f"{nt['sensitivity']:.5f}", f"{nt['precision']:.5f}", f"{f1(lo):.5f}",
            lo["fusion"], lo["split"], f"{f1(sp['donor']):.5f}", f"{f1(sp['acceptor']):.5f}",
            gt["tp"], gt["fp"], gt["fn"], s["predicted_transcripts"])))
