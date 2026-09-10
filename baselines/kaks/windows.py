#!/usr/bin/env python3
"""Sliding-window KA/KS coding classifier on one fetched alignment window,
scored against the reference annotation at the nucleotide level.

Standard library only (Python 3.11).  Reads the ``<stem>.*`` files that
``scripts/data/fetch_window.py`` writes and reuses the MAF painting of
``scripts/data/cut_windows.py`` (imported by path, so a checkout is enough).

What it does
------------
For one informant species, the reference sequence is paired base by base
with the informant row of the alignment (first block in file order wins,
``cut_windows.paint_informants``).  Windows of ``W`` reference bases at
step ``S`` are read in all three frames on both strands; in each frame the
window's codon pairs are those whose three reference bases are aligned to
A, C, G or T of the informant with no informant insertion between them.
``kaks.kaks`` gives dN, dS and the z-test; a window is *called* coding in
a frame when dN/dS is below one and dN is significantly below dS
(``kaks.is_coding``), and the window's call is the most significant frame.
Reference stop codons in a frame are counted (``ref_stops``) and, with
``--max-ref-stops``, can veto the frame, since a coding frame has none
inside an exon.

The per-base prediction is the union of called windows: a reference base
is predicted coding on a strand when any called window in a frame of that
strand covers it.  Truth is the union of CDS segments of the window's
coding transcripts per strand, which is ``docs/benchmark.md`` section 4.1
(a CDS on the other strand is a false positive and a false negative).
Sensitivity, precision, F1 and MCC are computed over (base, strand)
pairs; the same four are also reported restricted to the bases at which
the informant is aligned, and sensitivity is stratified by the length of
the CDS segment a base belongs to, which is where the 2002 paper says the
test fails (short exons).  Windows are also tallied by status
(``called``, ``tested`` but not called, ``identical``, ``saturated``,
``too_few_codons``), because a window the test cannot decide is not a
negative.

Coordinates are those of the fetched window (0-based, half-open); called
segments are written as GFF3 on the window's chromosome so that a
whole-chromosome run could be scored by ``benchmark/score.py``.

Examples::

    python3 baselines/kaks/windows.py --stem /tmp/win/Adh/Adh_124 \
        --informants droSim2,droYak3,droAna3,droPse3 --windows 90,150,300 \
        --out /tmp/kaks/Adh
    python3 baselines/kaks/windows.py --self-test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "scripts", "data"))
sys.path.insert(0, HERE)
import cut_windows as cw  # noqa: E402
import pairwise_codons as pc  # noqa: E402
import kaks  # noqa: E402

TOOL_VERSION = "0.1"
CODE = "ACGT"
STATUSES = ("called", "tested", "identical", "saturated", "too_few_codons")
LENGTH_STRATA = ((0, 60), (60, 100), (100, 200), (200, 400), (400, 10 ** 9))


def revcomp(s: str) -> str:
    return s.translate(cw.COMP)[::-1]


def informant_string(inf: bytearray, ins: bytearray) -> tuple[str, list[bool]]:
    """Informant letters over the window ('-' gap, '.' unaligned, 'N' other)
    and, per position, whether an informant insertion follows it."""
    letters = []
    for v in inf:
        letters.append(CODE[v] if v < 4 else ("-" if v == cw.INF_GAP else ("." if v == cw.INF_UNALIGNED else "N")))
    return "".join(letters), [x > 0 for x in ins]


def frame_pairs(ref: str, inf: str, ins: list[bool], lo: int, hi: int, frame: int, strand: str):
    """Codon pairs of window [lo, hi) in one frame on one strand, with the
    count of reference stops in the frame and of codons dropped."""
    pairs = []
    ref_stops = 0
    dropped = 0
    p = lo + frame
    while p + 3 <= hi:
        r = ref[p:p + 3]
        q = inf[p:p + 3]
        if strand == "-":
            r, q = revcomp(r), revcomp(q)
        if r in kaks.tables(1).stops or "N" in r:
            ref_stops += 1 if r in kaks.tables(1).stops else 0
        ok = all(ch in "ACGT" for ch in q) and not ins[p] and not ins[p + 1] and all(ch in "ACGT" for ch in r)
        if ok:
            pairs.append((r, q))
        else:
            dropped += 1
        p += 3
    return pairs, ref_stops, dropped


def classify_window(ref: str, inf: str, ins: list[bool], lo: int, hi: int, alpha: float, min_codons: int,
                    max_ref_stops: int | None, table=1, bonferroni: bool = True) -> dict:
    """Best frame of one window and the status list of all six frames.
    With ``bonferroni`` (the default) each frame is tested at ``alpha / 6``,
    because a window is called if *any* of six frames passes."""
    best = None
    statuses = []
    a = alpha / 6.0 if bonferroni else alpha
    for strand in "+-":
        for frame in range(3):
            pairs, stops, dropped = frame_pairs(ref, inf, ins, lo, hi, frame, strand)
            r = kaks.kaks(pairs, table=table, min_codons=min_codons)
            r.update({"strand": strand, "frame": frame, "ref_stops": stops, "codons_dropped": dropped})
            vetoed = max_ref_stops is not None and stops > max_ref_stops
            r["called"] = kaks.is_coding(r, a) and not vetoed
            statuses.append(r["status"])
            key = (0 if r["called"] else 1, r["p"] if r["p"] is not None else 2.0, -r["codons"])
            if best is None or key < best[0]:
                best = (key, r)
    r = best[1]
    if r["called"]:
        status = "called"
    elif "tested" in statuses:
        status = "tested"
    elif "identical" in statuses:
        status = "identical"
    elif "saturated" in statuses:
        status = "saturated"
    else:
        status = "too_few_codons"
    return {"lo": lo, "hi": hi, "status": status, "strand": r["strand"], "frame": r["frame"],
            "codons": r["codons"], "ref_stops": r["ref_stops"], "dS": r["dS"], "dN": r["dN"],
            "ratio": r["ratio"], "z": r["z"], "p": r["p"], "called": r["called"]}


def truth_arrays(transcripts: list[dict], start: int, end: int) -> tuple[bytearray, bytearray, list[int], list[int]]:
    """Per-base CDS truth on + and - strands (union over coding transcripts)
    and, per base, the length of the shortest CDS segment covering it."""
    n = end - start
    plus, minus = bytearray(n), bytearray(n)
    seg_plus, seg_minus = [0] * n, [0] * n
    for t in transcripts:
        arr, seg = (plus, seg_plus) if t.get("strand", "+") == "+" else (minus, seg_minus)
        for a, b in t.get("cds", []):
            L = b - a
            for i in range(max(a, start), min(b, end)):
                arr[i - start] = 1
                if seg[i - start] == 0 or L < seg[i - start]:
                    seg[i - start] = L
    return plus, minus, seg_plus, seg_minus


def confusion(pred: bytearray, truth: bytearray, mask=None) -> dict:
    tp = fp = fn = tn = 0
    for i in range(len(pred)):
        if mask is not None and not mask[i]:
            continue
        p, t = pred[i], truth[i]
        if p and t:
            tp += 1
        elif p:
            fp += 1
        elif t:
            fn += 1
        else:
            tn += 1
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn}


def metrics(c: dict) -> dict:
    tp, fp, fn, tn = c["tp"], c["fp"], c["fn"], c["tn"]
    sens = tp / (tp + fn) if tp + fn else None
    prec = tp / (tp + fp) if tp + fp else None
    f1 = (2 * sens * prec / (sens + prec)) if sens and prec else (0.0 if sens is not None and prec is not None else None)
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / den) if den else None
    return {**c, "sensitivity": sens, "precision": prec, "f1": f1, "mcc": mcc}


def add(c: dict, d: dict) -> dict:
    return {k: c[k] + d[k] for k in ("tp", "fp", "fn", "tn")}


def merged_segments(pred: bytearray, start: int) -> list[tuple[int, int]]:
    out = []
    i = 0
    n = len(pred)
    while i < n:
        if pred[i]:
            j = i
            while j < n and pred[j]:
                j += 1
            out.append((start + i, start + j))
            i = j
        else:
            i += 1
    return out


def run_one(ref: str, inf: str, ins: list[bool], start: int, W: int, S: int, alpha: float, min_codons: int,
            max_ref_stops: int | None, truth, table=1, bonferroni: bool = True) -> tuple[dict, list[dict]]:
    n = len(ref)
    plus, minus, seg_plus, seg_minus = truth
    pred_plus, pred_minus = bytearray(n), bytearray(n)
    windows = []
    counts = {s: 0 for s in STATUSES}
    lo = 0
    while lo + W <= n:
        w = classify_window(ref, inf, ins, lo, lo + W, alpha, min_codons, max_ref_stops, table, bonferroni)
        counts[w["status"]] += 1
        if w["called"]:
            arr = pred_plus if w["strand"] == "+" else pred_minus
            for i in range(lo, lo + W):
                arr[i] = 1
        windows.append(w)
        lo += S
    aligned = [c in "ACGT" for c in inf]
    c_plus, c_minus = confusion(pred_plus, plus), confusion(pred_minus, minus)
    both = add(c_plus, c_minus)
    a_both = add(confusion(pred_plus, plus, aligned), confusion(pred_minus, minus, aligned))
    strata = {}
    for lo_len, hi_len in LENGTH_STRATA:
        tp = fn = 0
        for arr, pred, seg in ((plus, pred_plus, seg_plus), (minus, pred_minus, seg_minus)):
            for i in range(n):
                if arr[i] and lo_len <= seg[i] < hi_len:
                    if pred[i]:
                        tp += 1
                    else:
                        fn += 1
        strata[f"{lo_len}-{hi_len if hi_len < 10 ** 9 else ''}"] = {"tp": tp, "fn": fn,
                                                                    "sensitivity": (tp / (tp + fn)) if tp + fn else None}
    frame_ok = frame_wrong = 0
    # of called windows entirely inside one CDS segment on the truth strand: is the frame right?
    for w in windows:
        if not w["called"]:
            continue
        arr = plus if w["strand"] == "+" else minus
        if all(arr[i] for i in range(w["lo"], w["hi"])):
            frame_ok += 1
    result = {"window": W, "step": S, "alpha": alpha, "bonferroni_over_frames": bonferroni,
              "min_codons": min_codons, "max_ref_stops": max_ref_stops,
              "windows": len(windows), "window_status": counts,
              "called_windows_inside_cds_on_truth_strand": frame_ok,
              "nucleotide": metrics(both),
              "nucleotide_plus": metrics(c_plus), "nucleotide_minus": metrics(c_minus),
              "nucleotide_aligned_only": metrics(a_both),
              "aligned_fraction": sum(aligned) / n if n else None,
              "sensitivity_by_cds_segment_length": strata,
              "predicted_segments": {"+": merged_segments(pred_plus, start), "-": merged_segments(pred_minus, start)}}
    return result, windows


def load(stem: str, informant: str, drop: set[str], transcript_types: str = "benchmark"):
    manifest, ref, chrom, start, end, seq, blocks, trees = pc.load_window(stem)
    rows, dist = pc.informant_rows(manifest, ref, blocks, trees)
    with open(stem + ".annotation.json") as fh:
        ann = json.load(fh)
    transcripts = ann.get("transcripts", ann if isinstance(ann, list) else [])
    transcripts, _ = cw.select_transcripts(transcripts, transcript_types)
    transcripts = [t for t in transcripts if t.get("cds")]
    return manifest, ref, chrom, start, end, seq, blocks, rows, dist, transcripts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stem", help="fetched window stem, e.g. /tmp/win/Adh/Adh_124")
    ap.add_argument("--informants", help="comma-separated informant names (default: every row)")
    ap.add_argument("--drop-species", default="", help="comma-separated rows to ignore")
    ap.add_argument("--windows", default="90,150,300", help="window sizes in reference bases")
    ap.add_argument("--step", type=int, default=None, help="step in bases (default: window / 3, at least 3)")
    ap.add_argument("--alpha", type=float, default=0.05, help="family-wise level for the six frames of a window")
    ap.add_argument("--no-bonferroni", action="store_true", help="test every frame at --alpha instead of --alpha / 6")
    ap.add_argument("--min-codons", type=int, default=15, help="fewest aligned codon pairs for a frame to be tested")
    ap.add_argument("--max-ref-stops", type=int, default=None,
                    help="veto a frame with more reference stop codons than this (default: no veto)")
    ap.add_argument("--table", type=int, default=1, help="genetic code: 1 standard, 6 ciliate")
    ap.add_argument("--transcript-types", default="benchmark")
    ap.add_argument("--out", help="output directory; <name>.<informant>.json, .windows.tsv and .gff3 per run")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not (a.stem and a.out):
        ap.error("--stem and --out are required")
    drop = {s.strip() for s in a.drop_species.split(",") if s.strip()}
    manifest, ref, chrom, start, end, seq, blocks, rows, dist, transcripts = load(a.stem, None, drop, a.transcript_types)
    rows = [r for r in rows if r not in drop]
    wanted = [s.strip() for s in a.informants.split(",")] if a.informants else rows
    missing = [r for r in wanted if r not in rows]
    if missing:
        raise SystemExit(f"informants not in the window's rows: {missing}; rows are {rows[:20]}...")
    inf_all, ins_all, block_stats = cw.paint_informants(blocks, ref, start, end, wanted)
    truth = truth_arrays(transcripts, start, end)
    os.makedirs(a.out, exist_ok=True)
    name = os.path.basename(a.stem)
    checks = {}
    for ext in (".maf", ".fa", ".annotation.json", ".manifest.json"):
        with open(a.stem + ext, "rb") as fh:
            checks[ext] = hashlib.sha256(fh.read()).hexdigest()
    sizes = [int(x) for x in a.windows.split(",")]
    summary_rows = []
    for k, informant in enumerate(wanted):
        inf, ins = informant_string(inf_all[k], ins_all[k])
        runs = []
        with open(os.path.join(a.out, f"{name}.{informant}.windows.tsv"), "w") as tsv, \
                open(os.path.join(a.out, f"{name}.{informant}.gff3"), "w") as gff:
            tsv.write("window\tstep\tlo\thi\tstatus\tstrand\tframe\tcodons\tref_stops\tdS\tdN\tratio\tz\tp\tcalled\n")
            gff.write("##gff-version 3\n")
            for W in sizes:
                S = a.step or max(3, W // 3)
                res, windows = run_one(seq, inf, ins, start, W, S, a.alpha, a.min_codons, a.max_ref_stops, truth,
                                       a.table, not a.no_bonferroni)
                runs.append(res)
                for w in windows:
                    tsv.write("\t".join(str(w[x]) if w[x] is not None else "" for x in
                                        ("lo", "hi")).join([f"{W}\t{S}\t", ""]) +
                              "\t" + "\t".join("" if w[x] is None else (f"{w[x]:.4g}" if isinstance(w[x], float) else str(w[x]))
                                               for x in ("status", "strand", "frame", "codons", "ref_stops", "dS", "dN", "ratio", "z", "p", "called")) + "\n")
                for strand in "+-":
                    for s0, e0 in res["predicted_segments"][strand]:
                        gff.write(f"{chrom}\tkaks\tCDS\t{s0 + 1}\t{e0}\t.\t{strand}\t.\tID=kaks_{informant}_W{W}_{s0};window={W}\n")
                m = res["nucleotide"]
                summary_rows.append((informant, dist.get(informant), W, S, res["window_status"], m, res["aligned_fraction"]))
        out = {"tool": "baselines/kaks/windows.py", "version": TOOL_VERSION, "stem": a.stem, "name": name,
               "reference": ref, "chrom": chrom, "start": start, "end": end, "informant": informant,
               "distance": dist.get(informant), "alignment_track": manifest.get("alignment_track") or manifest.get("track"),
               "transcripts": [t.get("id") for t in transcripts], "cds_bases": {"+": sum(truth[0]), "-": sum(truth[1])},
               "input_sha256": checks, "block_stats": {k: block_stats[k] for k in ("overlapping_blocks", "reference_minus_strand_blocks_flipped", "duplicate_rows_discarded")},
               "runs": runs}
        with open(os.path.join(a.out, f"{name}.{informant}.json"), "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
    print(f"{name}: {len(transcripts)} coding transcripts, CDS bases + {sum(truth[0])} / - {sum(truth[1])}, {len(seq)} bp")
    print("informant\tdist\tW\tstep\taligned\tcalled\ttested\tident\tsat\tfew\tsens\tprec\tF1\tMCC")
    for informant, d, W, S, st, m, af in summary_rows:
        fmt = lambda v: "" if v is None else f"{v:.3f}"
        print(f"{informant}\t{fmt(d)}\t{W}\t{S}\t{fmt(af)}\t{st['called']}\t{st['tested']}\t{st['identical']}\t{st['saturated']}\t{st['too_few_codons']}\t"
              f"{fmt(m['sensitivity'])}\t{fmt(m['precision'])}\t{fmt(m['f1'])}\t{fmt(m['mcc'])}")
    return 0


# ------------------------------------------------------------------ tests

def _synthetic_pair(rng, n_codons: int, purifying: bool, rate: float):
    """Reference codons and an informant copy: purifying keeps the amino acid
    (third-position synonymous changes only where possible), neutral changes
    any base; stops are avoided."""
    T = kaks.tables(1)
    ref, alt = [], []
    for _ in range(n_codons):
        c = rng.choice([x for x in T.sense if x not in ("ATG", "TGG")])
        ref.append(c)
        if purifying:
            syn = [c[:2] + b for b in "ACGT" if c[:2] + b != c and T.code[c[:2] + b] == T.code[c]]
            alt.append(rng.choice(syn) if syn and rng.random() < rate else c)
        else:
            while True:
                m = "".join(rng.choice("ACGT") if rng.random() < rate / 3 else b for b in c)
                if m not in T.stops:
                    break
            alt.append(m)
    return "".join(ref), "".join(alt)


def self_test() -> int:
    import random
    rng = random.Random(20260910)
    checks = 0
    # a 1,200 bp reference: 300 bp neutral, a 450 bp plus-strand CDS, 150 bp neutral,
    # a 240 bp minus-strand CDS, 60 bp neutral; the informant differs at 30% of
    # synonymous sites in CDS and at ~25% of bases elsewhere
    parts = []
    r0, a0 = _synthetic_pair(rng, 100, False, 0.25)
    r1, a1 = _synthetic_pair(rng, 150, True, 0.3)
    r2, a2 = _synthetic_pair(rng, 50, False, 0.25)
    r3, a3 = _synthetic_pair(rng, 80, True, 0.3)
    r4, a4 = _synthetic_pair(rng, 20, False, 0.25)
    ref = r0 + r1 + r2 + revcomp(r3) + r4
    inf = a0 + a1 + a2 + revcomp(a3) + a4
    assert len(ref) == 1200 and len(inf) == 1200
    ins = [False] * 1200
    start = 1000
    transcripts = [{"id": "plus", "strand": "+", "start": 1300, "end": 1750, "exons": [[1300, 1750]], "cds": [[1300, 1750]]},
                   {"id": "minus", "strand": "-", "start": 1900, "end": 2140, "exons": [[1900, 2140]], "cds": [[1900, 2140]]}]
    truth = truth_arrays(transcripts, start, start + 1200)
    assert sum(truth[0]) == 450 and sum(truth[1]) == 240 and truth[2][300] == 450 and truth[3][900] == 240
    checks += 1
    res, windows = run_one(ref, inf, ins, start, 90, 30, 0.05, 15, None, truth)
    m = res["nucleotide"]
    print("synthetic 90/30:", {k: round(v, 3) if isinstance(v, float) else v for k, v in m.items()}, res["window_status"])
    assert m["sensitivity"] > 0.9 and m["precision"] > 0.8, m
    # called windows inside the plus CDS are mostly on the plus strand in frame 0, inside the
    # minus CDS mostly on the minus strand: with changes only at third positions the
    # antisense frame whose third positions coincide (the "shadow" frame) also passes, and
    # the two are separated only by the few changes that are nonsynonymous in the shadow frame
    plus_in = [w for w in windows if w["called"] and w["lo"] >= 300 and w["hi"] <= 750]
    minus_in = [w for w in windows if w["called"] and w["lo"] >= 900 and w["hi"] <= 1140]
    assert plus_in and minus_in
    assert sum(w["strand"] == "+" and w["frame"] == 0 for w in plus_in) * 2 > len(plus_in), [(w["strand"], w["frame"]) for w in plus_in]
    assert sum(w["strand"] == "-" for w in minus_in) * 2 > len(minus_in), [(w["strand"], w["frame"]) for w in minus_in]
    checks += 2
    # strand matters: predicting the minus CDS on the plus strand would be wrong
    # (the plus-strand false positives are the shadow frame of the minus CDS; a window is 90 bp)
    assert res["nucleotide_minus"]["tp"] > 200 and res["nucleotide_plus"]["fp"] <= 180, (res["nucleotide_minus"], res["nucleotide_plus"])
    checks += 1
    # unaligned informant: nothing can be called there, and the aligned-only metrics exclude it
    inf2 = inf[:600] + "." * 600
    res2, _ = run_one(ref, inf2, ins, start, 90, 30, 0.05, 15, None, truth)
    assert res2["nucleotide_minus"]["tp"] == 0 and res2["nucleotide_aligned_only"]["fn"] < res2["nucleotide"]["fn"]
    assert abs(res2["aligned_fraction"] - 0.5) < 1e-9
    checks += 2
    # an insertion inside a codon drops that codon pair in every frame that spans it
    ins3 = [False] * 1200
    ins3[400] = True
    p_a, _, d_a = frame_pairs(ref, inf, ins, 300, 390 + 30, 0, "+")
    p_b, _, d_b = frame_pairs(ref, inf, ins3, 300, 390 + 30, 0, "+")
    assert len(p_b) == len(p_a) - 1 and d_b == d_a + 1
    checks += 1
    # identical sequences are 'identical', not negatives
    res4, w4 = run_one(ref, ref, ins, start, 90, 30, 0.05, 15, None, truth)
    assert res4["window_status"]["identical"] == res4["windows"] and res4["nucleotide"]["tp"] == 0
    checks += 1
    # stop veto: a frame with reference stops can be vetoed
    res5, w5 = run_one(ref, inf, ins, start, 90, 30, 0.05, 15, 0, truth)
    assert all(w["ref_stops"] == 0 for w in w5 if w["called"])
    checks += 1
    # merged segments and strata
    segs = merged_segments(bytearray([0, 1, 1, 0, 1]), 10)
    assert segs == [(11, 13), (14, 15)]
    assert set(res["sensitivity_by_cds_segment_length"]) == {"0-60", "60-100", "100-200", "200-400", "400-"}
    checks += 2
    # metrics arithmetic
    mm = metrics({"tp": 8, "fp": 2, "fn": 2, "tn": 88})
    assert abs(mm["f1"] - 0.8) < 1e-9 and abs(mm["mcc"] - (8 * 88 - 4) / math.sqrt(10 * 10 * 90 * 90)) < 1e-9
    checks += 1
    print(f"self-test passed ({checks} checks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
