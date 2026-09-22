# score-final/: fit-cpu-v3 scored on chr I and chr V, 2026-09-22T15:32Z

The pre-registered launcher was run unchanged after the fit exited 0. All
four workloads exit 0 (`run_final.out`), the gate passed
(`gate: fit-cpu-v3 complete, checkpoint present`), and `summarize.py`
reproduces every table below from the committed JSON and GFF3.

Checkpoint: `best.pt` sha256 `8a32c93e…`, step 3,000 of 3,000, dev NLL
28.640 on the fixed 256-window subsample (v2: step 900, 45.335).
Measurement commit `85d534c`; the fit ran at `6c7aac1`. The only code
difference between them inside `model/` and `benchmark/` is the cosine
config validation added for engels-0105 plus its tests — no inference,
decoder or scorer path changed (`git diff 6c7aac1 85d534c -- model/
benchmark/`).

## Cost (one core, `taskset -c 2`, `/usr/bin/time -v`, float64, 19 segments)

| chromosome | bases | CPU-s | CPU-s/Mb v3 | v2 | encoder | decode | peak RSS |
|---|---|---|---|---|---|---|---|
| S. cerevisiae chr I (NC_001133.9) | 230,218 | 2.85 | **12.38** | 12.54 | 1.29 | 1.52 | 0.96 GiB |
| C. elegans chr V (NC_003283.11) | 20,924,180 | 185.80 | **8.88** | 8.88 | 77.51 | 105.48 | 2.05 GiB |

Cost is unchanged: same architecture, same tiling, only the weights
differ. The CPU verdict of section 5 therefore stands as measured before,
and **no positive CPU allowance for B follows from this run**. Decode is
still the larger stage on chr V (105.48 of 185.80 CPU-s).

## Accuracy, v3 against v2

| metric | chr I v2 | chr I v3 | chr V v2 | chr V v3 |
|---|---|---|---|---|
| nucleotide F1 | 0.825 | **0.891** | 0.538 | **0.598** |
| nucleotide MCC | 0.746 | **0.841** | 0.468 | **0.549** |
| nt sensitivity / precision | 0.931 / 0.742 | 0.912 / 0.870 | 0.531 / 0.546 | 0.533 / 0.680 |
| predicted chains | 214 | 107 | 13,410 | 6,156 |
| median predicted CDS span (b) | 346.5 | 1,149 | 294 | 621 |
| locus TP / FP / FN | 86 / 128 / 8 | 79 / 28 / 15 | 3,804 / 9,606 / 1,191 | 3,361 / 2,794 / 1,634 |
| locus F1 | 0.558 | **0.786** | 0.413 | **0.603** |
| fusion / split | 2 / 1 | 3 / 0 | 4 / 1,827 | **312** / 910 |
| exact transcripts (sens) | 69 (0.734) | 64 (0.681) | 101 (0.015) | 70 (0.010) |
| exon F1 (exact) | 0.407 | 0.582 | 0.011 | 0.030 |
| donor F1 | 0 | 0 | 0.014 | 0.070 |
| acceptor F1 | 0 | 0 | 0.025 | 0.081 |
| GT-AG introns TP / FP / FN | 0 / 1 / 3 | 0 / 0 / 3 | 68 / 97 / 22,552 | 433 / 470 / 22,187 |
| predicted introns outside GT-AG/GC-AG | 27 | 16 | 2,118 | 3,998 |
| of which `other` class (tp + fp) | 27 | 16 | 2,114 | 3,987 |
| predicted introns (median length) | 28 (21 b) | 16 (81.5 b) | 2,290 (34 b) | 4,931 (615 b) |

What the longer fit with the decayed rate changed, and what it did not:

- **Precision, not sensitivity.** Nucleotide sensitivity is flat on both
  chromosomes (0.931 → 0.912, 0.531 → 0.533); precision carries the whole
  improvement (0.742 → 0.870, 0.546 → 0.680). Chain count halves on both
  and the median predicted CDS span roughly doubles, so the gain comes
  from emitting fewer, longer chains rather than from finding more coding
  sequence.
- **Intron length moved from too short to too long, and the splice
  sites are still wrong.** The v2 pile-up at the 20 b floor is gone
  (introns of exactly 20 b: 538 of 2,290 = 23.49% → 88 of 4,931 = 1.78%)
  but the median predicted intron on chr V goes from 34 b to 615 b
  against a reference q50 of 56 b, and 47.60% of predictions now exceed
  the reference q90 of 695 b (v2: 0.83%); reference decile cuts
  `[45, 46, 48, 51, 56, 95, 199, 370, 695]` from `splice` in the same
  score JSON. Length is **not** fixed — it overshoots (engels-0110).
  Correct GT-AG introns rise 68 → 433, but donor F1 is 0.070 and
  acceptor F1 0.081: 22,187 reference GT-AG introns are still missed and
  3,998 predicted introns fall outside GT-AG/GC-AG (3,987 in the `other`
  class — 3,986 fp plus 1 tp — and 11 AT-AC), more in absolute number
  than v2's 2,118 because v3 predicts more introns overall. Splice
  placement and duration calibration both remain unsolved.
- **A new fusion regime on chr V.** Splits fall 1,827 → 910, but fusions
  rise 4 → 312 and exact transcripts fall 101 → 70 (sens 0.015 → 0.010).
  Longer chains merge neighbouring genes. On chr I, fusions go 2 → 3 and
  exact transcripts 69 → 64 while nucleotide F1 rises, the same trade in
  miniature. The locus-level F1 gain (0.413 → 0.603 on chr V, 0.558 →
  0.786 on chr I) is real but is bought partly with fusion.
- **The accuracy target is still missed.** Exact-transcript sensitivity
  on chr V is 0.010. Nothing here accepts T-human-014.

Both v2 and v3 locus FP counts mix two kinds: predictions with no
same-strand CDS overlap at all, and overlapping predictions that failed
to match a reference locus. The v2 split (chr I 127 / 1, chr V 6,376 /
3,230) was reproduced for stalin-0106; the equivalent split for v3 has
not been computed and is not claimed here.

## Host conditions

One core (core 2), `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, on the same
24-core Intel Core Ultra 9 285K used for v1 and v2, immediately after the
fit released that core; load average 0.97 at launch. An idle core is an
operator precondition, not something the launcher's gate can establish
(engels-0108); the cost columns above are comparable with v2's on that
basis.

## Files

`run_final.sh` (pre-registered, with the gate and failure propagation),
`run_final.out`, `declaration.yaml`, the two seqid lists, `measure_*.json`
and `score_*.json`, `/usr/bin/time -v` reports for all four workloads,
the two predicted GFF3s (chr V gzipped), `sha256.txt` over all of them,
and `summarize.py` (copied unchanged from v2's directory), which prints
every number above from this directory.
