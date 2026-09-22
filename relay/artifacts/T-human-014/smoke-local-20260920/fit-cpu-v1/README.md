# fit-cpu-v1: bounded local CPU fit of candidate A (finished 2026-09-22T03:59Z)

Run of record for the first fitted checkpoint (a-pilot section 3.3). Code
`4c43819` (PR #38), config `config.json`, launcher `run.sh`
(leakage check → source MD5 → `taskset -c 2 /usr/bin/time -v python -m
model.a.train train`), one thread on one core of an Intel Core Ultra 9 285K.

- Sources: *S. cerevisiae* (dev chr I `NC_001133.9`) + *C. elegans* (dev
  chr V `NC_003283.11`); `leakage_check.out` 0 violations; `source_md5.txt`
  matches the committed summaries (loader-gated).
- 24,601 windows: 19,708 train / 4,893 dev (548 gene-free 2,048-base
  background tiles); 1,500 Adam steps, batch 8, lr 3e-4, eval every 100
  steps on 256 seeded dev windows.
- Finished: `train.out` (15 evaluations), `run_manifest.json` (`actual`:
  `best_step` **900**, `best_dev_nll` **38.04**, full `history`, `timing`),
  `train_time.txt` (`/usr/bin/time -v`: user 5,631.9 s + sys 201.4 s =
  **5,833 s = 1.62 CPU-h**, wall 1:37:17, 99 % of one core, peak RSS
  **2.14 GiB**), `run.log` (start/end, commit, exit 0). Timing inside the
  manifest: load 20.7 s, fit 5,811 s of which evaluations 467 s.
- Dev NLL on the 256 selected windows by step: 129 / 60 / 42 / 45 / 66 / 58
  / 49 / 46 / **38** / 48 / 51 / 52 / 70 / 59 / 44 (steps 100–1,500); train
  NLL per window bounces 16–71. Both curves are noisy at batch 8 and the
  run covers ~0.6 of one pass (12,000 draws, 26.26 Mb sampled).
- `best.pt` (step 900, 1.8 MB) stays on the host; sha256 in
  `best_pt.sha256`
  (`181adad7c9da846ee9c6faa4d85a3ab5aea7896cd55c97e2cd0bbfac4c4cf88b`).

## score-final/: development-chromosome scores under the finished `best.pt`

See `score-final/README.md` (filled in the same tick the fit finished).

## score-dry/: scoring-pipeline dry run (interim checkpoint, no accuracy claim)

While the fit ran, the accuracy pipeline was exercised end to end on the
**interim** `best.pt` (step 300, copied as `best_interim.pt`, sha256 in
`sha256.txt`) on *S. cerevisiae* chr I, core 3: `measure --profile
chromosome --segments 19 --overlap 4096 --window 12288 --dtype float64
--gff-out` → `benchmark/score.py --seqids seqids_chrI.txt --genome ...
--declaration declaration.yaml`. Both exit 0 (`run_dry.sh`).

- Cost row (unfitted-vs-fitted encoder cost is checkpoint-independent):
  230,218 bases, 76 tiles / 38 windows, **12.79 CPU-s/Mb** (encoder 1.30 s,
  decode 1.60 s, preprocess 0.04 s, io 0.009 s, output 0.003 s), peak RSS
  0.96 GiB; `score.py` 0.17 CPU-s, 33 MB.
- Interim step-300 numbers on chr I (94 reference transcripts; these are a
  pipeline check on an unfinished fit, not the reported accuracy):
  nucleotide sensitivity 0.963 / precision 0.456 (F1 0.619, MCC 0.439);
  locus 84 TP / 578 FP / 10 FN (662 predicted vs 94 reference, 7 fusions);
  transcript F1 0.069; 320 predicted introns against 3 reference, 316 of
  them non-GT-AG. The decoder at step 300 over-predicts CDS on both strands
  and opens non-canonical introns freely; whether the finished fit closes
  that is the next tick's question.
- `declaration.yaml` is the section 3.3 submission declaration for this
  run (training species listed, both dev chromosomes are *seen* species;
  no pretraining, protein DB, alignment, informants or RNA-seq).
