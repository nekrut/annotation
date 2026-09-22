# fit-cpu-v3: longer bounded local CPU fit with a decayed learning rate (2026-09-22T11:07Z to 15:27Z, exit 0)

The longer fit proposed in a-pilot section 3.3 after v2 sharply reduced
the fusions (4 remain on chr V, 2 on chr I, from 238 and 2 under v1;
[v2 score table](../fit-cpu-v2/score-final/README.md)) but left splicing
and precision missing. The GPU grant requested in `lenin-0083` is still
open, so this is the local CPU fallback named in the same section, on one
core as before.

Unchanged from `../fit-cpu-v2/`: species, dev reservations
(*S. cerevisiae* NC_001133.9, *C. elegans* NC_003283.11), loader
(`context: 512`, `background_windows: 800`, `background_length: 2048`,
`max_window: 12288`), loss kernel, decoder, optimizer (Adam, weight decay
0, grad clip 1.0), batch size 8, base `lr` 3e-4, seed 0, `dev_windows_max`
256, launcher (core 2, `OMP_NUM_THREADS=1`, `/usr/bin/time -v`),
leakage check before the run and MD5-gated sources.

Three config changes, all declared in `config.json` before the run
(engels-0105 / stalin-0110 P3: `eval_every` was previously described as
neutral bookkeeping, and it is not — see below):

- `steps: 3000` (v2: 1,500) — 24,000 planned draws against 19,975 train
  windows, rather than 12,000. The ratio 24,000/19,975 = 1.2015 is
  **draw-equivalent epochs, not pool coverage**: the sampler draws
  uniformly *with replacement*, so replaying the declared seed and batch
  shape gives 14,027 distinct windows among the 24,000 planned draws
  (70.22 %), leaving 5,948 train windows never attempted; v2's 12,000
  draws reached 9,035 distinct windows (45.23 %). These are planned
  index-exposure counts, not observed accepted-window counts. This is the
  "longer fit" of the proposed revision.
- `lr_schedule: "cosine"`, `warmup_steps: 150`, `lr_min_factor: 0.05`
  (v2: constant 3e-4) — linear ramp over the first 150 steps, then a half
  cosine from 3e-4 down to 1.5e-5 at step 3,000 (`model.a.train.lr_at`,
  code `6c7aac1`, 10 unit tests after the engels-0105 endpoint fix). v1
  and v2 both selected step 900 of
  1,500 and then oscillated: v2 dev NLL went 45.3 (step 900) → 74.1
  (1,100) → 48.6 (1,300) → 58.7 (1,400) → 50.7 (1,500) at a constant step
  size, so more steps at 3e-4 would very likely have kept bouncing rather
  than converged. The decay is therefore part of "longer", not an
  independent accuracy intervention; it is still a changed variable in
  its own right.
- `eval_every: 150` (v2: 100). Not only reporting cadence: only evaluated
  steps can become `best.pt`, so v3 selects from 20 candidate checkpoints
  on a 150-step grid where v2 selected from 15 on a 100-step grid, and
  over the shared first 1,500 steps the two grids offer 10 and 15
  candidates at different steps. The fixed dev subsample keeps individual
  NLL values comparable; it does not remove this checkpoint-selection
  difference.

A v2 → v3 difference is therefore attributable to the combined
step-budget, learning-rate-schedule and checkpoint-cadence revision.

**Correction (engels-0105 / stalin-0110).** An earlier version of this
file claimed `eval_every: 150` keeps the evaluation share of fit CPU at
v2's ~8 %. It does not. From the
[v2 manifest](../fit-cpu-v2/run_manifest.json), evaluation cost
`E` = 659.18535381 CPU-s against inclusive fit `F` = 7,938.915226944
CPU-s, i.e. `E/F` = 8.303 %. Holding per-step and per-evaluation unit
costs fixed while doubling fitting updates and taking 20/15 as many
evaluations projects `(20/15·E)/(2·(F−E)+20/15·E)` = **5.693 %**. That is
a projection under an unverified constant-unit-cost assumption; v3's
actual share is reported here from its own manifest after it exits.

The dev subsample is drawn by the same seed from the same v2 dev pool, so
v2 and v3 raw dev NLLs *are* comparable to each other (unlike v1 vs v2);
`score-final/` on chr I and chr V remains the reported comparison.

## Result (run finished 15:27:57Z, `train exit=0`, `/usr/bin/time` `Exit status: 0`)

3,000 steps in **15,562.6 CPU-s inclusive of evaluations = 4.32 CPU-h**
on one core (`/usr/bin/time` user+system 15,588.05 s = 4.33 CPU-h,
wall 4:19:59), 5.188 CPU-s/step against v2's 5.293, peak RSS 2.11 GiB
(2,209,192 kB). Loading cost 24.16 CPU-s, as in v2. Both the launcher's
recorded `train exit=0` and the GNU-time `Exit status: 0` are present;
per stalin-0113 the launcher's own status alone would not establish this.

- **Actual evaluation share 868.75 / 15,562.64 = 5.582 %**, against the
  5.693 % projected above from v2's 8.303 % under the constant-unit-cost
  assumption. The projection was close but not exact, as flagged.
- **Best checkpoint is the last step**: step 3,000 of 3,000, dev NLL
  **28.640** on the fixed 256-window subsample, against v2's 45.335 at
  step 900 of 1,500. `best.pt` sha256 `8a32c93e…`.
- Dev NLL trace: 176.0 (150), 66.6, 46.2, 48.6, 46.4, 52.0, 50.4, 50.2,
  41.6, 42.7, 38.0 (1,650), 48.3, 49.0, **32.4 (2,100), 32.4, 29.9,
  29.92, 31.1, 30.1, 28.64 (3,000)**, at `lr` 3.000e-04 → 1.500e-05. The
  ~10-NLL oscillation of the first half disappears once the rate falls
  below ~1e-4; the last seven evaluations sit in a 28.6–32.4 band. The
  monotone tail and the best point landing on the final step mean the
  budget, not convergence, ended this fit — a still longer schedule is
  not ruled out by the trace.
- Accepted-draw composition: 24,000 accepted of 24,000 attempted,
  71,218,348 sampled bases = **40.0 % CDS + 29.6 % intron + 30.3 %
  intergenic** (19,194,637 context bases, 955 background draws), matching
  v2's 39.9 / 29.6 / 30.6 at twice the volume. The supervision balance
  is unchanged; only the number of updates and the rate schedule are.

Dev NLL is not the acceptance criterion. The chromosome scores are in
[`score-final/`](score-final/README.md): nucleotide F1 0.825 → 0.891 on
chr I and 0.538 → 0.598 on chr V, locus F1 0.558 → 0.786 and 0.413 →
0.603, cost unchanged at 12.38 and 8.88 CPU-s/Mb — but exact-transcript
sensitivity on chr V falls 0.015 → 0.010, chr V fusions rise 4 → 312,
and donor/acceptor F1 reach only 0.070/0.081. A does not meet the
accuracy target under v3, and no positive CPU allowance for B follows.

`score-final/` held the scoring inputs, committed before fit completion
and final checkpoint selection, and before any v3 chromosome scoring:
v2's launcher with the v2 -> v3 paths changed, v2's `declaration.yaml`
with only the `model:` line changed, and v2's seqid lists. They were then
run unchanged after the fit exited 0, and that directory now holds the
measurements. The earlier "before the checkpoint exists" wording was
wrong — the training loop writes `best.pt` on every improving development
evaluation, so interim checkpoints predate the commit (engels-0108,
stalin-0113); what is fixed in advance is the scoring configuration, not
the absence of a model. The launcher has since gained an explicit
completion gate and workload failure propagation, with every workload
argument still identical to v2's; see
[`score-final/README.md`](score-final/README.md).

## Files

- `run.sh` — the launcher, as run.
- `config.json` — the declared `TrainConfig`.
- `leakage_check.out`, `source_md5.txt` — pre-run gates (leakage exit 0,
  0 violations).
- `run.log` — start/end stamps, commit, core, CPU model.
- `score-final/` — the pre-registered scoring inputs and the chr I /
  chr V measurements run under `best.pt`.
- `train.out`, `train_time.txt`, `run_manifest.json`, `best_pt.sha256` —
  the completed run: per-evaluation trace, `/usr/bin/time -v` report,
  manifest (hyperparameters with the per-step rate, history, timing,
  composition, source digest) and the checkpoint checksum.
