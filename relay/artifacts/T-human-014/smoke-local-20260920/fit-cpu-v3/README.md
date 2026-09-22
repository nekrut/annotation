# fit-cpu-v3: longer bounded local CPU fit with a decayed learning rate (started 2026-09-22T11:07Z, running)

The longer fit proposed in a-pilot section 3.3 after v2 removed the
fusions but left splicing and precision missing. The GPU grant requested
in `lenin-0083` is still open, so this is the local CPU fallback named in
the same section, on one core as before.

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

Projected cost at v2's measured 5.3 CPU-s/step: ~4.4 CPU-h on one core,
finishing around 15:30Z. Nothing is claimed from it until it exits 0 and
the manifest, history, timing and `best.pt` checksum are collected here.

## Files

- `run.sh` — the launcher, as run.
- `config.json` — the declared `TrainConfig`.
- `leakage_check.out`, `source_md5.txt` — pre-run gates (leakage exit 0,
  0 violations).
- `run.log` — start/end stamps, commit, core, CPU model.
- `train.out`, `train_time.txt`, `run_manifest.json`, `best_pt.sha256` —
  added when the run finishes.
