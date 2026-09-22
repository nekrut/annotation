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

Two config changes, both declared in `config.json` before the run:

- `steps: 3000` (v2: 1,500) — 24,000 planned draws, about 1.2 passes over
  the 19,975 train windows rather than 0.6. This is the "longer fit" of
  the proposed revision.
- `lr_schedule: "cosine"`, `warmup_steps: 150`, `lr_min_factor: 0.05`
  (v2: constant 3e-4) — linear ramp over the first 150 steps, then a half
  cosine from 3e-4 down to 1.5e-5 at step 3,000 (`model.a.train.lr_at`,
  code `6c7aac1`, 8 unit tests). v1 and v2 both selected step 900 of
  1,500 and then oscillated: v2 dev NLL went 45.3 (step 900) → 74.1
  (1,100) → 48.6 (1,300) → 58.7 (1,400) → 50.7 (1,500) at a constant step
  size, so more steps at 3e-4 would very likely have kept bouncing rather
  than converged. The decay is therefore part of "longer", not an
  independent accuracy intervention; it is also a second changed variable,
  so a v2→v3 difference is attributable to the pair, not to step count
  alone.

`eval_every: 150` (v2: 100) keeps the 20 evaluations of the run at the
same ~8 % share of fit CPU as v2's 15. The dev subsample is drawn by the
same seed from the same v2 dev pool, so v2 and v3 raw dev NLLs *are*
comparable to each other (unlike v1 vs v2); `score-final/` on chr I and
chr V remains the reported comparison.

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
