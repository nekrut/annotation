# score-final/: **prepared, not run** (2026-09-22T14:07Z)

Nothing here is a measurement yet. The fit-cpu-v3 run is still in flight
(step 1,950 of 3,000 at 14:05Z), so no `best.pt` is final and no scoring
has been executed. This directory holds only the inputs, committed before
the run so that the comparison against
[v2's scores](../../fit-cpu-v2/score-final/README.md) is fixed in advance
rather than chosen after seeing the v3 checkpoint:

- `run_final.sh` — byte-for-byte
  [v2's launcher](../../fit-cpu-v2/score-final/run_final.sh) with only
  `fit-cpu-v2` -> `fit-cpu-v3` and `score-final-v2` -> `score-final-v3`
  substituted (verified by substituting back and diffing). Same two
  chromosomes, same `measure --profile chromosome --window 12288
  --overlap 4096 --segments 19 --dtype float64`, same core 2, one thread,
  `/usr/bin/time -v`, then `benchmark/score.py` with `--seqids`,
  `--genome` and this `--declaration`.
- `declaration.yaml` — identical to v2's except the `model:` line, which
  names the v3 schedule (3,000 steps, cosine 3e-4 -> 1.5e-5 after a
  150-step warmup, `eval_every` 150) and code commit `6c7aac1`.
- `seqids_chrI.txt`, `seqids_chrV.txt` — copied from v2 unchanged.

The launcher pins core 2, the core the fit occupies, so it runs only
after the fit exits 0; the cost rows here and v2's are then measured on
the same idle core.

When the run happens this file is replaced by the usual cost and accuracy
tables with `run_final.out`, the `measure_*` / `score_*` JSON, timings,
predicted GFF3s, `sha256.txt` and `summarize.py`, and the v2 columns
alongside. Until then, no cost, accuracy or CPU-allowance claim follows
from this directory.
