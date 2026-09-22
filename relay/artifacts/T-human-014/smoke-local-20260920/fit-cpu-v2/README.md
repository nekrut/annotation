# fit-cpu-v2: bounded local CPU fit with intergenic context (launched 2026-09-22T05:13Z, running)

The controlled revision proposed in a-pilot section 3.3 after fit v1 missed
the accuracy target. Same species, dev reservations, optimizer, step count,
batch, seed, dev subsample and launcher as `../fit-cpu-v1/`; the only
changes are in the training data, both declared in `config.json`:

- `context: 512` — every clean chain window is widened by up to 512
  annotated-intergenic bases per side, clipped at the nearest other gene's
  span and the sequence ends (`model.a.dataset.iter_windows`, code
  `61a9480`), so the added bases are supervised as intergenic `U` by the
  unchanged support mask; no loss or decoder change.
- `background_windows: 800` per species (v1: 400), 2,048-base gene-free
  tiles. *S. cerevisiae* still yields only 148 such tiles, *C. elegans* 800.

Loader probe on the same sources (train pool after the dev split; the
manifest's `actual.composition` records what the draws actually sampled):

| context | train windows | mean length | intergenic share of pool bases |
|---:|---:|---:|---:|
| 0 (v1) | 19,708 | 2,202 | 3.0 % |
| 256 | 19,673 | 2,654 | 20.2 % |
| **512 (v2)** | 19,975 incl. 948 background | 2,970 | **30.3 %** |
| 1,024 | 19,572 | 3,405 | 39.1 % |

512 was chosen because the mean window then matches the 3,072-base encoder
core and the intergenic share approaches the genomes' own (~28 % of
*S. cerevisiae*, ~45 % of *C. elegans* is intergenic), at ~1.35× the
per-step cost of v1. Windows over `max_window` 12,288 after widening are
skipped (456 vs 381 in v1).

Records: `run.log` (start, commit, exit), `leakage_check.out` (0
violations before the run), `source_md5.txt` (matches the pinned
summaries), `config.json`, `run.sh`. To follow on finish: `train.out`,
`run_manifest.json` (`actual` with `composition`, `history`, `timing`),
`train_time.txt`, `best_pt.sha256`, then `score-final/` on chr I and chr V
under `best.pt` exactly as for v1.
