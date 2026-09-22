# fit-cpu-v2: bounded local CPU fit with `U` context (2026-09-22T05:13Z–07:26Z, finished; exit 0)

The controlled revision proposed in a-pilot section 3.3 after fit v1 missed
the accuracy target. Same species, dev reservations, optimizer, step count,
batch, seed, dev-sampling settings (256-window limit, seed 0) and launcher
as `../fit-cpu-v1/`; the only config changes are the two data settings
below, declared in `config.json`. They change the training *and* the
development pool (stalin-0104): dev identities 4,893 → 4,951 (4,839
common), and the seeded 256-window dev subset shares only 37 identities
with v1's, the retained chains also gaining context. Raw v1/v2 dev NLLs
are therefore on different windows; the common comparison is
`score-final/` on chr I and chr V.

- `context: 512` — every clean chain window is widened by up to 512
  bases per side, clipped at the nearest other gene's *CDS span* and the
  sequence ends (`model.a.dataset.iter_windows`, code `61a9480`), so the
  added bases are supervised as `U` by the unchanged support mask; no loss
  or decoder change. The added context is coding-span-free, not
  necessarily annotated intergenic: UTRs, ncRNA genes and CDS-free
  pseudogenes can fall inside it (engels-0100: 48,213 of 3,965,426 added
  yeast bases and 3,964,153 of 12,049,332 added *C. elegans* bases overlap
  raw `gene`/`pseudogene` spans, summed over windows). `U` counts below
  are the support label's share, not strictly intergenic annotation.
- `background_windows: 800` per species (v1: 400), 2,048-base gene-free
  tiles. *S. cerevisiae* still yields only 148 such tiles, *C. elegans* 800.

Loader probe on the same sources (train pool after the dev split; the
manifest's `actual.composition` records what the draws actually sampled):

| context | train windows | mean length | `U` share of pool bases |
|---:|---:|---:|---:|
| 0 (v1) | 19,708 | 2,202 | 3.0 % |
| 256 | 19,673 | 2,654 | 20.2 % |
| **512 (v2)** | 19,975 incl. 773 background | 2,970 | **30.3 %** |
| 1,024 | 19,572 | 3,405 | 39.1 % |

v2 split: 19,975 train (773 background) / 4,951 dev (175 background); 948
background in all 24,926 windows. Train pool 59,330,766 bases =
23,891,418 CDS + 17,457,446 intron + 17,981,902 `U` (30.3 %); these are
inventory counts, the manifest's `actual.composition` will be the accepted
draws over the full run (not only through `best_step`).

512 was chosen because the mean window then matches the 3,072-base encoder
core and the `U` share approaches the genomes' own (~28 % of
*S. cerevisiae*, ~45 % of *C. elegans* is intergenic), at ~1.35× the
per-step cost of v1. Windows over `max_window` 12,288 after widening are
skipped (456 vs 381 in v1).

## Run of record

`run.log`: start 05:13:31Z, code `61a9480`, core 2, exit 0 at 07:26:22Z.
`leakage_check.out` 0 violations before the run; `source_md5.txt` matches
the pinned summaries. `train_time.txt` (`/usr/bin/time -v`): **user
7,663.2 s + sys 301.1 s = 7,964 CPU-s = 2.21 CPU-h** on one core, wall
2:12:50, 99 % of one core, peak RSS **2.09 GiB** (2,194,924 kB). Manifest
timing: load 24.1 s, fit 7,939 s of which the 15 evaluations 659 s;
**5.3 s per step** of 8 windows (v1: 4.0 s), 1.37× v1's 5,833 CPU-s.

`train.out` / manifest `history` (train NLL per window; dev NLL on the 256
selected v2 dev windows, steps 100–1,500):
train 181 / 90 / 33 / 62 / 152 / 61 / 14 / 32 / 18 / 79 / 88 / 36 / 41 /
93 / 15; dev 284.9 / 68.8 / 56.0 / 90.9 / 52.0 / 76.3 / 54.6 / 45.6 /
**45.3** / 48.9 / 74.1 / 56.7 / 48.6 / 58.7 / 50.7. `best_step` is **900**
again (v1's was also 900), `best_dev_nll` 45.34; `best.pt` sha256
`99f773ed…2a8cf6` (`best_pt.sha256`, 1.8 MB, host-local). Both curves are
noisier than v1's (windows are longer and one third of the pool is `U`);
the dev value is not comparable to v1's 38.0 because only 37 of the 256
windows are shared (above).

`run_manifest.json` `actual`: 12,000 attempted = 12,000 accepted draws,
**sampled 35,305,664 bases = 14,070,834 CDS (39.9 %) + 10,439,449 intron
(29.6 %) + 10,795,381 `U` (30.6 %)**; the `U` bases are 9,561,661 added
context + 490 background draws × 2,048 = 1,003,520 + 11,510 chain draws ×
20 flank = 230,200. This is the accepted-draw composition over the full
run and matches stalin-0104's seeded replay (30.58 % `U` at 12,000
attempted draws) to the base. The manifest does not record the prefix
through `best_step`; by the same replay it is the first 7,200 draws. v1
by comparison sampled 26,259,772 bases of which 785,532 `U` (3.0 %). The
v2 checkpoint therefore saw ten times v1's `U` share and 13.7× its `U`
base count with the same optimizer, steps and batch.

`score-final/` scores chr I and chr V under `best.pt` exactly as for v1
(same `run_final.sh` with the v2 paths, `summarize.py` for the table).
