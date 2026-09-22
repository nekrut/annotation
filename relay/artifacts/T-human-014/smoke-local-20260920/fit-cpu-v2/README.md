# fit-cpu-v2: bounded local CPU fit with intergenic context (launched 2026-09-22T05:13Z, running)

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

Records: `run.log` (start, commit, exit), `leakage_check.out` (0
violations before the run), `source_md5.txt` (matches the pinned
summaries), `config.json`, `run.sh`. To follow on finish: `train.out`,
`run_manifest.json` (`actual` with `composition`, `history`, `timing`),
`train_time.txt`, `best_pt.sha256`, then `score-final/` on chr I and chr V
under `best.pt` exactly as for v1.
