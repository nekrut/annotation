---
id: T-human-009
title: Cost baseline of existing tools and a compute budget for ours
status: review
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-10T15:54:19Z
depends_on: [T-human-006]
touches: [docs/cost-baseline.md]
pr: https://github.com/nekrut/annotation/pull/20
---

## Goal

Put numbers on "current tools are slow and expensive", and set the budget
the new model must fit in.

## What to produce

- `docs/cost-baseline.md`: for EGAPx, Tiberius, BRAKER3, AUGUSTUS, Helixer,
  and GeneMark-ETP at minimum: documented or measured wall clock, CPU or GPU
  hours, peak memory, and failure rate per genome, normalized per megabase.
  Use the Galaxy EGAPx analysis in nekrut/scalingPaper (1,409 jobs, 290
  genomes, about 416,000 CPU-hours, 45% failure, 24% CPU efficiency) as one
  primary source and cite others.
- Where you can, run one tool on one small benchmark genome on a laptop and
  record the actual numbers. If a tool cannot run within the laptop
  constraint, say so; that is itself a finding. Measuring EGAPx or Tiberius
  end to end will need the cluster: post an `alert` with an estimate, as
  the charter describes, and run there once granted.
- A target budget: what "fast" must mean for the new model in seconds per
  megabase on one CPU core and on one consumer GPU, and the parameter and
  memory ceiling, with justification.

## Definition of done

Pull request merged after one `review`. Every number has a source.

## Log

- 2026-09-09 human: created.
- 2026-09-10 marx: claimed.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: claimed once T-human-006 went `done` (the draft had been
  built on `work/T-human-009-marx` over the previous ticks, see marx-0026,
  0031, 0032). This tick: merged `main` into the branch (clean), rewrote the
  status line, cited the merged synthesis instead of the open PR, added the
  Helixer consumer-GPU memory note from `benchmark/validation/README.md`
  (default batch 32 exhausts 16 GB at the vertebrate window; the validation
  declarations for Tiberius and Helixer record no time or memory, checked),
  posted the cluster `alert` 20260910T115307Z-marx-0033 for EGAPx on
  *C. intestinalis* (32 CPU, 256 GB, 70 to 100 CPU-h) and Tiberius on the
  panel (one GPU, about 3 GPU-h), and opened PR #20
  (https://github.com/nekrut/annotation/pull/20, head 86a9091, 7 files,
  +532). Status -> review. Definition of done needs one `review` from
  another agent and a merge. Next: answer review comments; run the two
  cluster rows only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: lenin's `review` (lenin-0033) received, verdict accept
  after six fixes; all six numbers re-derived here and confirmed (panel
  13.959 Gb; covered Tiberius species 10.602 Gb; `cds_fraction_pct` sums
  145,477 human transcripts). Pushed 01985d9 on `work/T-human-009-marx`
  (PR #20): panel 14 Gb and 1.9 h; Tiberius row lists the eleven species,
  10.60 Gb, 5.6 GPU-h at 1.9 GPU-s/Mb; §5.3 conclusion 2 now cites a
  measured genomic CDS union (new `docs/cost-baseline/cds_union.py` and
  `cds_union.tsv`: human 1.34%, mouse 1.38%, maize 1.96%, yeast 72.5%,
  md5-checked RefSeq GFFs) instead of the isoform-summed column; one-core
  peak stated with two FMA ports and conclusion 1 shown to hold at 9e10;
  GPU ceiling normalized to tensor throughput (4090 165 vs A100 312 TFLOPS
  BF16, sources cited); CPU ceiling bound to the §3.1 runner with the
  AUGUSTUS ratio as the portable form. §6 records lenin's two suggested rows
  (gigabase AUGUSTUS: 35 to 100 CPU-h, outside one tick; timed Tiberius and
  Helixer reruns: need the RTX 5080 box) as not measured, with reasons. The
  corrected cluster estimate is restated in marx's note this tick because
  marx-0033 cannot be edited. Status stays `review`; the merge by `human`
  is what remains. Next: any further review comments; cluster rows only on
  a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: lenin's re-review (lenin-0034) received, verdict accept,
  one non-blocking finding: the human CDS union counted GRCh38.p14's alt
  and patch scaffolds, which are not in the panel's `genome_bp`. Landed it
  now rather than deferring, since the PR is unmerged and the fix is small:
  `cds_union.py` fetches the NCBI assembly report, keeps only
  primary-assembly nuclear seqids (roles assembled-molecule, unlocalized-
  and unplaced-scaffold; unit not non-nuclear) and checks their lengths sum
  to `genome_bp` (exact for all four species); `--all` restores the old
  count. Re-derived: human primary union 36,126,376 bp = 1.17% (lenin's
  36,137,717 includes chrM; the 11,341 bp difference is chrM CDS); the
  unfiltered run still gives 41,405,486, so primary plus excluded
  (5,279,110) reproduces the old table exactly. Isoform inflation on the
  primary assembly is 7.6x (273.37/36.13), not the 8.3x in lenin-0034.
  Other rows: yeast 72.35% (chrM dropped), mouse 1.38%, maize 1.95%.
  `cds_union.tsv` regenerated with two new columns (`sequences`,
  `excluded_cds_union_bp`); §3.4 explains the filter; §5.3 conclusion 2
  reads 1.2% human, 8.8% isoform sum, 7.6x, candidate stage 1.7 to 4.3x.
  Pushed as one commit on `work/T-human-009-marx` (PR #20). Status stays
  `review`; nothing else in lenin-0034 asked for a change; the merge by
  `human` is what remains. Next: further review comments if any; cluster
  rows only on a `decision`.
