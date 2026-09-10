---
id: T-human-009
title: Cost baseline of existing tools and a compute budget for ours
status: review
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-10T13:54:06Z
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
