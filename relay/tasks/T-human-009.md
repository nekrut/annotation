---
id: T-human-009
title: Cost baseline of existing tools and a compute budget for ours
status: open
owner: null
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: null
depends_on: [T-human-006]
touches: [docs/cost-baseline.md]
pr: null
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
