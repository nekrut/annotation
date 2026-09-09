---
id: T-human-010
title: KA/KS comparative baseline on the benchmark
status: open
owner: null
created_by: human
created: 2026-09-09T12:00:00Z
lease_until: null
depends_on: [T-human-007, T-human-008]
touches: [baselines/kaks/]
pr: null
---

## Goal

Reproduce the simplest comparative signal (Nekrutenko, Makova, Li 2002:
KA/KS of pairwise alignment windows) as a coding versus non-coding
classifier, run it on the benchmark, and establish the floor every later
model must beat.

## What to produce

- `baselines/kaks/`: code that takes alignment windows from the data
  inventory's fetch script, computes KA/KS in all three frames on both
  strands, and calls coding windows; a sweep over window size and
  evolutionary distance of the partner species.
- Results on at least four benchmark species pairs at different distances,
  reported with the benchmark's metrics.
- A short write-up: where the signal fails (short exons, fast-evolving
  genes, distant pairs, non-coding conserved elements) and what that says
  about the inputs the new model needs.

## Definition of done

Pull request merged after one `review`. Results reproducible from a single
command in a fresh clone.

## Log

- 2026-09-09 human: created.
