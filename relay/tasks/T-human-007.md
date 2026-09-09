---
id: T-human-007
title: Design the benchmark: species panel, references, metrics, leakage rules
status: open
owner: null
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: null
depends_on: []
touches: [docs/benchmark.md, benchmark/]
pr: null
---

## Goal

A benchmark that a model can be judged on and that rewards cross-clade
generalization rather than memorization.

## What to produce

- `docs/benchmark.md`: a species panel of 12 to 20 genomes covering
  mammals, birds, fish, insects, nematodes, plants, fungi, and at least one
  protist, chosen to span intron length distributions (median and tail),
  genome size, GC content, and gene density. Include the numbers that
  justify each pick. State which reference annotation is truth for each
  (RefSeq, Ensembl, MANE, community curated) and its known weaknesses.
- Metrics: nucleotide, exon, and gene level sensitivity and specificity;
  splice-site accuracy; BUSCO completeness of predicted proteomes; wall
  clock and peak memory per megabase. Define them precisely.
- Held-out rules: which species are never used in training, how to detect
  leakage through alignments that include held-out species, and a minimum
  phylogenetic distance rule.
- `benchmark/`: a manifest (`panel.tsv`) with assembly accessions, annotation
  versions, and checksums, plus a fetch script that does not commit data.

## Definition of done

Pull request merged after one `review`. The fetch script runs for at least
one species in a fresh clone.

## Log

- 2026-09-09 human: created.
