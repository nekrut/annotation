---
id: T-human-008
title: Inventory usable data: alignments, conservation, expression, annotations
status: open
owner: null
created_by: human
created: 2026-09-09T12:00:00Z
lease_until: null
depends_on: []
touches: [docs/data-sources.md, scripts/data/]
pr: null
---

## Goal

Know exactly what data the new model can be trained on, how to get it, and
how to cut it into training examples.

## What to produce

- `docs/data-sources.md`: for UCSC (multiz and Cactus whole-genome
  alignments and their trees, phyloP and phastCons, RNA-seq, regulation,
  variation tracks), Ensembl Compara, NCBI RefSeq, Zoonomia, and any others
  you find: what assemblies are covered, formats (MAF, HAL, bigWig, GFF3),
  sizes, access methods (rsync, bigDataUrl, APIs), licences and rate limits.
- A section on non-coding alignment: how well introns, UTRs, and
  intergenic sequence are aligned in each source and at what evolutionary
  distances alignment breaks down.
- `scripts/data/fetch_window.py`: given an assembly, a locus, and a
  flanking size, fetch the multiple alignment block with its Newick tree
  and the reference annotation for that window, in a format a model can
  consume. Demonstrate on one human locus and one non-mammal locus.
- A short assessment of which sources are suitable for training versus
  held-out evaluation given the benchmark's leakage rules.

## Definition of done

Pull request merged after one `review`. The window script runs in a fresh
clone for both demonstrated loci without committing data over 5 MB.

## Log

- 2026-09-09 human: created.
