---
id: T-human-010
title: KA/KS comparative baseline on the benchmark
status: done
owner: null
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: null
depends_on: [T-human-007, T-human-008]
touches: [baselines/kaks/]
pr: https://github.com/nekrut/annotation/pull/14
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
- 2026-09-10 marx: claimed.
- 2026-09-10 marx: -> in_progress (work/T-human-010-marx).
- 2026-09-10 marx (tick 1 on this task): claimed after decision 20260910T004401Z-human-0009 (T-human-007 and T-human-008 done). Also did the decision's Zoonomia follow-up as PR #10 (branch `work/T-human-008-zoonomia-marx`: Cactus 241/447 marked excluded in `docs/data-sources.md` tables, section 5 and section 7; the 6.3 duplicate-row measurement kept as a recorded finding; no script defaulted to a Cactus track, the docstring/README example changed). Started `baselines/kaks/` on `work/T-human-010-marx` (a9c003e): `kaks.py` (Nei-Gojobori pathway counting, Jukes-Cantor, one-sided z-test of dN < dS, standard and ciliate codes, 22 self-test checks) and `windows.py` (reads a `fetch_window.py` window, pairs the reference with one informant, six frames per window at alpha/6, calls projected to bases, strand-aware nucleotide metrics per section 4.1 plus aligned-only metrics, sensitivity by CDS segment length, undecidable windows counted apart: identical / saturated / too few codons; 12 checks). The original test's codeml likelihood-ratio form is replaced by the counting z-test because codeml is outside the stack; full text of the paper is not open (abstract only; publisher PDF is a page image). First live numbers on fly Adh (124-way) and human TP53 (100-way) in the README: chimp (0.01) and simulans (0.10) are too close (identical windows), yakuba (0.23) F1 0.77 at 300 bp, mouse (0.50) F1 0.50 with 30% of the window aligned, chicken/zebrafish 7%/4% aligned, ananassae and Anopheles absent from the track at Adh; 90 bp windows call nothing at alpha/6. Next: the single-command driver over a species-pair table with a seeded gene sample per pair (human-chimp/mouse/chicken, fly-yakuba/pseudoobscura, mouse-rat on mm39 35-way, worm-briggsae on ce11 135-way), the window-size by distance sweep, then the write-up and PR.
- 2026-09-10 marx: -> in_progress.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/14).
- 2026-09-10 marx (tick 2): -> review, PR #14 (`work/T-human-010-marx`, a43e51f). Added `run.py` (single command: seeded gene sample per pair via `sample_genes.py`, windows fetched into a cache via `fetch_window.py`, sweep over informants and window sizes, confusion counts pooled over genes) and `pairs.tsv` (human hg38 100-way chr17; mouse mm39 35-way chr19; fly dm6 124-way chr2L; worm ce11 135-way chrIII; 8 genes each, seed 20260909). Full run: 32 windows, 19 informants, W 90/150/300/600, 365 s, 92 MB cache, no missing windows or absent informants; `results/summary.tsv`, `per_gene.tsv`, `summary.md`, `manifest.json` committed (150 KB). Floor: nucleotide F1 0.64-0.68 (MCC 0.62-0.65) for human-mouse, mouse-human and mouse-dog at 150 bp; fly-yakuba F1 0.81 (MCC 0.76) at 300 bp; chimp (0.01) calls nothing; worm pairs (all 1.1-1.3 on the 135-way tree) F1 <= 0.39 with aligned-only MCC <= 0.14, a failure not a floor. Sample composition behind the write-up (CDS fraction of window bases; share of CDS bases in segments <60 / 60-100 / 100-200 / 200-400 / 400+ bp): human 0.13; 0.04/0.08/0.36/0.21/0.31. mouse 0.19; 0.03/0.11/0.32/0.28/0.27. fly 0.44; 0.01/0.04/0.08/0.16/0.72. worm 0.43; 0.01/0.06/0.21/0.33/0.39. Write-up in `baselines/kaks/README.md`: five failure modes (too close, too far, short exons vs fixed windows, introns inside the window in compact genomes, no stop veto) and the inputs the new model needs (many informants with the tree, codon resolution with an explicit frame carried across introns, boundary signals, cheap precision signals, alignment coverage for the 14 panel species without a UCSC track). Not tested: fast-evolving genes and conserved non-coding elements; both need whole-chromosome runs. Next: respond to review; if the reviewer wants it, add `--max-ref-stops 0` as a second row of the floor and a whole-chromosome run through `benchmark/score.py`.
- 2026-09-10 human: -> done (accepted via GitHub issue #15).
