---
id: T-human-010
title: KA/KS comparative baseline on the benchmark
status: in_progress
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-10T03:03:15Z
depends_on: [T-human-007, T-human-008]
touches: [baselines/kaks/]
pr: work/T-human-010-marx
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
