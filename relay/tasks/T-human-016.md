---
id: T-human-016
title: Controlled encoder comparison: DNA-only, MSA, tree tokens, patristic bias
status: open
owner: null
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: null
depends_on: [T-human-015]
touches: [model/, docs/design/, docs/results/]
pr: null
---

## Goal

Phase 4 milestone 4, the experiment the project exists to run: the
controlled encoder comparison of proposal section 7, items 1 to 3 and 5.

- Arms: DNA-only (A); same MSA, no tree; same MSA, tree tokens; same MSA,
  patristic bias; and, only if `nekrut/axomeme` or an equivalent
  implementation becomes readable, an audited MDS/Tree-RoPE arm. Every
  comparative arm shares the frozen A and decoder from T-human-014 and
  trains only its residual. Same informants, masks, support, split and
  effective capacity across arms.
- Check the representation symmetries of section 7.3 at emission-score
  level before comparing labels: taxon permutation, rotation and sign,
  removing taxa, rescaling branch lengths, K=1 fallback, all-masked
  attention, reverse complement, chunk seams.
- Evaluate on the held-out panel with `benchmark/score.py` and
  `benchmark/report.py`: exact transcript, exon, splice, start/stop and
  locus metrics, fusions and splits, GC/intron/exon strata; paired and
  cross-clade aggregates kept separate; Tetrahymena reported and excluded
  from ranking. Report aligned and unaligned subsets with the whole-genome
  denominator. At least three seeds per arm.
- Run `benchmark/leakage_check.py` before every held-out evaluation.

This is the largest compute item; it runs through gagarin under the Phase 4
cap, and the coordinator may split it across several grants.

## Definition of done

- `docs/results/encoder-comparison.md` with the per-arm, per-species tables,
  seed uncertainty, the symmetry checks, and one paragraph per arm on
  whether geometry helped, hurt or did nothing, with the cost of each arm
  in the same table.
- The KA/KS floor from `baselines/kaks/` and the strongest external tool
  measured in T-human-009 as rows in the same table.
- A `proposal` to `human` stating which arm, if any, should become the
  model, and what Phase 5 should be.
- Pull request, at least two `review` messages, coordinator merge, and a
  coordinator `decision`.

## Log

- 2026-09-15 human: created.

