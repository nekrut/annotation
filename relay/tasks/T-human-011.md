---
id: T-human-011
title: Design proposal for the geometric gene prediction model
status: in_progress
owner: stalin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-15T02:37:38Z
depends_on: [T-human-006, T-human-007, T-human-008, T-human-009]
touches: [docs/design/]
pr: null
---

## Goal

Turn the review, data inventory, benchmark, and cost baseline into two or
three concrete candidate designs in the spirit of HyphAeon, and recommend
one for prototyping.

## What to produce

`docs/design/proposal.md` with, for each candidate:

- Inputs: which alignment source, how many taxa, window length, whether
  the tree is given (as in Tree-RoPE) or inferred, which extra tracks
  (conservation, RNA-seq, splice-site motifs) if any.
- Geometry: what symmetries and coordinates are encoded rather than learned
  (codon periodicity and frame, strand, phylogenetic distance, splice-site
  position relative to exon boundary, intron length as a clade-conditioned
  prior), and how. Say explicitly what the network no longer needs to
  learn because of each choice.
- Output and decoding: per-base labels, a differentiable HMM or CRF layer,
  or direct structure prediction, and how long-range intron consistency is
  enforced.
- Parameter budget (target under 5M), training data plan against the
  leakage rules, expected seconds per megabase, and the three most likely
  ways the design fails.
- How it degrades when alignments are absent (a genome with no close
  relatives) and how it handles very long or very short introns.

End with a ranked recommendation and post a `proposal` message to `human`
requesting a `decision` on which candidate to prototype in Phase 4.

## Definition of done

Pull request merged after `review` messages from at least two other agents,
and a `decision` from the coordinator.

## Log

- 2026-09-09 human: created.
- 2026-09-15 stalin: claimed.
- 2026-09-15T00:41:59Z stalin: first design pass in `relay/artifacts/T-human-011/design-notes.md`, with reproducible matrix/attention cost scenarios in `budget-arithmetic.tsv`. Compared a compact DNA/CRF control, optional narrow comparative refinement, and a splice-graph decoder; fixed split/label/fallback rules and separated patristic-distance bias from unaudited MDS/Tree-RoPE. The sparse comparative scenario leaves conditional budget headroom but dense support exceeds the CPU target; these are arithmetic estimates, not measured inference. Next: complete decoder/state and full-layer resource accounting, refine candidate-density policy, then submit the work-branch proposal for two reviews and a human decision. No training or cluster job run.
