---
id: T-human-011
title: Design proposal for the geometric gene prediction model
status: review
owner: stalin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-15T08:42:56Z
depends_on: [T-human-006, T-human-007, T-human-008, T-human-009]
touches: [docs/design/]
pr: https://github.com/nekrut/annotation/pull/26
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

- 2026-09-15T01:42:01Z stalin: revised `design-notes.md` in place with a codon-prefix CRF specification that preserves split start/stop codons across introns, exact delayed entry for the minimum intron length, partial/ambiguous-base rules, and explicit traceback resource regimes. Derived 165 active delayed-entry states and up to 1.648 billion canonical-input transition candidates/Mb for both strands plus replay; these are design arithmetic, not timings. Added emission-regeneration scenarios to `budget-arithmetic.tsv`: the sparse B core alone then exceeds the conditional CPU target, so the proposed first regime spools emissions and charges scratch I/O. Rechecked the NCBI genetic-code source, all seven encoder scenario calculations, decoder storage formulas and hand examples. Next: full stem/head parameter inventory, candidate-density/support policy and label-ceiling counts, then the work-branch PR and two review requests. Task remains in progress; no model implementation, training or cluster run.

- 2026-09-15T02:47:37Z stalin: replaced A/B parameter reservations with explicit layer inventories in `design-notes.md` (455,841 and 495,021 learned scalars including the pooled decoder); specified pooling-aligned A chunks and per-orientation B tiles with a 1.53125 halo factor. Recomputed all seven `budget-arithmetic.tsv` scenarios including stems, projections and output heads: sparse B at K=8 and mean allocated support 0.05 now costs 327.1901 GFLOP/Mb, or 10.9063 conditional CPU-s/Mb at the assumed throughput; emission regeneration exceeds the CPU target even for A. These are arithmetic estimates, not measured inference. Independently enumerated layer counts, reproduced every TSV value from the embedded formula, and checked all core/frame halo dependencies and local links. Processed the cost-baseline cluster-request broadcast; no question required an answer. Next: seed/support density and label-ceiling audit, a cap/fallback schedule that accounts for buffering, then C edge accounting and the work-branch review PR. Task remains in progress; no model implementation, training or cluster run.

- 2026-09-15T03:44:55Z stalin: revised `design-notes.md` in place with an A-score support scan, a 98,304-base group buffer, a deterministic per-group CPU tile quota and separate A-CPU/B-capped-CPU/B-gated-GPU/B-full-GPU regimes. Audited four checksum-verified train-species GFFs with the accepted scorer filters: CDS-only oracle tile fractions are 0.451476 yeast, 0.242737 worm, 0.125079 fly and 0.020037 mouse, after explicit strand accounting; the first three exceed the conditional CPU ceiling if every retained CDS receives refinement. Added reproducible single-path interval-packing bounds, explicitly distinguished from benchmark transcript F1, and stem/decoder/I/O sensitivity arithmetic. Read both inbox notes; no question required an answer. Verified tile orientation, tiny exhaustive interval packing, quota/buffer arithmetic, all seven budget TSV rows and local links; rechecked axomeme public access (404). Next: actual longest-CDS/conflict-mask label audit and C edge accounting, then prepare the work-branch review PR. Task remains in progress; no model implementation, training or cluster job run.

- 2026-09-15T04:45:15Z stalin: revised `design-notes.md` in place with an exact longest-CDS selection audit using the existing cutter and four checksum-verified train GFFs, a conservative full-span same-strand conflict mask, and all/selected/unmasked boundary counts. The topology mask excludes 138/132/530/370 representatives in yeast/worm/fly/mouse; retained representatives still require the declared partial/phase/translation exception checks. Reproduced the inbox quota finding and replaced forfeited per-group slots with forward credit on the same sequence/orientation: mouse oracle CDS tiles served increase from 168,501 to 282,530 of 284,713 at the identical total slot allowance; no measured gate recall or runtime is claimed. Added reproducible flank and offline-pooling comparisons, reran the embedded audit, checked interval masks against 500 independent all-pairs cases and carried quotas against 20,736 exhaustive cases, and verified local links and arithmetic. Processed the one unread note; no question required an answer. Next: specify partial-CDS/grammar exception handling and C edge accounting, then prepare the work-branch review PR. Task remains in progress; no model implementation, training or cluster job run.

- 2026-09-15T05:45:44Z stalin: revised `design-notes.md` in place with a structured-label admission contract for partial CDS, raw phase/row checks, conservative exception masks and independent auxiliary-site supervision. Extended the checksum-verified four-train-GFF reproduction to raw CDS metadata before scorer merging: after topology masking, 6/28/455/102 yeast/worm/fly/mouse representatives carry exclusion/review flags; 5,858/19,811/12,980/21,711 have none of these flags, but still need FASTA/grammar validation. Specified C's 577-scalar edge head (495,598 total), endpoint-first range pruning with E<=32V, direct exon-scan O(P*W) accounting, top-path/output costs and pilot-frozen fallback limits. Folded in both independent inbox arithmetic checks; no question needed an answer. Reproduced all existing annotation/quota tables and the new metadata counts, checked 18 strand-aware raw-CDS fixtures, graph arithmetic and local links. Next: final consistency pass and prepare the work-branch proposal PR for two reviews and the human decision. Task remains in progress; no model implementation, training or cluster run.
- 2026-09-15 stalin: -> review (https://github.com/nekrut/annotation/pull/26).

- 2026-09-15T06:42:56Z stalin: completed the consistency pass, revised the owned artifact in place, and submitted `docs/design/proposal.md`, `annotation-audits.md` and `budget-arithmetic.tsv` on `work/T-human-011-stalin` at `d8655e0` in PR #26 (https://github.com/nekrut/annotation/pull/26). The proposal ranks B with mandatory A first and C later, fixes the initial frozen-A/decoder fitting order and C's synthetic partial-edge rules, incorporates both independent inbox checks, and ends with bounded Phase 4 milestones and the coordinator decision request. Reproduced all seven budget TSV rows and A/B/C scalar totals, checked the summary time arithmetic, Python-listing syntax, local links/anchors, reference labels, unchanged appendix code and whitespace. No question needed an answer. Requested formal reviews from marx and lenin through relay and a decision from human; task is in review. Next: address the formal reviews and coordinator feedback. No model implementation, training, FASTA admission audit or cluster job run.
