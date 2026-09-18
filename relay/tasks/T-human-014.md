---
id: T-human-014
title: Implement and measure candidate A end to end on pilot chromosomes
status: in_progress
owner: engels
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-18T19:20:50Z
depends_on: [T-human-013]
touches: [model/, docs/design/, benchmark/]
pr: null
---

## Goal

Phase 4 milestone 2. Implement candidate A (section 3 of the proposal:
compact DNA encoder, at most 5 M parameters, the decoder from T-human-013)
and measure it end to end on train-only pilot chromosomes.

- Train A on development chromosomes of the train species only, under the
  Phase 4 training cap in `relay/TASK.md`. Select checkpoints on train
  development chromosomes only.
- Measure preprocessing, encoder, decoder, traceback, scratch I/O and
  output separately, with peak host and device memory and the hardware,
  using `/usr/bin/time -v` and the same conventions as
  `docs/cost-baseline/measured.tsv`. Report the CPU regime on one core and
  the GPU regime with the multi-worker decoder accounting of section 6.1,
  worker count and aggregate memory included.
- Compare against the accepted budgets (15 CPU-s/Mb, 0.5 GPU-s/Mb, 8 GB)
  and against AUGUSTUS on the same machine using the S. pombe
  normalization of section 6.1. Score with `benchmark/score.py` on the
  development chromosomes.
- If A misses a target, report the failed regime and propose the revision
  before any positive CPU allowance is assigned to B.

GPU work goes through gagarin: post an `alert` with
`relay/templates/compute-request.md` filled in, wait for the `decision`.

## Definition of done

- `model/a/` with training and inference entry points, a config that
  reproduces the reported run, and a run manifest (data, seed, commit,
  hardware).
- New rows in `docs/cost-baseline/measured.tsv` for A on at least
  S. pombe and one metazoan development chromosome, both regimes, and a
  `docs/design/a-pilot.md` that puts them against the budget with the
  section 6.1 sensitivity table filled with measured numbers.
- Actual GPU-hours and CPU-hours recorded in this task's log.
- Pull request, at least two `review` messages from other agents (this is
  the measurement everything else depends on), coordinator merge, and a
  `decision` on whether A meets the budget.

## Log

- 2026-09-15 human: created.
- 2026-09-18 engels: claimed.
- 2026-09-18 engels: -> in_progress.
- 2026-09-18T12:23Z engels: claimed after T-human-013 merged; status in_progress.
  Audited the merged grammar and admitted-label interfaces. All ten manifest
  checksums, panel GFF checksums, summary counts and composite transcript
  identities agree (179,227 admitted / 17,451 masked); this rechecks committed
  results, not the original FASTAs. Existing grammar suite: 36 tests passed
  in 4.062 s on the available Python 3.14.4, with Python 3.11 verification
  still required. Detailed evidence, reproduction and implementation plan:
  [engels-0042](../messages/20260918T122324Z-engels-0042.md).
  Findings: the scalar grammar still needs a differentiable chain loss;
  manifests omit exon geometry and auxiliary-site coordinates, so the loader
  must rejoin checksummed raw GFF3/FASTA; trainable partial-end families and
  motif tables must match the proposal inventory. S. pombe stays held out,
  used only after freezing for the required runtime-normalization evaluation;
  fitting, checkpoint selection and development scoring use train species.
  No training or model benchmark run; training/cluster CPU-hours 0, GPU-hours
  0 (local verification CPU-hours unmetered). Next tick: implementation on
  work/T-human-014-engels, then a concrete compute request through gagarin.

- 2026-09-18T13:25Z engels: renewed lease and processed all three unread notes; no
  question required an answer. Verified an additional structured-target
  requirement: emission/phase/site masks alone count an undeclared E0 path
  for a complete gene touching the true 5' edge. Both decoders produce the
  analytically predicted excess log(7/6); explicit gold initial/terminal
  support removes it without changing the shared edge weights. An in-memory
  adapter check passed 464 synthetic cases (928 decoder configurations),
  including both strands, genetic codes 1/6, partial-end combinations,
  missing-prefix phases, acceptor conversion, and a seam at every base.
  Reproduction and limitations: [engels-0043](../messages/20260918T132457Z-engels-0043.md).
  Latest sweep: 0.1754 CPU seconds; previous equivalent sweep: 0.1793 CPU
  seconds (Python 3.14.4), total 0.3547 seconds (0.00009853 local CPU-hours).
  Other inspection/counterexample CPU time unmetered. Training/cluster
  CPU-hours 0; GPU-hours 0. No source-FASTA audit or pilot measurement.
  Implementation-file creation awaits the operator's clarification of the
  supplied AGENTS.md restrictions versus the work-branch rule; all writes
  this tick remain within owned relay metadata and new messages. Task stays
  in_progress. Next: implement the boundary-aware loader and reserved train
  chromosome split on work/T-human-014-engels once path scope is resolved;
  retain masks/crop state and prepare the executable before requesting compute.

- 2026-09-18T14:24Z engels: renewed lease and processed all three unread notes; no
  question required an answer. Checked the masked-label chain-loss contract:
  unknown targets remove support constraints while retaining shared model
  scores and latent gradients. Reference/delayed partitions agreed on 48
  synthetic cases with seams at every base; 288 finite-difference checks of
  emission gradients passed (maximum absolute error 4.631317551684333e-11).
  A closed-form example shows that zeroing numerator potentials can produce
  a negative loss for an entirely unknown target. Reproducer and limitations:
  [engels-0044](../messages/20260918T142348Z-engels-0044.md).
  Measured local CPU 0.296213548 seconds (0.0000822815 CPU-hours), wall
  0.29623782890848815 seconds, Python 3.14.4; other inspection unmetered.
  Training/cluster CPU-hours 0; GPU-hours 0. No fitting or source-FASTA audit.
  Scope clarification for implementation-file creation was requested from the
  operator and remains pending; no implementation files, PR or compute request.
  Task remains in_progress. Next: implement the boundary-aware loader and
  reserved train split on work/T-human-014-engels when file scope is resolved,
  then check the differentiable loss against this emission-gradient oracle.

- 2026-09-18T15:23Z engels: renewed lease; processed three unread notes, no
  question required an answer. Verified training-crop intron conditioning:
  retain mandatory-prefix age, S/E codon state, and survival-weighted latent
  duration components. Reference seeds and delayed checkpoint continuation
  match closed-form weights on 924 crops (798 mandatory-prefix / 126 tail,
  42 beginning at an acceptor), maximum absolute error 4.440892098500626e-16.
  Three counterexamples show that resetting age admits a subminimum intron,
  resetting component weights changes the target, and using true-edge J at
  an immediate acceptor removes a legal path. Reproducer and scope:
  [20260918T152338Z-engels-0045](../messages/20260918T152338Z-engels-0045.md).
  Local measured CPU 0.067340822 s (0.0000187058 CPU-hours), wall
  0.067355989 s, Python 3.14.4; other inspection unmetered. Training/cluster
  CPU-hours 0; GPU-hours 0. No learned loss, real-data crop, or benchmark run.
  Implementation file-scope clarification remains pending; only owned relay
  metadata and a new note changed. Task stays in_progress, no PR or compute
  request. Next: implement the boundary-aware adapter with the crop-state
  contract on work/T-human-014-engels when implementation writes are allowed.

- 2026-09-18T16:23Z engels: renewed lease and processed both unread notes; no question
  required an answer. Checked crop-local delayed pending-donor initialization
  and conditional intron-censored targets: 924 completed + 252 censored crops,
  2,352 decoder configurations and 2,352 masked-target rejections pass.
  Maximum score error 8.881784197001252e-16; 3,024 duration-gradient finite
  differences pass, maximum error 2.2283952461066292e-10. Identified a crop
  integration hazard: reusing whole-sequence pending-exit sums with negative
  rebased donor coordinates can double-count crop emissions (0.42 vs 0.15).
  Reproducer and limitations: [engels-0046](../messages/20260918T162315Z-engels-0046.md).
  These are in-memory scalar checks, not a public crop API or autograd test.
  Python 3.14.4; measured local CPU 0.485892064 s (0.000134970018 CPU-hours),
  wall 0.4860360769 s; other inspection unmetered. Training/cluster CPU-hours
  0; GPU-hours 0. No fitting, source-data evaluation or benchmark run.
  Implementation-file scope clarification remains pending; only owned relay
  metadata and a new note changed. Task stays in_progress, no PR or compute
  request. Next: implement the boundary-aware adapter and conditional loss
  on work/T-human-014-engels when implementation writes are permitted.

- 2026-09-18T17:23Z engels: renewed lease; processed both unread notes, no
  question required an answer. Audited all ten train admission manifests for
  crop-length exposure, with manifest MD5s, panel GFF checksums, counts and
  composite identities verified. Of 179,227 admitted representatives, 72,350
  (40.37% pooled; 36.51% equal-species mean) span more than the 3,072-base
  encoder core; 57,369 of those have CDS length at or below 3,072. Even
  49,152-base blocks cannot contain 9,658 complete spans. This is a length-only
  exclusion floor for a hypothetical whole-gene sampler, not an implemented
  sampler or model measurement. Evidence and reproducer:
  [engels-0047](../messages/20260918T172256Z-engels-0047.md).
  Python 3.14.4; measured local CPU 0.375624981 s (0.000104340273 CPU-hours),
  wall 0.375660879 s; other inspection unmetered. Training/cluster CPU-hours
  0; GPU-hours 0. No raw FASTA/GFF audit, fitting or held-out evaluation.
  Implementation-file scope clarification requested from the operator remains
  pending; only owned relay metadata and a new note changed. Task remains
  in_progress, no PR or compute request. Next: implement the boundary-aware
  loader and conditional loss on work/T-human-014-engels when implementation
  writes are allowed; do not turn encoder core length into a gene-length cap.
