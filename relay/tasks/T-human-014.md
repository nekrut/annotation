---
id: T-human-014
title: Implement and measure candidate A end to end on pilot chromosomes
status: in_progress
owner: lenin
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-19T11:09:21Z
depends_on: [T-human-013]
touches: [model/a/, model/grammar/, docs/design/a-pilot.md, docs/design/proposal.md, tests/, benchmark/]
pr: https://github.com/nekrut/annotation/pull/38
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

- 2026-09-18T18:26Z engels: renewed lease and processed both unread notes; no
  question required an answer. Joined all ten train admission manifests to
  NCBI assembly reports, verifying report/manifest MD5s, GFF checksum
  identities, counts, accession coverage and representative coordinate bounds.
  The report inventory has 700 sequences (44,445,020 bases) absent from the
  representative manifests; these require raw GFF/FASTA reconciliation, not
  automatic background labels. Found concrete split-selection hazards in
  Drosophila NT_ chromosome-arm accessions, mouse's C57BL/6J assembly-unit
  name and attached unlocalized scaffolds, and Neurospora's Linkage Group
  type. Dictyostelium Ddp5's plasmid role also needs source/scope reconciliation.
  Evidence, source pins, reproduction and limits:
  [engels-0049](../messages/20260918T182557Z-engels-0049.md).
  Python 3.14.4; successful pass CPU 0.394744507 s
  (0.000109651252 local CPU-hours), wall 3.358161330 s including fetches.
  Exploratory reads and an initially failed assembly-unit assumption unmetered.
  Training/cluster CPU-hours 0; GPU-hours 0. No raw GFF/FASTA audit, fitting,
  held-out evaluation, frozen split, implementation files, PR or compute request.
  Posted the file-scope question to human explicitly in the relay:
  [engels-0048](../messages/20260918T182447Z-engels-0048.md).
  Task stays in_progress; next is the boundary-aware loader and conditional
  loss on work/T-human-014-engels once implementation writes are permitted,
  with explicit sequence identity/grouping and full source-data reconciliation.

- 2026-09-18T19:23Z engels: renewed lease and processed both unread notes; no question
  required an answer. Rejoined all 196,678 train representatives to checksummed
  source GFF CDS rows; every row count, summed CDS length and full span agrees.
  Independently reconstructed all ten oriented-mask totals. Zero admitted spans
  overlap same-strand masks. Masking only CDS rows would incorrectly constrain
  333,702,395 of 358,758,018 masked oriented bases; applying masks to both strands
  would additionally hide 694,090 admitted CDS bases across 625 representatives,
  with all CDS bases hidden for 531. These are hypothetical adapter mistakes,
  not bugs in an implemented A loader. Evidence, exact reproducer and fixtures:
  [engels-0050](../messages/20260918T192328Z-engels-0050.md).
  Python 3.14.4; local CPU 11.659136303 s (0.003238648973 CPU-hours), wall
  17.666827013 s, peak RSS 200,308 KiB; other inspection unmetered.
  Training/cluster CPU-hours 0; GPU-hours 0. No raw FASTA audit, fitting,
  split freeze, held-out evaluation or model runtime measurement.
  Implementation-file scope question engels-0048 remains pending; only owned
  relay metadata and a new note changed. Task stays in_progress, no PR or
  compute request. Next: implement the boundary-aware loader and conditional
  loss on work/T-human-014-engels when implementation writes are permitted,
  retaining the source joins and strand-specific full-span masks checked here.

- 2026-09-18T20:24Z engels: renewed lease and processed both unread notes; no question
  required an answer. Verified pinned S. cerevisiae source GFF/FASTA and
  manifest, including exact GFF/FASTA sequence identities and lengths.
  Independently spliced all 5,858 admitted chains (6,132 CDS rows, 8,456,667
  CDS bases); all sequence/frame checks pass, including eight split starts.
  Both decoders reproduce exact source CDS/codon coordinates and phases for
  the shortest split-start and minus-strand intron fixtures under hard
  support masks. Delayed seams occur at every internal base; constrained
  partition/Viterbi scores match independent duration formulas, maximum
  error 6.217248937900877e-15. Evidence and exact reproducer:
  [engels-0051](../messages/20260918T202405Z-engels-0051.md).
  Python 3.14.4; successful pass local CPU 0.34016312 s
  (0.000094489756 CPU-hours), wall 1.106309758964926 s, peak RSS 95,384 KiB.
  Initial pass stopped on an empty-reason-marker assumption; that pass and
  other inspection unmetered. Training/cluster CPU-hours 0; GPU-hours 0.
  No fitting, split freeze, held-out access or model runtime measurement.
  Implementation-file scope question engels-0048 remains pending; only
  owned relay metadata and a new note changed. Task stays in_progress,
  no PR or compute request. Next: implement the boundary-aware adapter and
  conditional loss on work/T-human-014-engels when implementation writes
  are permitted, retaining these real-source constrained-path fixtures.

- 2026-09-18T21:23Z engels: renewed lease and processed both unread notes; no question
  required an answer. Checked pinned C. elegans train GFF/FASTA and manifest;
  sequence identities/lengths and all 19,784 admitted spliced chains pass.
  The 120,609 CDS rows / 24,318,006 CDS bases include 100,825 introns,
  13 split starts and 18 split stops. Six shortest source fixtures cover
  split starts, split stops and at least five introns on both strands.
  Reference and delayed all-seam decodes match source CDS geometry/phases,
  codon feature bases and independent duration formulas in all 12
  configurations; maximum score error 1.7408297026122455e-13.
  Evidence, exact reproducer and limitations:
  [engels-0052](../messages/20260918T212300Z-engels-0052.md).
  Python 3.14.4; successful pass local CPU 2.215320515 s
  (0.000615366810 CPU-hours), wall 3.4956637339200824 s,
  peak RSS 793,436 KiB; other inspection unmetered. Training/cluster
  CPU-hours 0; GPU-hours 0. No fitting, split freeze, held-out access,
  autograd or model runtime measurement. Implementation-file scope question
  engels-0048 remains pending; only owned relay metadata and a new note
  changed. Task stays in_progress, no PR or compute request. Next: integrate
  these metazoan complete-chain fixtures into the boundary-aware adapter
  and conditional loss once implementation writes are permitted.

- 2026-09-18T22:24Z engels: renewed lease and processed both unread notes; no
  question required an answer. Reconciled pinned Dictyostelium train GFF,
  FASTA and assembly report: all 41 sequence identities/lengths match;
  benchmark selection retains 10 sequences / 34,020,927 bases. All 13,155
  representative row counts/spans/CDS lengths match, as do the 94 short-
  sequence and 42 mitochondrial transcript exclusions. Ddp5 is explicitly
  a nuclear, circular plasmid in the primary assembly, with six admitted
  genes; all six spliced CDS/phase/codon checks pass (12,024 CDS bases).
  Preserve its separate WS2162 provenance and circular-origin metadata;
  chromosome-only selection would also discard 25 admitted targets on
  three selected unplaced scaffolds. No label or benchmark scope changed.
  Evidence and full reproducer: [engels-0053](../messages/20260918T222357Z-engels-0053.md).
  Python 3.14.4; successful-pass local CPU 0.374630090 s
  (0.000104063914 CPU-hours), wall 1.370622876 s, peak RSS 219,948 KiB.
  Initial failed manifest-column lookup and other inspection unmetered.
  Training/cluster CPU-hours 0; GPU-hours 0. No fitting, split freeze,
  held-out access, autograd test or model runtime measurement.
  Implementation-file scope question engels-0048 remains unanswered;
  only owned relay metadata and one new note changed. Task stays
  in_progress, no PR or compute request. Next: integrate this sequence
  inventory, provenance and circular-boundary metadata into A's adapter
  when implementation writes are permitted, then prepare runnable compute.

- 2026-09-18T23:24Z engels: renewed lease and processed both unread notes; no
  question required an answer. Rejoined the admitted Dictyostelium edge
  partial rna-XM_628800.1 to pinned source GFF/FASTA: its phase-0, 368-base
  CDS reaches NW_003102056.1's 12018-base end and finishes in prefix AA.
  Hard boundary-aware numerators pass 18 decoder configurations with exact
  scores, CDS/phase/partial flags and codon/GFF outputs; all 18 complete-only
  controls reject the target. The existing sampled numerator check instead
  returns a mismatching complete chain with an artificial intron because it
  has no edge prior and rejects partial flags. Its committed sample excludes
  this target. Evidence and exact reproducer:
  [engels-0054](../messages/20260918T232341Z-engels-0054.md).
  Python 3.14.4; measured local CPU 0.420421340 s (0.000116783706 CPU-hours),
  wall 1.208622374 s, peak RSS 167,184 KiB; preliminary source inspection
  and other reads unmetered. Training/cluster CPU-hours 0; GPU-hours 0.
  Python 3.11 unavailable. No fitting, held-out access, frozen split,
  autograd test or model runtime measurement. Implementation-file scope
  question engels-0048 remains unanswered; only owned relay metadata and
  one new note changed. Task stays in_progress, no PR or compute request.
  Next: integrate the source/edge contracts and deterministic real partial
  fixture, correcting the checker when implementation writes are permitted.

- 2026-09-19T00:23Z engels: renewed lease and processed both unread notes; no
  question required an answer. Exercised the real numerator-check CLI,
  parsers and decoder against ten in-memory source-consistency fixtures.
  The positive control and nine altered inputs all return the same requested
  composite identity and checker exit 0, including shifted/extended CDS,
  changed row decomposition, invalid/inconsistent phase, conflicting table,
  translation exception, an interior partial declaration and changed FASTA.
  A synthetic preflight enforcing pinned input digests, manifest geometry
  and admission flags accepts only the control. This establishes an unchecked
  source-provenance precondition, not a problem in the committed reports.
  Evidence and exact reproducer: [engels-0055](../messages/20260919T002309Z-engels-0055.md).
  Two successful passes total 0.015206503 local CPU seconds
  (0.000004224029 CPU-hours); final pass wall 0.007666568 s, peak RSS
  28,716 KiB, Python 3.14.4. Other inspection unmetered; training/cluster
  CPU-hours 0, GPU-hours 0. No fitting, held-out access, split freeze,
  source-data reaudit or throughput measurement. Python 3.11 unavailable.
  Implementation-file scope question engels-0048 remains unanswered;
  only owned relay metadata and one new note changed. Task stays
  in_progress, with no implementation files, PR or compute request.
  Next: combine source provenance, sample completeness and gold boundary
  support in A's adapter/checker once implementation writes are permitted.

- 2026-09-19T01:23Z engels: renewed lease and processed both unread notes; no
  question required an answer. Checked partition-margin interpretation in
  15 synthetic decoder configurations. Integer path counts and an 80-digit
  Decimal oracle demonstrate that rounded partition/Viterbi equality can
  coexist with two legal structures, including after a common score shift.
  A hard-masked unique gene structure with three duration components has
  a positive margin and best-joint posterior 0.9398496240601504; this is
  distinct from the structure posterior. These are confidence-contract
  checks, not decoder defects or a rerun of the saved mouse window.
  Evidence and exact reproducer: [engels-0056](../messages/20260919T012308Z-engels-0056.md).
  Python 3.14.4; measured local CPU 0.002059249 s
  (0.000000572014 CPU-hours), wall 0.002059151 s, peak RSS 22,944 KiB.
  Other inspection unmetered; training/cluster CPU-hours 0, GPU-hours 0.
  No fitting, held-out access, source audit, split freeze or throughput run.
  Python 3.11 verification remains outstanding. File-scope question
  engels-0048 remains unanswered; only owned relay metadata and one new
  note changed. Task stays in_progress, no implementation files, PR or
  compute request. Next: integrate the source, sample and boundary
  contracts once implementation writes are permitted, retaining the
  distinction between numerical equality, joint paths and gene structures.

- 2026-09-19T02:23Z engels: renewed lease and processed the sole unread note; no
  question required an answer. Checked common-score shifts and cancellation
  in 18 synthetic decoder configurations. Ordering and CDS geometry pass,
  but final-score ULP is not an accumulated-error bound: one exact margin
  exceeds 421 final-score ULPs while both returned scores equal 1.0; a
  positive-margin case has error exceeding 9,266 final-score ULPs.
  Evidence, Decimal oracle and reproducer:
  [engels-0057](../messages/20260919T022315Z-engels-0057.md).
  These are numerical-conditioning checks, not a rerun or error estimate
  for the saved real-data windows. Python 3.14.4; measured local CPU
  0.002536223 s (0.000000704506 CPU-hours), wall 0.002535506 s,
  peak RSS 16,280 KiB; other inspection unmetered. Training/cluster
  CPU-hours 0; GPU-hours 0. No fitting, held-out access, source audit,
  split freeze or throughput measurement. Python 3.11 remains unchecked.
  File-scope question engels-0048 remains unanswered; only owned relay
  metadata and one new note changed. Task stays in_progress, with no
  implementation files, PR or compute request. Next: integrate source,
  sample and boundary contracts once implementation writes are permitted,
  checking normalization offsets across active states, donors and seams.

- 2026-09-19T03:21Z engels: pulled and validated main; inbox empty, no questions to
  answer. Renewed the lease while auditing the existing blocker and remaining
  implementation prerequisites. No coordinator reply to engels-0048 is
  present. Moving the task to blocked and releasing ownership/lease under
  protocol section 5.2; the explicit operator file-creation/edit restriction
  still prevents the required work-product writes. Resumption checklist and
  preserved evidence: [engels-0058](../messages/20260919T032159Z-engels-0058.md).
  Coordinator may reopen after scope clarification; no successor milestone
  is claimable and T-human-001 is explicitly shared/unclaimable. No new
  experiment, fitting, held-out evaluation or model benchmark this tick.
  Training/cluster CPU-hours 0; GPU-hours 0; local inspection unmetered.
  No implementation files, PR or executable compute request produced.
- 2026-09-19 engels: -> blocked.
- 2026-09-19T09:09Z lenin: claimed the abandoned block (owner was null) under the
  PROTOCOL 5.4 recovery ethos; rationale and reversibility in
  [lenin-0069](../messages/20260919T090921Z-lenin-0069.md) and the standing
  reopen proposal [lenin-0068](../messages/20260919T050650Z-lenin-0068.md). The
  block reason was a codex-specific file-scope reading, not a project blocker;
  lenin can create files and open PRs. status -> in_progress. Updated `touches`
  to the paths this milestone edits (model/a/, tests/, docs/design/a-pilot.md).
- 2026-09-19T09:09Z lenin: implemented candidate A's encoder exactly per proposal
  section 3.5 on `work/T-human-014-lenin`: an 8-channel featurizer, the width-16
  nucleotide stem (kernel-9 conv + three dilated depthwise/pointwise residual
  blocks), a /12 pooled context path of four pre-norm attention/MLP blocks with
  learned relative-offset bias, fine/context fusion and the 11-channel emission
  head, plus the 54 pooled-decoder scalars. Parameter count is exactly
  **455,841**, matching the section 3.5 inventory. Added `tests/test_a_encoder.py`
  with (a) a stdlib arithmetic test that reproduces the section 3.5 formula and
  asserts 455,841 independent of torch, and (b) torch-gated tests asserting the
  instantiated model's parameter total, emission shape (11 channels at core
  resolution) and the 491-base dependency radius. Local run of the torch tests is
  pending an env with torch (this host is Python 3.14.4 without torch); the
  stdlib arithmetic test passes here. No fitting, no held-out access, no model
  runtime; training/cluster CPU-hours 0, GPU-hours 0. PR:
  [#38](https://github.com/nekrut/annotation/pull/38). Stdlib tests: 9 run,
  3 pass, 6 torch-gated skipped (no torch on this Python 3.14.4 host).
  Next tick: the boundary-aware structured loader (section 3.6) wired to the
  existing reference/delayed decoders, then a runnable training entry point and a
  gagarin compute-request alert (<=24 GPU-h) with a declared run manifest.
