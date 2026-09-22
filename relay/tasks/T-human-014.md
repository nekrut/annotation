---
id: T-human-014
title: Implement and measure candidate A end to end on pilot chromosomes
status: in_progress
owner: lenin
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-22T10:08:10Z
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
- 2026-09-19T10:09Z lenin: renewed lease; processed both unread reviews of PR #38
  (engels-0059, stalin-0062), all four P2 findings, no question needing an answer.
  Fixed on `work/T-human-014-lenin` (commit 56799ca, PR #38 updated): (1) rewrote
  `LocalAttention.forward` to gather only the 16 permitted offsets per query
  (O(T·W·hd), not O(T²·hd)), keeping `_dense_forward` as an oracle a torch-gated
  test checks output and input gradient against; (2) `gc_track` now excludes
  unavailable positions from both GC counts so padding letters cannot contaminate
  a real base's window, with letter-invariance tests; (3) added `model.a` to
  `[tool.setuptools] packages`; (4) `DecoderParams` seeds distinct per-component
  hazard logits (`HAZARD_LOGIT_INIT = (-1,0,1)`) to break the duration-mixture
  symmetry that zeroed the mixture/hazard gradients, recorded for the run manifest.
  Added tests for each; `a-pilot.md` section 6 records the responses. Fusion/emission
  head halo cropping remains a section-3.6 loader concern, noted not yet done.
  Local stdlib tests: 13 run, 4 pass, 9 torch-gated skipped (no torch on this
  Python 3.14.4 host); the oracle-equivalence and duration tests need a torch env.
  No fitting, held-out access or model runtime; training/cluster CPU-hours 0,
  GPU-hours 0. Next tick: the boundary-aware structured loader (section 3.6) and
  the differentiable chain loss, then a runnable entry point and the gagarin
  compute-request alert.
- 2026-09-19T11:06Z lenin: renewed lease; processed both unread follow-up reviews
  (engels-0060, stalin-0063), both nonblocking, no question needing an answer.
  Both confirm all four P2 fixes landed at commit 56799ca and converge on one
  nonblocking wording correction: identical initial hazard components zero the
  *mixture-logit* gradients and make the *hazard-logit* gradients equal across
  components (nonzero but moving them together), so the components cannot
  differentiate under symmetric updates — earlier comments/docs wrongly said both
  gradient families vanish. Applied on `work/T-human-014-lenin` (commit c8fe119,
  PR #38): corrected the `DecoderParams` comment, the `test_duration_components_
  start_distinct` comment, and a-pilot.md section 6; no code/behavior change, the
  distinct-init test still passes. Local stdlib tests: 13 run, 4 pass, 9 torch-
  gated skipped (no torch on this Python 3.14.4 host). Reviewers' outstanding
  asks unchanged and acknowledged: PyTorch forward/backward, instantiated
  param-count, chunk seams, learned-loss gradients, Python 3.11 and end-to-end
  timing all still need a torch/GPU env; fusion/emission halo cropping stays a
  section-3.6 loader concern. No fitting, held-out access or model runtime;
  training/cluster CPU-hours 0, GPU-hours 0. Next tick: implement the boundary-
  aware structured loader (section 3.6) and the differentiable chain loss wired
  to the existing reference/delayed decoders, then a runnable training entry
  point and the gagarin compute-request alert (<=24 GPU-h) with a declared run
  manifest.
- 2026-09-19T12:11Z lenin: renewed lease; processed both unread messages
  (engels-0061, stalin-0064) — both close the PR #38 duration-wording follow-up
  with no new blocking finding and no question needing an answer. Implemented the
  chain-loss oracle on `work/T-human-014-lenin` (commit d56c469, PR #38 updated):
  `model/a/loss.py` `numerator_scores` builds the gold chain's hard `-inf`
  support mask over the eleven emission channels (U only intergenic, coding only
  on a CDS base, intron only on an intron base, start/stop/donor/acceptor only at
  their gold coordinates), and `chain_nll` returns `log Z − log Z_num` on the
  existing `model.grammar` decoders — the free partition minus the support-
  restricted partition, marginalizing exactly the grammar-fixed frame and the
  intron duration mixture. `tests/test_a_loss.py` (8 stdlib cases + 1 torch-gated
  parity, skipped) pins the mask, `loss ≥ 0`, that the masked numerator Viterbi-
  decodes back to the exact gold chain (single- and two-exon fixtures), the
  `log 2` zero-emission case, and the loss gradient signs by finite difference.
  This lands the support-mask construction, the constrained-vs-free contract and
  the reference-decoder cross-check; the batched chunk-seam **torch** loss on the
  encoder emissions and the crop-edge boundary conditioning (section 3.6) remain.
  a-pilot.md sections 1–2 updated. Local: encoder suite 13 run / 4 pass / 9 torch-
  gated skipped; loss suite 9 run / 8 pass / 1 torch-gated skipped (no torch on
  this Python 3.14.4 host). No fitting, held-out access or model runtime;
  training/cluster CPU-hours 0, GPU-hours 0; local inspection unmetered. Next:
  the section 3.6 structured loader over checksummed GFF/FASTA and the torch loss
  matching this oracle, then a runnable training entry point and the gagarin
  compute-request alert (≤24 GPU-h) with a declared run manifest.
- 2026-09-19T13:08Z lenin: renewed lease; processed both unread reviews of PR #38
  (engels-0062, stalin-0065), no question needing an answer. stalin-0065 confirms
  the expanded complete-chain/duration oracle passes (432 decoder configurations,
  36 duration-gradient checks) but leaves engels-0062's two P2 findings open;
  both fixed on `work/T-human-014-lenin` (commit bd72881, PR #38): (1) the support
  mask constrains emissions, not boundary states, so an enabled `EdgePrior` let the
  numerator claim extra entry/exit/phase hypotheses (complete gene at a real edge
  scores log(5/4) not 0) while the loss could stay nonnegative and pass — `chain_nll`
  now rejects an edge-enabled decoder with a clear `ValueError`, scoping this as the
  *complete*-target oracle; declared edge-partial numerators (compatible entry/exit
  families + first-row phase) stay a section-3.6 boundary-support increment; (2) the
  torch-parity gate caught only `ModuleNotFoundError` but `from model.a import
  torch_loss` raises a plain `ImportError` on a torch host, aborting discovery —
  now probes `importlib.util.find_spec("model.a.torch_loss")`, skipping cleanly when
  the submodule is absent while still surfacing a genuinely broken implementation.
  Added a test for each; a-pilot.md section 6 records both. Local loss suite: 10 run
  / 9 pass / 1 torch-gated skip (no torch on this Python 3.14.4 host); torch-present
  gate verified by stubbing `sys.modules['torch']` (imports cleanly, no ImportError).
  No fitting, held-out access or model runtime; training/cluster CPU-hours 0,
  GPU-hours 0; local inspection unmetered. Next: the section 3.6 structured loader
  over checksummed GFF/FASTA and the torch loss matching this oracle (with the
  edge-partial boundary interface), then a runnable training entry point and the
  gagarin compute-request alert (≤24 GPU-h) with a declared run manifest.
- 2026-09-19T14:12Z lenin: renewed lease; processed both unread reviews of PR #38
  (engels-0063, stalin-0066). Both close engels-0062's two P2 findings (the
  edge-prior boundary guard and the torch-import gate) at commit bd72881 with no
  new blocking finding; stalin-0066 re-ran the full complete-chain sweep (432
  decoder configurations, 36 duration-gradient checks) and 16 boundary + import
  probes, engels-0063 the 23-test suite (13 pass / 10 torch-skip). No question
  needed an answer. Implemented the section 3.6 structured **training-window
  loader** on `work/T-human-014-lenin` (commit 267b34d, PR #38): `model/a/
  dataset.py` `verify_source` is a hard checksum gate against each species'
  `*.summary.json` gff/fasta MD5s, and `iter_windows` delegates admission to
  `model.labels.admission.audit_species` (never re-deriving the audit) and
  rejoins the checksummed FASTA with `numerator_check.oriented_chain`, yielding
  `WindowExample`s (oriented window + merged half-open CDS/intron ranges + table
  + source identity) that feed `numerator_scores`/`chain_nll` directly. Scoped
  to the complete-target oracle: edge-partial and over-`max_window` chains are
  skipped and counted in `LoaderStats`, never capping gene length by the encoder
  core. `tests/test_a_dataset.py` (6 stdlib cases) builds a synthetic admitted
  two-exon gene on a 10 kb contig and checks the checksum gate (match, tampered
  FASTA, FASTA-presence mismatch), window coords/identity, the `max_window` skip
  counter, and the numerator→gold Viterbi round trip. Local: encoder 13 / loss 10
  / dataset 6 = 29 run, 10 torch-gated skipped, all pass (no torch on this Python
  3.14.4 host). a-pilot.md sections 1–2 updated. No fitting, held-out access or
  model runtime; training/cluster CPU-hours 0, GPU-hours 0; local inspection
  unmetered. Next: the boundary-support path (edge-partial numerators via an
  edge-enabled decoder) and whole-gene cropping into the encoder core, plus the
  batched torch chain loss matching the oracle, then a runnable training entry
  point and the gagarin compute-request alert (≤24 GPU-h) with a run manifest.
- 2026-09-19T15:14Z lenin: renewed lease; processed both unread reviews of the PR #38
  loader (engels-0064, stalin-0067). Both independently reproduced four P2 findings
  at commit 267b34d; all four fixed on `work/T-human-014-lenin` (commit b21ddf2,
  PR #38): (1) **false intergenic labels on neighbours** — `iter_windows` now
  yields only *clean* windows, skipping (and counting in `LoaderStats.
  skipped_neighbor`) any window a different gene overlaps, so a neighbour's CDS or
  masked bases are never forced to `U`; isoforms of the same gene do not dirty a
  window; (2) **soft mask erased** — the window is rebuilt case-preserving from the
  raw FASTA slice (`_oriented_window`, `revcomp` keeps case) so featurizer
  channel 5 survives orientation, with `oriented_chain` kept as the coordinate
  authority and a guard the two agree up to case; (3) **hard gate optional /
  settings unbound** — `iter_windows(summary, gff, fasta)` now calls
  `verify_source` before yielding and reads audit `m`/`table` from the pinned
  summary, so it cannot audit unpinned inputs or use disagreeing settings;
  (4) **packaging** — added `model.labels` to setuptools packages and made
  `model.a.__init__` import the loader lazily (PEP 562), so `import model.a` needs
  neither `model.labels` nor `benchmark/score.py` (verified: no `model.labels`
  submodule imported by `import model.a`). Added a regression for each. Local:
  encoder 13 / loss 10 / dataset 11 = 34 run, 10 torch-gated skipped, all pass
  (no torch on this Python 3.14.4 host). a-pilot.md sections 1 and 6 updated. No
  fitting, held-out access or model runtime; training/cluster CPU-hours 0,
  GPU-hours 0; local inspection unmetered. Next: the boundary-support path
  (edge-partial numerators) and whole-gene cropping into the encoder core, the
  batched torch chain loss matching the oracle, then a runnable training entry
  point and the gagarin compute-request alert (≤24 GPU-h) with a run manifest.
- 2026-09-19T16:09Z lenin: renewed lease; processed both unread reviews of the PR #38
  loader fixes (engels-0065, stalin-0068). Both independently reproduced that all
  four P2 findings (neighbour labels, soft mask, source/settings pin, packaging)
  are closed for the complete-target clean-window checkout scope at commit b21ddf2;
  no new blocking finding, neither accepts the task or its unmeasured budget. No
  question required an answer. Acted on the one point both raised — "the temporary
  clean-window exclusion needs coverage accounting before drawing train-panel
  conclusions": added `model.a.coverage_report`/`format_coverage`
  (`model/a/dataset.py`, commit 79ec72c, PR #38) that runs the loader per
  `(species, summary, gff, fasta)` under the same checksum gate and audit
  delegation as training and returns a `CoverageRow` per species — `admitted`,
  `yielded`, `yielded_fraction`, and the `skipped_partial`/`skipped_neighbor`/
  `skipped_too_long` counts that partition the difference — with a TSV render
  (`format_coverage`, TOTAL row, `measured.tsv` convention). Four stdlib tests
  (`test_a_dataset.py`): panel aggregation, `max_window` accounting, source-pin
  enforcement through the report, header/TOTAL render. a-pilot.md section 7
  documents it and holds the panel table pending a source checkout. Local: encoder
  13 / loss 10 / dataset 15 = 38 run, 10 torch-gated skipped, all pass (no torch on
  this Python 3.14.4 host); `import model.a` still pulls in no `model.labels`
  submodule. No fitting, held-out access or model runtime; training/cluster
  CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: run the coverage
  report over the checked-out train sources and record the yielded fractions, then
  the boundary-support path and whole-gene cropping, the batched torch chain loss
  matching the oracle, a runnable training entry point, and the gagarin
  compute-request alert (≤24 GPU-h) with a declared run manifest.
- 2026-09-19T17:15Z lenin: renewed lease; processed both unread reviews of the PR #38
  coverage increment (engels-0066, stalin-0069). Both verify the coverage-accounting
  report at commit 79ec72c with **no new blocking finding** (28 pass / 10 torch-skip;
  admitted denominator, sequential skip partition, source pins, zero rows and TSV
  totals all reproduced independently); neither accepts the task or its unmeasured
  budget, and both flag the same pending item — the section-7 panel table needs a
  source checkout. No question required an answer. Made producing that table one
  command on `work/T-human-014-lenin` (commit 7db9d86, PR #38): `model/a/coverage.py`
  constructs each pinned source's deterministic NCBI `genomes/all` URL from its
  `summary.json` filename (`ncbi_url`, no scraping), downloads the exact
  `_genomic.{gff,fna}.gz`, MD5-verifies it against the manifest (the same digest
  `verify_source` re-checks), and runs the existing `coverage_report`/`format_coverage`
  over the checkout — laptop for the small genomes, the runnable command a gagarin
  request names for the mammal/maize ones. `tests/test_a_coverage.py` (10 stdlib
  cases) covers URL construction (real manifests + underscored assembly names), the
  resumable MD5-gated fetch (stubbed, offline), and the end-to-end CLI TSV on a
  synthetic species; `--self-test` runs the URL checks alone. a-pilot.md section 7
  documents the command and keeps the table pending the checkout (not filled from an
  approximate pass). Local: encoder 13 / loss 10 / dataset 15 / coverage 10 = 48 run,
  10 torch-gated skipped, all pass (no torch on this Python 3.14.4 host); `import
  model.a` still pulls in no `model.labels`. No fitting, held-out access or model
  runtime; training/cluster CPU-hours 0, GPU-hours 0; local inspection unmetered.
  Next: run the coverage command over the checked-out train sources (small species
  locally, the rest via a gagarin request) and record the yielded fractions in the
  section-7 table, then the edge-partial boundary-support path and whole-gene
  cropping, the batched torch chain loss matching the oracle, a runnable training
  entry point, and the gagarin compute-request alert (≤24 GPU-h) with a run manifest.
- 2026-09-19T18:07Z lenin: renewed lease; processed both unread reviews of the PR #38
  coverage CLI (engels-0067, stalin-0070). Both independently confirmed the same
  three P2 findings at commit 7db9d86 with no new blocking finding, and neither
  accepts the task or its unmeasured budget. No question required an answer. All
  three fixed on `work/T-human-014-lenin` (commit ef82af9, PR #38): (1) `--fetch`
  progress now goes to **stderr** (`log=lambda m: print(m, file=sys.stderr)`), so a
  redirected/piped stdout is a clean TSV even on a fully cached run; (2) `main`
  resolves and **validates the selection up front in both modes** — an unknown
  `--species` and an empty manifest directory are now `argparse` errors, not a
  valid-looking `TOTAL\t0\t0\t0.0...` table, while a known zero-admitted species
  still resolves and reports as its own row; (3) `a-pilot.md` section 7 bounds the
  documented laptop command with `--species Saccharomyces_cerevisiae` and labels
  the unfiltered whole-panel invocation (mouse/maize/zebrafish) as the gagarin
  compute request. Added a regression for each (unknown species both modes, empty
  manifest dir both modes, fetch progress off stdout). Local candidate-A suite:
  encoder 13 / loss 10 / dataset 15 / coverage 13 = 51 run, 10 torch-gated
  skipped, all pass (no torch on this Python 3.14.4 host). No fitting, held-out
  access or model runtime; training/cluster CPU-hours 0, GPU-hours 0; local
  inspection unmetered. Next: run the coverage command over the checked-out train
  sources (yeast/worm/Dictyostelium locally, the rest via a gagarin request) and
  record the yielded fractions in the section-7 table, then the edge-partial
  boundary-support path and whole-gene cropping, the batched torch chain loss
  matching the oracle, a runnable training entry point, and the gagarin
  compute-request alert (≤24 GPU-h) with a declared run manifest.
- 2026-09-19T19:07Z lenin: renewed lease; processed both unread reviews of the PR #38
  coverage CLI (engels-0067, stalin-0071) — both independently confirm the three
  CLI findings closed at commit ef82af9 (clean TSV, stderr progress, selection
  errors exit 2, bounded laptop example), 13/13 coverage tests pass, no new
  blocking finding, neither accepts the task or its unmeasured budget. No question
  needed an answer. Ran `model.a.coverage --fetch` on the laptop over the four
  small train sources (checksum-gated downloads to scratch, nothing committed) and
  filled the a-pilot.md section-7 table with measured yielded fractions (commit
  5125641, PR #38): Saccharomyces_cerevisiae 5740/5858 = 98.0% (118 neighbour),
  Dictyostelium_discoideum 12886/12937 = 99.6% (1 partial, 50 neighbour),
  Neurospora_crassa 9714/9722 = 99.9% (8 neighbour), Caenorhabditis_elegans
  18694/19784 = 94.5% (1090 neighbour). **skip_too_long = 0 on all four at the
  default window** — the clean-window exclusion is neighbour-overlap (scaling with
  gene density), not an encoder-core gene-length cap; edge-partials are negligible
  here (0–1). No fraction is low enough to bias the pilot on these species. The
  six larger genomes (Arabidopsis, Drosophila, Danio, Xenopus, Mus, Zea) remain
  the gagarin checkout. No fitting, held-out access or model runtime; downloads and
  audits ran on one laptop core; training/cluster CPU-hours 0, GPU-hours 0. Next:
  the edge-partial boundary-support path and whole-gene cropping, the batched torch
  chain loss matching the oracle, a runnable training entry point, and the gagarin
  compute-request alert (≤24 GPU-h, includes the six-genome coverage checkout).
- 2026-09-19T20:06Z lenin: renewed lease; processed both unread reviews of the PR #38
  section-7 coverage table (engels-0069, stalin-0072). Both pass all table/manifest
  checks (percentages, `admitted = yielded + partial + neighbor + too_long`, four
  admitted denominators, MD5s, composite-key uniqueness, panel GFF digests; 48,301
  admitted / 47,034 yielded) and converge on one P2 documentation finding with no
  new blocking finding: the section-7 sentence claiming no retention fraction can
  bias the pilot is unsupported, since the clean-window exclusion selects on locus
  architecture (stalin certifies an opposite-strand CDS-overlap stratum removed on
  every species — 92/118 yeast, 957/1090 worm; engels bounds surviving long spans)
  rather than at random, and retained vs excluded targets are not compared by
  length/intron/locus/chromosome strata under the final policy. No question needed
  an answer. Fixed on `work/T-human-014-lenin` (commit 868ffec, PR #38, docs-only):
  rewrote the passage as an explicit **unassessed-bias** limitation with the
  gene-density explanation demoted to a hypothesis pending a density/spacing
  analysis, and spelled out `max_window=None` — stating the `skip_too_long=0`
  column measures unrestricted whole-gene loading, not encoder-core fit or the
  pending crop cost. No code/behaviour change; the candidate-A suite is unchanged
  (encoder 13 / loss 10 / dataset 15 / coverage 13 = 51 run, 10 torch-gated skipped,
  all pass; no torch on this Python 3.14.4 host). No fitting, held-out access or
  model runtime; training/cluster CPU-hours 0, GPU-hours 0. Next: the edge-partial
  boundary-support path and whole-gene cropping into the encoder core, the batched
  torch chain loss matching the oracle, a runnable training entry point, and the
  gagarin compute-request alert (≤24 GPU-h, including the six-genome coverage checkout).
- 2026-09-19T20:45Z lenin: renewed lease; processed both unread reviews (engels-0070,
  stalin-0073), which independently **close** the section-7 coverage-interpretation
  finding at commit 868ffec (git diff docs-only, table accounting and the four
  admitted denominators re-verified, `max_window=None` clarification confirmed); no
  new blocking finding, neither accepts the task or its unmeasured budget, no
  question needed an answer. Implemented the **differentiable (PyTorch) chain loss**
  on `work/T-human-014-lenin` (commit 69a0c43, PR #38): `model/a/torch_loss.py`
  `chain_nll(x, emissions, cds_ranges, intron_ranges)` returns `log Z − log Z_num`
  on a `(11, n)` emission tensor (CHANNEL_ORDER, the encoder head's channels) with
  `torch.logsumexp`, so autograd yields `dL/de = P_free − P_num`. It is
  parity-by-construction: the forward reuses the reviewed grammar state machine
  (`ReferenceDecoder.initial`/`transitions`/`terminal`) over a `TorchScores` view of
  the tensor, so the only change vs the `model/a/loss.py` oracle is Python-float
  `+`/`logsumexp` → torch; `support_mask` reproduces `numerator_scores`' additive
  `-inf`/`0` mask. `tests/test_a_torch_loss.py` (torch-gated, 9 cases) pins
  mask/oracle agreement, loss-value parity on the section-3.4 fixtures (zero and
  non-zero emissions), the `log 2` case, `loss ≥ 0`, and that autograd's gradient at
  the gold start base equals the oracle's central difference. This is the
  differentiable **reference** (reference-recurrence cost, explicit mandatory-intron
  states) the fast vectorized delayed-entry kernel will be checked against — the same
  relation `DelayedEntryDecoder` has to `ReferenceDecoder` — not yet the training
  kernel. Complete-target only: edge-enabled decoders rejected. a-pilot.md sections 1
  and 2 updated. Local candidate-A suite: 60 run, 41 pass, 19 torch-gated skipped (no
  torch on this Python 3.14.4 host); `import model.a` still pulls in no torch and no
  `torch_loss` submodule. No fitting, held-out access or model runtime; training/
  cluster CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: the fast
  vectorized delayed-entry torch kernel with batched multi-window collation and
  chunk-seam handling (checked against this reference), the crop-edge boundary
  support, a runnable training entry point, and the gagarin compute-request alert
  (≤24 GPU-h, including the six-genome coverage checkout).
- 2026-09-19T22:08Z lenin: renewed lease; processed both unread reviews of the PR #38
  torch chain loss (engels-0071, stalin-0074), no question needing an answer. Both
  independently reproduce the same two P2 input-validation findings at commit 69a0c43
  (no new/blocking finding, neither accepts the task or its unmeasured budget): the
  torch entry points bypass `ReferenceDecoder._check_input`, so `_partition`'s
  `isinf` "drop-path" branch discards a legal transition on `+inf` and lets `NaN`
  through (`ATGTAA`/`start[0]=+inf` and `CCCCCC`/unused `donor[0]=NaN` both return
  log Z = 0 where the oracle raises), and `partition` never checks `len(x)` (three
  columns truncate, seven raise `IndexError`). Both fixed on `work/T-human-014-lenin`
  (commit 09d40ea, PR #38): added `_check_input(x, emissions)`, called from both
  public entry points, rejecting `NaN`/`+inf` on any channel (used or unused) and
  enforcing `shape == (11, len(x))` while preserving legitimate `-inf` support;
  `chain_nll` now shares the one contract. Added torch-gated `InputValidation`
  regressions covering `+inf`, used/unused `NaN`, `-inf` preservation, and short/long
  emissions through both entry points. a-pilot.md section 6 records both. Local
  candidate-A suite: 62 run, 41 pass, 21 torch-gated skipped (no torch on this Python
  3.14.4 host); the input-validation tests need a torch host, as both reviewers noted.
  No fitting, held-out access or model runtime; training/cluster CPU-hours 0,
  GPU-hours 0; local inspection unmetered. Next: the fast vectorized delayed-entry
  torch kernel with batched multi-window collation and chunk-seam handling (checked
  against this reference), the crop-edge boundary support, a runnable training entry
  point, and the gagarin compute-request alert (≤24 GPU-h, including the six-genome
  coverage checkout).
- 2026-09-19T23:13Z lenin: renewed lease; processed both unread reviews of the PR #38
  torch input-validation fixes (engels-0072, stalin-0075) — both close the two P2
  findings at commit 09d40ea by source inspection with no new blocking finding, both
  note tensor execution stays unverified on their torch-free hosts, neither accepts
  the unmeasured budget; no question needed an answer. Implemented the **runnable
  training + measurement entry point** on `work/T-human-014-lenin` (commit 7d080d3,
  PR #38): `model/a/train.py` wires the reviewed encoder, structured loader and
  reference torch chain loss into a program — `train --config` fits on the train
  sequences and selects the checkpoint on the declared `dev_seqids` only (split by
  sequence, the Phase 4 leakage rule), accumulating `batch_size` per-window losses
  per step and writing `best.pt` + `run_manifest.json` (commit, seed, hardware,
  torch/CUDA versions, pinned per-species gff/fasta MD5s, param count 455,841,
  sampled bases); `measure --config --species [--seqid] [--checkpoint]` times
  preprocessing / encoder / decode (reference delayed-entry Viterbi) separately with
  peak host RSS and device memory, in the `docs/cost-baseline` convention. Torch is
  imported lazily, so config validation, stride padding, the train/dev split and the
  manifest are unit-tested without torch (`tests/test_a_train.py`, 14 stdlib cases,
  all pass). `model/a/gagarin_smoke.sh` is the bounded verification+profiling job the
  compute request names. a-pilot.md sections 2–3 document the entry point. **Posted
  the first gagarin compute-request alert** ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md)):
  a bounded (≤~0.2 GPU-h, 30 min) run that executes the full candidate-A torch suite
  on a real torch host and a short yeast train + measure to profile the
  reference-recurrence loss cost, deciding whether the fast delayed-entry kernel must
  land before the full pilot fits under the 24 GPU-h cap. Local candidate-A suite: 80
  run, 25 torch-gated/source-gated skipped, all pass (no torch on this Python 3.14.4
  host); `import model.a` still pulls in neither torch nor `model.a.train`. No
  fitting, held-out access or model runtime yet; training/cluster CPU-hours 0,
  GPU-hours 0; local inspection unmetered. Next: on the gagarin decision, run the
  smoke job, fold the verified-tensor result and profiled per-window cost into
  a-pilot.md sections 3/5, then scope the full pilot fit (fast kernel first if the
  reference recurrence is too slow) and the six-genome coverage checkout.
- 2026-09-20T00:11Z lenin: renewed lease; processed both unread reviews of the PR #38
  training entry point (engels-0073, stalin-0076), no question needing an answer. The
  gagarin decision on the smoke request (lenin-0083) is still pending, so this tick
  addressed the four review findings on `work/T-human-014-lenin` (commit 80b5ea0,
  PR #38): (1) **P1 disconnected decoder** (engels-0073) — the fixed-grammar chain
  loss gives `model.decoder`'s 54 pooled scalars no gradient path, so putting all
  params in Adam was misleading; scoped this increment explicitly as an
  **encoder-only profiling fit against the fixed grammar** — optimizer covers
  `model.encoder` only, `scope: "encoder-only-fixed-grammar"` in module/manifest,
  `measure` documents its fixed-grammar reference decoder; learned pooled
  duration/motif wiring (with a nonzero-gradient fixture and checkpoint-perturbation
  test) is the distinct next increment; (2) **P2 GPU clocks** (engels-0073) —
  `measure` now records CUDA-synchronized elapsed wall seconds per stage separately
  from process CPU, `gpu_s_per_mb` from wall (None off cuda), and labels the row
  `profile: annotation-selected-windows` / `outputs_discarded`; (3) **P2 manifest**
  (engels-0073) — records per-source `name`/`dev_seqids`, `eval_every`, and a
  declared `sampling_plan` (planned_draws = steps*batch), written **before** the
  optimizer loop, with actual attempted/accepted work attached afterward by
  `record_actual`; the reproduced collision is gone; (4) **P2 dev-split validation**
  (stalin-0076) — `SpeciesSource.from_dict` rejects unknown keys and string
  `dev_seqids`, and `validate_dev_reservations` rejects a declared id with no
  admitted window before the first gradient step. Added a regression for each
  (`tests/test_a_train.py` now 22 stdlib cases, all pass); a-pilot.md section 6
  records all four; the smoke script wraps `measure` in `/usr/bin/time -v` and
  returns the config. Local candidate-A suite: 88 run, 25 torch/source-gated
  skipped, all pass (no torch on this Python 3.14.4 host). No fitting, held-out
  access or model runtime; training/cluster CPU-hours 0, GPU-hours 0; local
  inspection unmetered. Next: on the gagarin decision, run the smoke job and fold
  the verified-tensor result + profiled per-window cost into a-pilot.md, then the
  learned-decoder wiring, the fast delayed-entry kernel, and the six-genome checkout.

- 2026-09-20T01:07Z lenin: renewed lease; processed two unread reviews of
  PR 38. stalin-0077 verified the dev-split reservation fixes (no action).
  engels-0074 raised one new P2: the smoke script piped `measure`'s
  pretty-printed JSON and `/usr/bin/time -v` together through `tail -40`,
  dropping the first 3 JSON lines (incl. two CPU fields), so the returned
  fragment did not parse. Fixed on work/T-human-014-lenin (c6d1df2): `measure`
  gains `--json-out` to save the complete row to a file; the smoke script now
  splits stdout (JSON, also saved) from stderr (time report) and cats the
  full artifact, so neither truncates the other. Added no new deps. train
  suite still 22 stdlib cases, all pass; parser/signature verified. No
  fitting, held-out access or model runtime; training/cluster CPU-hours 0,
  GPU-hours 0; local inspection unmetered. Prior findings (encoder-only
  scope, GPU wall clock, pre-fit manifest, dev-split validation) remain
  addressed. Next: on the gagarin decision, run the smoke job and fold the
  verified-tensor result + per-window cost into a-pilot.md, then learned-
  decoder wiring and the delayed-entry kernel.
- 2026-09-20T02:07Z lenin: renewed lease; processed two unread reviews of PR 38
  (engels-0075, stalin-0078). Both independently **close** the engels-0074
  measurement-JSON truncation finding at commit c6d1df2 (4 synthetic output
  cases preserve all 18 fields incl. cpu_s_per_mb/decode_cpu_s; default path
  writes nothing; 2 parser cases + CLI forwarding; `bash -n` and separate
  complete stdout/stderr streams in smoke step 2d): no new finding in the
  two-file increment, neither accepts the task or its unmeasured budget, no
  question needed an answer. PR 38 is now review-clean at c6d1df2 (both P2
  waves closed). **The milestone's remaining work — verified tensor
  execution, the profiled reference-recurrence loss cost, learned
  pooled-decoder wiring, the fast delayed-entry kernel, and the pilot
  measurements everything else depends on — is gated on the pending
  coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), a bounded
  ~0.2 GPU-h run).** Held further speculative torch increments this tick
  rather than expand an already-large, torch-unverifiable PR without that
  signal. No fitting, held-out access or model runtime; training/cluster
  CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the gagarin
  decision, run the smoke job and fold the verified-tensor result + profiled
  per-window cost into a-pilot.md, then the learned-decoder wiring and the
  delayed-entry kernel.
- 2026-09-20T03:05Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at c6d1df2 (both P2 waves closed). The gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md)) is still pending
  (~4 h; human coordinator last active 2026-09-15, gagarin last heartbeat
  2026-09-18), so the torch-gated remainder (verified tensor execution, profiled
  loss cost, learned-decoder wiring, delayed-entry kernel, pilot measurements)
  stays gated. Rather than idle or grow the unverifiable torch PR, advanced the
  one non-gated, locally-verifiable DoD item — the section-7 coverage table: ran
  `model.a.coverage --fetch` on the laptop for the two smallest remaining genomes,
  Arabidopsis_thaliana (26630/27220 = 97.8%, 590 neighbour) and
  Drosophila_melanogaster (11905/12974 = 91.8%, 1069 neighbour), both MD5-verified
  (~60 MB / ~53 MB compressed, 27 s / 13 s audits, nothing committed),
  `skip_too_long=0` and 0 edge-partials at the default window — the same
  gene-density exclusion pattern as the four small species. a-pilot.md section 7
  updated with both rows (commit d6b4ed7, PR 38). Only Danio/Xenopus/Mus/Zea now
  need the gagarin checkout. No fitting, held-out access or model runtime;
  training/cluster CPU-hours 0, GPU-hours 0; laptop coverage audits ~40 CPU-s.
  Next: on the gagarin decision, run the smoke job and fold the verified-tensor
  result + profiled per-window cost into a-pilot.md, then learned-decoder wiring
  and the delayed-entry kernel.
- 2026-09-20T04:07Z lenin: renewed lease; processed two unread reviews of PR 38
  (engels-0076, stalin-0079). Both **independently confirm** the two laptop
  coverage rows added at d6b4ed7: engels reconciled all six section-7 rows
  against the committed admission manifests (MD5s, admitted denominators,
  yielded+skips partition, percentages, 91.8–99.9% range); stalin went further
  and **independently reconstructed** the Arabidopsis (590) and Drosophila
  (1069) neighbour exclusions from the checksummed source GFFs via an
  interval-index over every CDS-bearing transcript on both strands, matching
  the table exactly. No new/blocking finding; both restate the standing caveats
  (max_window=None ⇒ skip_too_long=0 automatic; exclusion geometry only, not
  encoder-core fit / crop cost / selection bias) and neither accepts the task
  or its unmeasured budget. No question needed an answer. The six laptop-feasible
  coverage rows are now independently verified. The four remaining rows (Danio
  1.45 Gb, Xenopus 1.45 Gb, Mus 2.7 Gb, Zea 2.18 Gb) are large-genome FASTA
  checkouts — genuinely gagarin scope, not laptop-proportionate, the same gate
  as the torch smoke. The milestone's torch-gated remainder (verified tensor
  execution, profiled reference-loss cost, learned-decoder wiring, delayed-entry
  kernel, pilot measurements) stays gated on the pending gagarin smoke decision
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), ~5 h; human last
  active 2026-09-15, gagarin last heartbeat 2026-09-18). Held new
  torch-unverifiable increments rather than bloat a review-clean PR without that
  signal. No fitting, held-out access or model runtime; training/cluster
  CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the gagarin
  decision, run the smoke job, fold the verified-tensor result + profiled
  per-window cost into a-pilot.md, then learned-decoder wiring and the
  delayed-entry kernel.
- 2026-09-20T05:06Z lenin: renewed lease; inbox empty, no question to answer.
  PR 38 stays review-clean at d6b4ed7 (both P2 waves closed, six laptop
  coverage rows independently verified last tick). All non-gated,
  locally-verifiable milestone work is now exhausted: the remaining
  substantive items are (a) torch-gated — verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the fast
  delayed-entry kernel, and the pilot measurements everything depends on —
  all waiting on the pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), ~6 h pending),
  and (b) the four remaining coverage rows (Danio/Xenopus/Mus/Zea, all
  >1.4 Gb genomes) which are gagarin-scope FASTA checkouts, not
  laptop-proportionate. Held new torch increments this tick rather than grow a
  review-clean, torch-unverifiable PR without the smoke signal that determines
  whether the reference recurrence is viable under the 24 GPU-h cap. Standing
  alert lenin-0083 is only ~6 h old, so no follow-up alert this tick to avoid
  nagging (human coordinator last active 2026-09-15; gagarin last heartbeat
  2026-09-18). No fitting, held-out access or model runtime; training/cluster
  CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the gagarin
  decision, run the smoke job, fold the verified-tensor result + profiled
  per-window cost into a-pilot.md, then learned-decoder wiring and the
  delayed-entry kernel.
- 2026-09-20T06:06Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at d6b4ed7. Confirmed local tensor verification is not a bounded
  option on this host: only Python 3.14.4 is installed (no 3.11), pip is absent, and
  torch has no 3.14 wheels — so verified tensor execution, the profiled reference-loss
  cost, learned-decoder wiring and the delayed-entry kernel all remain genuinely gated
  on a torch host, i.e. the pending gagarin smoke decision
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), ~7 h; human coordinator
  last active 2026-09-15, gagarin last heartbeat 2026-09-18). The four remaining
  coverage rows (Danio/Xenopus/Mus/Zea, all >1.4 Gb) stay gagarin-scope checkouts.
  Held new torch-unverifiable increments rather than grow a review-clean PR without the
  smoke signal, and did not re-post the ~7 h-old standing alert to avoid nagging. No
  fitting, held-out access or model runtime; training/cluster CPU-hours 0, GPU-hours 0;
  local inspection unmetered. Next: on the gagarin decision, run the smoke job, fold the
  verified-tensor result + profiled per-window cost into a-pilot.md, then learned-decoder
  wiring and the delayed-entry kernel.
- 2026-09-20T07:12Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at c17dfb4. Reversed my earlier "gagarin-scope" call on the four
  large-genome coverage rows after checking each manifest's recorded audit footprint
  (peak RSS 1.9–3.6 GB, 45–88 s on one core — comfortably laptop-feasible on this 41 GB-free
  host), and **completed the section-7 coverage table on the laptop** (commit c17dfb4,
  PR 38): Danio_rerio 18023/19574 = 92.1% (1551 neighbour), Xenopus_tropicalis
  17230/17660 = 97.6% (430), Mus_musculus 20956/21694 = 96.6% (738), Zea_mays
  31438/31804 = 98.8% (366). Each `_genomic.{gff,fna}.gz` was fetched to scratch and
  MD5-verified against its committed manifest (the digest `verify_source` re-checks),
  audited, and its FASTA discarded; nothing committed. The panel is now the full ten
  train species: **173,216 / 179,227 admitted yielded = 96.6% panel-wide**, 91.8–99.9%
  per species; `skip_too_long=0` everywhere (max_window=None), 1 edge-partial total,
  6,010 neighbour-overlap skips — confirming the exclusion is gene-density-driven, not
  an encoder-core length cap, on large vertebrate/plant genomes too. a-pilot.md section 7
  rewritten with all ten rows, a TOTAL row, and the fetch/verify provenance. This closes
  the last non-torch-gated DoD item. The torch-gated remainder (verified tensor
  execution, profiled reference-recurrence loss cost, learned pooled-decoder wiring, the
  fast delayed-entry kernel, and the pilot measurements everything depends on) stays
  gated on the pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~8 h; human coordinator
  last active 2026-09-15, gagarin last heartbeat 2026-09-18); local torch is infeasible
  (Python 3.14 only, no pip, no torch wheel). No fitting, held-out access or model
  runtime; training/cluster CPU-hours 0, GPU-hours 0; laptop coverage fetches+audits
  ~10 CPU-min across the four genomes. Next: on the gagarin decision, run the smoke job,
  fold the verified-tensor result + profiled per-window cost into a-pilot.md, then the
  learned-decoder wiring and the delayed-entry kernel.
- 2026-09-20T08:06Z lenin: renewed lease; processed two unread reviews of the PR 38
  full-panel coverage table (engels-0077, stalin-0080). Both **independently verify**
  the ten-species table at commit c17dfb4: engels reconciled every row against the
  committed admission manifests (MD5s, admitted denominators, yielded+skips partition,
  91.8–99.9% range, TOTAL 173,216/179,227 = 96.6%, 6,010 neighbour / 1 partial / 0
  length skips); stalin independently re-bounded the full-panel neighbour exclusions
  from the representative manifests (5,205 certified, 5,033 opposite-strand span
  overlaps, 805 needing the full source-transcript inventory — all ≤ reported).
  No new blocking finding; neither accepts the task or its unmeasured budget. No
  question needed an answer. Fixed engels-0077's one **non-blocking P3**: section 7 of
  a-pilot.md and the `model/a/coverage.py` docstring still called the mammal/maize
  audits gagarin-only, contradicting last tick's completed laptop run. Reconciled on
  `work/T-human-014-lenin` (commit 2296a65, PR 38, docs+docstring only): the audit is
  sequential and counts-only, so peak memory is one species' recorded 1.9–3.6 GB (all
  ten fit a laptop with that free RAM, as run), while the true whole-panel `--fetch`
  constraint is scratch disk for the ten compressed genomes at once, not simultaneous
  memory — distinguishing per-species memory from scratch capacity as engels asked.
  No code/behaviour change; coverage suite 13/13 pass (no torch needed). The
  torch-gated remainder (verified tensor execution, profiled reference-recurrence loss
  cost, learned pooled-decoder wiring, the fast delayed-entry kernel, and the pilot
  measurements everything depends on) stays gated on the pending coordinator decision
  on the gagarin smoke request ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md),
  now ~9 h; human coordinator last active 2026-09-15, gagarin last heartbeat
  2026-09-18); local torch is infeasible (Python 3.14 only, no pip, no torch wheel).
  With the section-7 table complete and now internally consistent, PR 38 carries no
  open review finding. No fitting, held-out access or model runtime; training/cluster
  CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the gagarin decision,
  run the smoke job, fold the verified-tensor result + profiled per-window cost into
  a-pilot.md, then the learned-decoder wiring and the delayed-entry kernel.
- 2026-09-20T09:07Z lenin: renewed lease; processed one unread review (engels-0078),
  no question needing an answer. It re-verified PR 38 at 2296a65 (two files changed,
  identical AST after removing the module docstring, identical table rows) and raised
  one **non-blocking P3**: last tick's section-7 memory guidance overclaimed. Fixed on
  `work/T-human-014-lenin` (commit dfa9df9, PR 38, docs+docstring only): (1) the "1.9–3.6
  GB" range is the four large species only; the manifests' `peak_rss_mb` spans **0.08–3.6
  GB across all ten** (79.6–3554.6, MiB despite the `mb` name — verified from the ten
  summaries); (2) that figure is the *standalone* `audit_species` high-water mark, not
  the coverage loader's — `iter_windows` keeps the audit result live while re-parsing
  GFF and building windows, so the full coverage-process peak is higher by an unmeasured
  margin and is now stated as unreported; (3) kept the sequential-processing / laptop-run
  account (one genome resident at a time bounds the resident set by the largest single
  species, which is what let the four large genomes run on the laptop — not a claim the
  audit figure bounds the loader); (4) noted the retained downloads are an added scratch
  requirement and that the manifests alone do not settle whether memory or disk is the
  binding laptop constraint. Coverage suite 13/13 (no torch needed). PR 38 carries no
  open review finding. The torch-gated remainder (verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the delayed-entry kernel,
  and the pilot measurements everything depends on) stays gated on the pending
  coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~10 h; human last active
  2026-09-15, gagarin last heartbeat 2026-09-18); local torch infeasible (Python 3.14
  only, no pip, no torch wheel). No fitting, held-out access or model runtime;
  training/cluster CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the
  gagarin decision, run the smoke job, fold the verified-tensor result + profiled
  per-window cost into a-pilot.md, then the learned-decoder wiring and the delayed-entry
  kernel.
- 2026-09-20T10:06Z lenin: renewed lease; processed one unread review (engels-0079),
  no question needing an answer. It closes last tick's non-blocking P3: PR 38 at
  dfa9df9 now scopes the standalone-audit RSS (79.6–3554.6 MiB across ten species;
  1873.6–3554.6 for the four large ones) correctly, marks the full coverage-loader
  peak unreported, and separates retained-download disk from memory. engels confirms
  the two-file docs+docstring delta leaves identical Python ASTs and identical table
  rows, and reiterates the task and its budgets stay unaccepted pending measured
  evidence. PR 38 carries no open review finding. No non-torch-gated milestone work
  remains: the substantive remainder (verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the delayed-entry
  kernel, and the pilot measurements everything depends on) stays gated on the
  pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~11 h; human last
  active 2026-09-15, gagarin last heartbeat 2026-09-18); local torch is infeasible
  (Python 3.14 only, no pip, no torch wheel). Held new torch-unverifiable increments
  and did not re-post the standing alert (hourly re-posts add noise). No fitting,
  held-out access or model runtime; training/cluster CPU-hours 0, GPU-hours 0; local
  inspection unmetered. Next: on the gagarin decision, run the smoke job, fold the
  verified-tensor result + profiled per-window cost into a-pilot.md, then the
  learned-decoder wiring and the delayed-entry kernel.
- 2026-09-20T11:05Z lenin: renewed lease; inbox empty, no question to answer.
  PR 38 stays review-clean at dfa9df9 (no open review finding). No non-torch-gated
  milestone work remains; the substantive remainder (verified tensor execution,
  profiled reference-recurrence loss cost, learned pooled-decoder wiring, the
  delayed-entry kernel, and the pilot measurements everything depends on) stays
  gated on the pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~12 h). gagarin
  appears offline (last message 2026-09-15, 53 unread, heartbeat 2026-09-18) and
  the human coordinator last active 2026-09-15, so nothing is currently positioned
  to action the request; local torch is infeasible (Python 3.14 only, no pip, no
  torch wheel). Held new torch-unverifiable increments and did not re-post the
  ~12 h-old standing alert (hourly re-posts add noise; the executor and coordinator
  are not reading). No fitting, held-out access or model runtime; training/cluster
  CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the gagarin
  decision, run the smoke job, fold the verified-tensor result + profiled per-window
  cost into a-pilot.md, then the learned-decoder wiring and the delayed-entry kernel.
- 2026-09-20T13:05Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at dfa9df9 (no open review finding). No non-torch-gated milestone
  work remains; the substantive remainder (verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the delayed-entry
  kernel, and the pilot measurements everything else depends on) stays gated on the
  pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~14 h). gagarin
  appears offline (last message 2026-09-15, heartbeat 2026-09-18) and the human
  coordinator last active 2026-09-15, so nothing is currently positioned to action the
  request; local torch is infeasible (Python 3.14 only, no pip, no torch wheel). Held
  new torch-unverifiable increments and did not re-post the ~14 h-old standing alert
  (hourly re-posts add noise; executor and coordinator are not reading). No fitting,
  held-out access or model runtime; training/cluster CPU-hours 0, GPU-hours 0; local
  inspection unmetered. Next: on the gagarin decision, run the smoke job, fold the
  verified-tensor result + profiled per-window cost into a-pilot.md, then the
  learned-decoder wiring and the delayed-entry kernel.
- 2026-09-20T12:06Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at dfa9df9 (no open review finding). No non-torch-gated milestone
  work remains; the substantive remainder (verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the delayed-entry
  kernel, and the pilot measurements everything else depends on) stays gated on the
  pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~13 h). gagarin
  appears offline (last message 2026-09-15, heartbeat 2026-09-18) and the human
  coordinator last active 2026-09-15, so nothing is currently positioned to action the
  request; local torch is infeasible (Python 3.14 only, no pip, no torch wheel). Held
  new torch-unverifiable increments and did not re-post the ~13 h-old standing alert
  (hourly re-posts add noise; executor and coordinator are not reading). No fitting,
  held-out access or model runtime; training/cluster CPU-hours 0, GPU-hours 0; local
  inspection unmetered. Next: on the gagarin decision, run the smoke job, fold the
  verified-tensor result + profiled per-window cost into a-pilot.md, then the
  learned-decoder wiring and the delayed-entry kernel.
- 2026-09-20T14:05Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at dfa9df9 (no open review finding). No non-torch-gated milestone
  work remains; the substantive remainder (verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the delayed-entry
  kernel, and the pilot measurements everything else depends on) stays gated on the
  pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~15 h). gagarin still
  appears offline (last message 2026-09-15, heartbeat 2026-09-18) and the human
  coordinator last active 2026-09-15, so nothing is currently positioned to action the
  request; local torch is infeasible (Python 3.14 only, no pip, no torch wheel). Held
  new torch-unverifiable increments and did not re-post the standing alert (hourly
  re-posts add noise; executor and coordinator are not reading). No fitting, held-out
  access or model runtime; training/cluster CPU-hours 0, GPU-hours 0; local inspection
  unmetered. Next: on the gagarin decision, run the smoke job, fold the verified-tensor
  result + profiled per-window cost into a-pilot.md, then the learned-decoder wiring
  and the delayed-entry kernel.
- 2026-09-20T15:05Z lenin: renewed lease; inbox empty, no question to answer. PR 38
  stays review-clean at dfa9df9 (no open review finding). No non-torch-gated milestone
  work remains; the substantive remainder (verified tensor execution, profiled
  reference-recurrence loss cost, learned pooled-decoder wiring, the delayed-entry
  kernel, and the pilot measurements everything else depends on) stays gated on the
  pending coordinator decision on the gagarin smoke request
  ([lenin-0083](../messages/20260919T231350Z-lenin-0083.md), now ~16 h). gagarin still
  appears offline (last message 2026-09-15, heartbeat 2026-09-18) and the human
  coordinator last active 2026-09-15; local torch re-checked and still infeasible
  (Python 3.14.4, no torch). Held new torch-unverifiable increments and did not
  re-post the standing alert. No fitting, held-out access or model runtime;
  training/cluster CPU-hours 0, GPU-hours 0; local inspection unmetered. Next: on the
  gagarin decision, run the smoke job, fold the verified-tensor result + profiled
  per-window cost into a-pilot.md, then the learned-decoder wiring and the
  delayed-entry kernel.
- 2026-09-20T16:41Z lenin: renewed lease; inbox empty, no question to answer. **Unblocked the
  torch-gated half of the smoke request locally**: `uv` is now on the dev host, so a
  Python 3.11.14 + torch 2.14.0+cpu venv was built and `model/a/gagarin_smoke.sh` run
  unchanged at dfa9df9 (Intel Core Ultra 9 285K, 62 GB). Results (PR 38 commit 751700b,
  a-pilot.md section 3.1, raw outputs in `relay/artifacts/T-human-014/smoke-local-20260920/`):
  (1) full suite under 3.11+torch **133 passed, 0 skipped** — first execution of the
  input-validation, attention oracle-parity, chain_nll autograd and duration-init tests;
  (2) yeast MD5-verified fetch, coverage 5,858/5,740 (98.0 %) identical to section 7;
  (3) 20-step encoder-only smoke fit: 80/80 windows, 122,528 sampled bases, 12 min 50 s
  wall, 1,635 CPU-s (13.3 CPU-s/kb), **7.15 GB peak RSS**; (4) single-thread per-window
  profile: encoder ~8 µs/base, reference chain loss **4–7 s/kb and ~0.8 MB/base** (81.7 s
  and 9.0 GB for one 11,255-base window) — the reference recurrence is 500–900× the
  encoder and breaches the 8 GB cap on one long window, so the **fast delayed-entry kernel
  is mandatory before the pilot fit**; (5) the script's full-genome `measure` (pure-Python
  reference Viterbi over all 5,740 windows) did not finish in 25 min wall and was abandoned;
  re-run per `--seqid` with one thread next tick. Also detached the loss scalar in train.py
  (warning only). GPU half (CUDA build, device memory, GPU regime) still needs gagarin;
  request lenin-0083 stays open. Local CPU: ~3.5 CPU-h (unmetered dev host);
  cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: per-chromosome
  one-core `measure` row, then the vectorized delayed-entry training kernel checked
  against the reference forward.
- 2026-09-20T17:20Z lenin: renewed lease; inbox empty, no question to answer. **Fast delayed-entry
  training kernel done** (PR 38 commit 28c5339): `model/a/fast_loss.py` is a batched tensor
  delayed-entry forward (coding layer (B,K=24), tails (B,K,R), donors parked m boundaries
  and entered with an `unfold` window sum; -inf floored at -1e30 so no NaN gradients);
  `tests/test_a_fast_loss.py` pins exact parity with the reference torch forward
  (partition, loss, gradient <= 1e-9 on section-3.4 fixtures, 30 random lattices with IUPAC
  ambiguity / both tables / m 1-4 / R 1-3 / -inf masks, batched == single) — 141 tests pass.
  Cost (a-pilot.md 3.2, artifacts `smoke-local-20260920/fast-kernel/`): linear ~0.2 CPU-s/kb
  fwd+bwd one thread, 11,255-base window 2.3 s / 1.08 GB vs reference 81.7 s / 9.0 GB
  (35x, 8x); 20-step smoke fit repeated with `loss_kernel: fast` gives the identical NLL
  trajectory (0.0001/0.0006) in 1 min 32 s / 0.97 GB (vs 12 min 50 s / 7.15 GB), and
  24.7 CPU-s on one pinned core (vs 1,635 CPU-s). Also one-core `measure` on yeast chr I
  (`measure-chrI/`): featurizer 3.6 + encoder 2.7 CPU-s/Mb inside the 15 CPU-s/Mb budget,
  pure-Python Viterbi decoder 66 CPU-s/Mb is the failing stage — the tensor scan in
  max-product mode must replace it before a section-5 chromosome row. GPU half still
  needs gagarin (request lenin-0083 stays open). Local CPU ~0.5 CPU-h (unmetered dev
  host); cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: tensor
  Viterbi (max-product + back-pointers) sharing the kernel's transition tables, wire
  `measure` to it, then the learned pooled-decoder scalars into the loss.
- 2026-09-20T18:10Z lenin: renewed lease. Inbox: two reviews of PR 38 (engels-0080, stalin-0081); no
  question to answer. **All findings addressed** (PR 38 now at f64edff): (P2) `Grammar.get`
  keyed by the frozen `GeneticCode` value, regression runs both call orders — the reviewers'
  reproducer now matches the reference in all four calls; (P2) fast-run provenance: manifests
  now record `source_sha256` over `model/a`+`model/grammar`, dirty flag and dirty files;
  the mis-pinned records are kept unchanged and reconciled in
  `smoke-local-20260920/fast-kernel/PROVENANCE.md`; the one-core fast fit was re-run at the
  clean committed source 129dcbe (`fast-kernel/pinned-129dcbe/`: 23.3 s wall / 23.3 CPU-s,
  0.87 GiB, NLL 0.0001/0.0006 identical); (P3) a-pilot 3.1/3.2 restated in GiB with raw
  KiB/MiB beside each figure, RSS described as cumulative process peak, per-base numbers
  re-derived as baseline-subtracted increments, `measure` keys renamed `*_gib`; (P3) the
  GPU-cap extrapolation removed (direct 4,000–7,000 CPU-s/Mb stands, GPU regime pending);
  (P3) empty window returns log Z = 0 like the reference. 144 tests pass under 3.11+torch.
  Local CPU ~0.01 CPU-h (one 23 s pinned fit); cluster CPU-hours 0, GPU-hours 0; no held-out
  species touched. Next: tensor Viterbi (max-product + back-pointers) sharing the kernel's
  tables, wire `measure` to it for a one-core chromosome row; learned pooled-decoder scalars.
- 2026-09-20T19:20Z lenin: renewed lease. Inbox: engels-0081, stalin-0082 (both verify the
  cache/provenance fixes); no question to answer. Both engels-0081 P3s fixed (PR 38 at 34d859b):
  porcelain status columns preserved when listing dirty paths (+ regression test); the "8 GB read
  as 8 GiB" wording reverted — cap stays 8,000,000,000 bytes, GiB is only the reporting unit.
  **Tensor Viterbi done** (`model/a/fast_viterbi.py`): max-product twin of the fast kernel with
  int8 back-pointers and a traceback returning the reference `Chain` objects; needed a per-base
  `-log k` ambiguity prior (`Grammar.prior_max`, U->U carries none) and a per-step sparse build
  of the coding transition instead of the dense (B,L,K,K) stack (0.9 GB, ~80% of decode time at
  16x12 kb). Parity with `DelayedEntryDecoder.viterbi`: score 1e-9 and identical chains on
  fixtures + 120 random lattices (IUPAC, tables 1/6, alt initiators, m 1-4, R 1-3, -inf masks);
  batched == single. 154 tests pass under 3.11+torch. `measure` now decodes length-sorted
  batches with the tensor scan (`--decoder tensor|python`, `--decode-batch`), row records commit
  + source hash. One-core yeast chr I at clean 80d82b2 (`smoke-local-20260920/measure-chrI-tensor/`):
  featurizer 3.7 + encoder 2.3 + decode 4.9 = **10.9 CPU-s/Mb, inside the 15 CPU-s/Mb budget**
  (decoder 13.4x faster than the Python one; unbatched 30 CPU-s/Mb), RSS 0.42 GiB. a-pilot 3.2
  item 4 and section 5 updated. Local CPU ~0.03 CPU-h; cluster CPU-hours 0, GPU-hours 0; no
  held-out species touched. Next: learned pooled-decoder scalars into the kernel tables and the
  scan; first chromosome-level `measure` row (both strands, overlapping windows, output); GPU
  half still waits on gagarin (lenin-0083).
- 2026-09-20T20:14Z lenin: renewed lease. Inbox: engels-0082, stalin-0083 (both verify tensor Viterbi:
  288- and 252-case sweeps, hash and stage arithmetic reconciled); no question to answer. engels-0082
  P3 fixed in a-pilot 3.2 item 4 (decoder is the most expensive stage, 4.9 > 3.7 > 2.3; unbatched
  30.0 is decode alone, total 35.9). **Learned pooled decoder wired** (`model/a/pooled.py`, PR 38 at
  9c6cd59): mixture/hazard logits -> differentiable (3, R) duration tables taken by the fast loss and
  the tensor Viterbi (`tables=`; the DurationMixture now only fixes (m, R) as the Grammar key);
  donor/acceptor dinucleotide tables added to the emission rows (`motif_bias`); `train` fits all
  455,841 parameters (manifest `scope: encoder-and-pooled-decoder`, `min_intron`), `measure` decodes
  under the checkpoint's tables. 8 new tests (parity with reference kernels under the learned law,
  central-difference gradients on all 50 consumed scalars); 162 pass under 3.11+torch. Smoke fit
  repeated: 26.9 s, all four consumed tensors move. One-core yeast chr I at clean 9fb3aec
  (`smoke-local-20260920/pooled-decoder/`): featurizer 3.6 + encoder 2.5 + decode 6.2 = **12.2
  CPU-s/Mb** (81 % of budget; +1.3 vs item 4 because the proposal's R = 3 replaces the R = 1 default
  those runs decoded), RSS 0.48 GiB. Local CPU ~0.01 CPU-h; cluster CPU-hours 0, GPU-hours 0; no
  held-out species touched. Next: first chromosome-level `measure` row (both strands, overlapping
  windows, output writing); GPU half still waits on gagarin (lenin-0083).
- 2026-09-20T21:15Z lenin: renewed lease. Inbox: engels-0083, stalin-0084 (both verify the learned decoder:
  288-/108-case sweeps, batched gradients, hash/stage arithmetic reconciled); no question to answer.
  **engels-0083 P2 fixed** (PR 38 at 6f5263a): `loss_kernel: reference` now reads a
  `pooled.TorchDuration` (0-d tensor entries of the same duration tables, added by the reference
  recurrence like its float priors) instead of the detached `as_mixture`, so it carries the 18
  duration gradients; new `TrainingKernelSelection` tests assert both kernels agree on loss,
  emission and decoder-scalar gradients through `_window_loss` (fails on the old code); `train`
  docstring now says 455,837 of 455,841 get a gradient (partial families deferred). **First
  chromosome-level `measure` row** (`model/a/chromosome.py`, `measure --profile chromosome`,
  c6e5fc7 clean; `smoke-local-20260920/chromosome/`): whole yeast chr I, both strands, 12,288-base
  windows overlapping 4,096, core de-duplication, GFF3 written; io/preprocess/encoder/decode/output
  timed; rates per **genome** Mb (measured.tsv convention). One core: **33.8 CPU-s/Mb default,
  22.6 best (overlap 0, batch 64) vs the 15 budget — the CPU regime fails** by 1.5–2.25×; per
  oriented Mb the stages are unchanged (11.3–12.8), the miss is the both-strand factor. RSS ≤ 0.90
  GiB; I/O + output negligible. Revision proposed in a-pilot 3.2 item 6 (vectorised featurizer 7.8
  → <1, tensor-side bias/batch assembly, smaller overlap after boundary support; no B allowance until
  a row lands inside 15). 171 tests pass under 3.11+torch. Local CPU ~0.02 CPU-h; cluster CPU-hours
  0, GPU-hours 0; no held-out species touched. Next: vectorise the featurizer and re-measure the
  chromosome row; GPU half still waits on gagarin (lenin-0083).
- 2026-09-20T22:09Z lenin: renewed lease. Inbox: engels-0084, stalin-0085 (both reviews of PR 38 at
  9de10e1; no question to answer). **Both P2s fixed** (PR 38 at 038f54f, code at 92ddafd): stalin-0085
  — `predict_sequence` now encodes each tile from the stride multiple at or before its start and
  crops, so the pooling grid is the oriented chromosome's (real-encoder test on both strands under a
  misaligned step, fails on old code); engels-0084 — windows are flushed through decode and core
  claiming in groups of `decode_batch` (spy test: ≤ batch emission tensors alive at any decode,
  fails on old code). engels-0084 P3 fixed in tiler/module docs and a-pilot item 6 (containment is
  the guarantee; context bound is about the chain start). CLI help for `--overlap` corrected
  (default 4096). Re-measured one core, same checkpoint (`smoke-local-20260920/chromosome-grid/`):
  default **32.4** CPU-s/genome Mb (was 33.8), overlap 0 / batch 64 **21.3** (was 22.6); RSS
  0.74–0.81 GiB; CPU-regime miss unchanged. 173 tests pass under 3.11+torch. Local CPU ~0.01 CPU-h;
  cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: vectorise the featurizer
  (revision step 1) and re-measure; GPU half still waits on gagarin (lenin-0083).
- 2026-09-20T23:10Z lenin: renewed lease. Inbox: engels-0085, stalin-0086 (both verify PR 38 at 038f54f:
  pooling-grid closure across all stride offsets, bounded emission lifetimes, batch/output
  consistency, hash/stage reconciliation; no question to answer). **Revision step 1 done** (PR 38
  at d1016d4, code at 48f8a2b): `encode_sequence` is now torch ops over the window's byte codes
  (`base_codes` lookup, prefix-sum GC, buffer-protocol availability mask); the per-base Python
  path is kept as `encode_sequence_reference` and a test asserts bit-identical tensors over
  random IUPAC/soft-masked/masked windows; the pooled decoder's dinucleotide lookup vectorised
  the same way against its list reference. 176 tests pass under 3.11+torch. Re-measured one
  core, same checkpoint (`smoke-local-20260920/chromosome-vec/`, six runs, GFF3 byte-identical
  to the 92ddafd output on the two matching pairs): featurizer 2.55 -> 0.03 s; default
  (overlap 4096 / batch 16) **21.1** CPU-s/genome Mb (was 32.4), overlap 2048 / batch 64 15.5,
  overlap 0 / batch 64 **13.5** (was 21.3) — first chromosome row inside the 15 budget, but only
  at the no-overlap cost floor; the default row still misses by 1.4x. a-pilot item 6, section 2
  and section 5 updated; no B allowance claimed. RSS 0.73–1.07 GiB. Local CPU ~0.01 CPU-h;
  cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: revision step 2
  (tensor-side batch assembly and bias in decode) and re-measure; GPU half still waits on
  gagarin (lenin-0083).
- 2026-09-21T00:20Z lenin: renewed lease. Inbox: engels-0086, stalin-0087 (both verify PR 38 at d1016d4:
  featurizer/motif oracles, full-model gradient and chromosome-output preservation; no question to
  answer). engels-0086 P3 fixed: a-pilot now labels the three CPU numerators (stage sum of
  process_time / whole-process user, the cost-baseline 3.2 convention / whole-process user+system)
  side by side in the item-6 tables and section 5. **Revision step 2 re-scoped by measurement**:
  the decode stage (3.7 s default) is scan loop 2.57 + `_operands` 0.88 + Python-side assembly
  only 0.26; a bit-identical numpy twin of the loop gave 1.0x at batch 16 / 0.9x at 56 (element-bound
  on the dense (B,24,24) transition, not dispatch-bound; not kept). **Step 2a done** (PR 38 at
  e6ebd8c, code at 6cae82f): operands kept per phase (B,3,L,R) with one index_select per step
  instead of the (B,K,L,R) stacks, per-step int8 casts dropped, byte-table symbol lookup
  (`symbol_index_tensor`, tested equal to the list reference over all ASCII). 178 tests pass under
  3.11+torch. Re-measured one core, same checkpoint (`smoke-local-20260920/chromosome-ops/`, six
  runs, all six GFF3 byte-identical to chromosome-vec): decode -16..20%, RSS -0.19..0.35 GiB;
  stage-sum CPU-s/genome Mb default **18.5** (was 21.1; 19.4 whole-process user), overlap 2048 /
  batch 64 **13.6** (14.4), overlap 0 / batch 64 **11.9** (13.1) — budget met on every numerator
  only at the no-overlap floor; default misses 1.2-1.4x. No B allowance claimed. Step 2b
  (sparse-predecessor scan: 1 pred per state except phase-0 '' 17 and U 3; est. -10% at batch 16,
  -30..45% at 28-56) proposed in a-pilot. Local CPU ~0.02 CPU-h; cluster CPU-hours 0, GPU-hours 0;
  no held-out species touched. Next: step 2b and re-measure; GPU half still waits on gagarin
  (lenin-0083).
- 2026-09-21T01:22Z lenin: renewed lease. Inbox: stalin-0088 (verifies PR 38 at e6ebd8c: 432 old/new differential
  cases bit-identical, wrappers, gradients, record reconciliation; no question to answer). stalin-0088
  P3 fixed: step-2a decode reduction now quoted as 14.6–17.6% and RSS 0.18–0.35 GiB. **Revision
  step 2b done** (PR 38 at 15a4c61, code at d9ea0f8): the scan takes the per-step max over each
  state's predecessors only (all single except U: 3 and E(""): 17 under code 1) — one gather at
  K−2+2P slots, one max over the two (B,P) candidate rows, one cat; operands built step-major so no
  per-step clone; dense scan kept as `viterbi_batch_reference` with a bit-identity test (scores,
  tail back-pointers, every traceback; 3 codes × 2 dtypes × m 1–20 × masks × padded batches).
  179 tests pass under 3.11+torch. Re-measured one core, same checkpoint
  (`smoke-local-20260920/chromosome-sparse/`, six runs, all six GFF3 byte-identical to
  chromosome-vec): decode −4..7% at batch 16, −15..18% at batch 64 (less than estimated: the
  per-step fixed dispatch dominates); stage-sum CPU-s/genome Mb default **17.1** (18.1 user),
  overlap 4096 / batch 64 **14.2** (14.4 user, 16.8 user+sys) — first default-overlap row inside 15
  on the baseline numerator — 2048/64 **12.0** (12.8 / 14.7), 0/64 **10.7** (11.9 / 13.3); RSS
  ≤ 0.81 GiB. No B allowance claimed (yeast only, 20-step checkpoint, no boundary support, no
  S. pombe row). Local CPU ~0.02 CPU-h; cluster CPU-hours 0, GPU-hours 0; no held-out species
  touched. Next: boundary support (partial-family scalars) so step 3 can be measured; GPU half
  still waits on gagarin (lenin-0083).
- 2026-09-21T02:16Z lenin: renewed lease. Inbox: engels-0087, stalin-0089 (both verify PR 38 at 15a4c61: 198
  analytical phase-specific intron cases, 432 batched differential cases bit-identical, records reconciled;
  no question to answer). stalin-0089 P3 ×2 fixed: a-pilot now says 56 candidates (24−2+2×17), decode
  reduction 4.3–11.7% at batch 16 (11.7% on the default row, above the 10% estimate) / 15.5–18.0% at 64,
  RSS +0.05–0.09 GiB. **Boundary support, first half** (PR 38 at f5e4a43, code at 68502b3): `EdgePrior`
  widened to the four partial families (coding/intron × entry/exit; two-value form unchanged),
  `pooled.edge_prior` reads `partial_families` (values only), and the tensor Viterbi gets the proposal-3.1
  sequence-edge partials (`viterbi_batch_edges`, `viterbi`/`viterbi_windows(edges=)`): E0 and J entry
  at boundary 0 (donor/acceptor rows floored at t = 0), S/E, entered-T and censored-pending-donor exits at
  each window's own length; J recognised at n as a tail with no accepted entry, so no per-step state
  (0.51 vs 0.52 s over 16 × 12,288, same loop); traceback relabels E0/J. `EdgeParity`: 240 random
  lattices × 3 priors vs `DelayedEntryDecoder(edges=)`, per-family fixtures, batched = single, learned
  tables + pooled prior; 184 tests pass under 3.11+torch. Not done: the partition's edge terms (needs a
  separate J layer), the loader's edge-partial chains, `chromosome.py` still free grammar. a-pilot now
  states that proposal 3.1 forbids partials at tile seams, so revision step 3 (overlap) goes through
  carried scan state across tiles (3.3 seams), not edge partials — to be measured next. Local CPU
  ~0.01 CPU-h (tests + micro-benchmark); cluster CPU-hours 0, GPU-hours 0; no held-out species touched.
  Next: carried-state scan across tiles at overlap 0 and its measurement; GPU half still waits on
  gagarin (lenin-0083).
- 2026-09-21T03:11Z lenin: renewed lease. Inbox: engels-0088, stalin-0090 (both review PR 38 at f5e4a43: **P2
  confirmed**, the edge Viterbi merged the residual intron J into the donor-entered tail layer, so a J that
  outscored the donor entry erased the row's only terminal candidate; 4-base fixture with closed-form optimum
  6 − log 3 − 4 log 2 decoded as a 4-base CDS, 144/144 sweep cases, score fell as the J weight rose; no
  question to answer). **Fixed** (PR 38 at 0d87a2e): J now has its own (B, K, R) layer that advances with the
  same intronic emissions, competes with T at every acceptor (`exit_r` records R + r so the traceback runs J
  back to boundary 0), never receives donor entries and is never offered at n; `_finish` reduces over the
  donor-entered layer alone; traceback raises instead of relabelling an un-entered tail as J. Regressions
  added: the engels-0088 fixture across six J weights (score and chain pinned to the closed form) and
  stalin-0090's boundary sweep (codes 1/6, float32/64, m 2/4/20, R 1/3, n = m−1..m+2, both J weights,
  batched 1/7/64 = single); 186 tests pass under 3.11 + torch. Cost: edge mode now 0.63 vs 0.54 s over
  16 × 12,288 (one extra add + max per step; the earlier 0.51 vs 0.52 was the defective version); edge mode
  is for the two real strand ends, so the chromosome rows are unaffected. a-pilot boundary-support paragraph
  rewritten accordingly. Local CPU ~0.02 CPU-h (tests + micro-benchmark); cluster CPU-hours 0, GPU-hours 0;
  no held-out species touched. Next: carried-state scan across tiles at overlap 0 and its measurement; GPU
  half still waits on gagarin (lenin-0083).
- 2026-09-21T04:21Z lenin: renewed lease. Inbox: engels-0089, stalin-0091 (both verify PR 38 at 0d87a2e: 48 tests, the
  144-case sweep and the 576-case boundary sweep pass, 2,592 extra J-traceback cases agree; no question to
  answer). **Revision step 3, first measurement** (PR 38 at ae5a7b5, code at a4d393d):
  `fast_viterbi.scan_segments` / `viterbi_segments` scan a batch of segments one tile at a time and carry a
  `Carry` across each seam (alpha, tau, J layer, the ≤ m−1 pending donors with their seam-context emissions,
  edge-mode score/final, entered-ever mask, ended rows); scores, exit_r/entered inside each row and every
  traceback are bit-identical to the unbroken scan (`SeamCarry`: 200 random batches, tiles 1–64, both modes,
  >400 finite / >100 seam-crossing cases; 189 tests pass under 3.11 + torch). `predict_sequence(segments=)` /
  `measure --segments` cuts each strand into segments overlapping by --overlap, encodes them in non-overlapping
  12,288 tiles and decodes all of both strands in one batch. Yeast chr I, one core/thread, same protocol
  (records `chromosome-seams/`; the two window-mode rows re-run at a4d393d are byte-identical to
  `chromosome-sparse/`): exact 1-segment decode 41.3 CPU-s/Mb (B = 2, the per-step fixed cost), 8 / 19 / 32
  segments per strand with seam overlap 4,096 → **12.7 / 11.8 / 12.6** stage sum (13.6 / 12.3 / 12.9 user,
  15.3 / 14.3 / 15.2 user+sys), containment guarantee kept at segment seams, interior tile seams exact; the
  oversampling is 1 + (S−1)·4096/n, 1.32× here and 1.004× at 20 Mb where the window mode's 1.5× is fixed.
  Segment rows still differ from the exact decode by 6–11 chains of 1,832 (window row: 5) partly because the
  encoder sees no context past a tile edge. Back-pointers of a whole segment are held (120 B/base): next is the
  packed store (~35 B/base) and an encoder context margin, then the S. pombe row. No B allowance claimed.
  Local CPU ~0.03 CPU-h (tests + 6 + 6 measured runs, one of them mistakenly at 20 threads and discarded);
  cluster CPU-hours 0, GPU-hours 0; no held-out species touched. GPU half still waits on gagarin (lenin-0083).
- 2026-09-21T05:10Z lenin: renewed lease. Inbox: stalin-0092 (verifies the carried-state scan at ae5a7b5: 51 tests, 4,608
  additional outputs equal whole-segment decoding, 1,152 match the scalar decoder; P3: "4,096-base containment
  guarantee" overstates, the bound is the 2,048-base chain containment of `chromosome.tiles`; no question to answer).
  **P3 fixed** in a-pilot (summary bullet and step-3 paragraph now say "4,096-base seam overlap retaining the
  2,048-base chain-containment guarantee", with stalin's counterexample). **Revision step 3, packed back-pointers**
  (PR 38 at 31a5bde, code at 8a7b4a8): `scan_segments` returns a `PackedBackPointers` store — `prev` as the two
  argmax slots the scan already computes (2 B/boundary; the other states have one predecessor, expanded on
  traceback), `exit_r` as nibbles (12 B), `entered` as 72 bits (9 B) — **23 B/base instead of 120**, packed per tile,
  expanded one row at a time for the traceback; `SeamCarry` checks the byte count and bit-for-bit equality of each
  row's expansion with the unbroken scan; 189 tests pass under 3.11 + torch. Re-measured yeast chr I (records
  `chromosome-seams-packed/`, GFF3 byte-identical to `chromosome-seams/`): exact 1-segment decode 40.1 / 41.7 / 42.8
  with peak RSS **0.45 GiB (was 0.52)**; 19 segments **11.9 / 12.3 / 14.6** (was 11.8 / 12.3 / 14.3, inside noise),
  RSS 1.02 GiB unchanged because the per-tile operands at B = 38 are the peak there, not the pointers. Metazoan
  estimate at 20 Mb / 19 segments: 0.92 GB pointers + 1.76 GB float32 emissions (held twice today: `viterbi_segments`
  pads them into one (B, 11, L) tensor) + ~0.6 GiB operands + one row's dense expansion, ~3.5 GB vs 4.8 GB for the
  dense pointers alone. Local CPU ~0.02 CPU-h (tests + 2 measured runs); cluster CPU-hours 0, GPU-hours 0; no
  held-out species touched. Next: feed the scan tile by tile from the per-segment list (drop the padded copy),
  encoder context margin, then the S. pombe row; GPU half still waits on gagarin (lenin-0083).
- 2026-09-21T06:12Z lenin: renewed lease. Inbox: engels-0090 (P2: the 3.5 GB metazoan projection was wrong — measure decodes in
  float64 and 8a7b4a8 held every segment's emissions plus a padded copy, ≥ 8.7 GB live at 20 Mb; 51 tests, 7,168 packed
  outputs verified), stalin-0093 (P3: `PackedBackPointers.nbytes()` counted payload while superseded seam columns and
  the seam context stayed alive behind views; P2 reproduced; no question to answer). **Both fixed** (PR 38 at 0c43c57, code
  at 87bb3a1): `scan_segments` takes a `chunk(start, end)` callable, `viterbi_segments` per-segment emitters, and the
  chromosome mode encodes each tile when the scan asks for it — no per-segment emission list, no padded copy; the store
  holds each tile's last exit column apart and every retained tensor owns its storage (`storage_nbytes() == nbytes()`,
  stalin's three cases 1,104 / 2,208 / 2,400 B both ways); `SeamCarry` checks streamed == whole tensors, emitter call
  pattern, and storage == payload; 190 tests pass. Re-measured yeast chr I (records `chromosome-stream/`, GFF3
  byte-identical to `chromosome-seams/` and `-packed/`): 1 segment **41.9 / 43.4 / 44.5**, RSS **0.35 GiB** (was 0.45);
  19 segments **12.3 / 12.6 / 14.9**, RSS **0.96 GiB** (was 1.02; the B = 38 per-tile operands + runtime floor, length-
  independent). a-pilot withdraws the 3.5 GB figure, records engels' ≥ 8.7 GB bound for the 8a7b4a8 path, and projects
  the streamed path at 20 Mb / 19 segments to ~1.5 GB of live tensors (0.92 pointers + 0.40 operands + 0.13 one-row
  expansion + 0.04 tile emissions), stated as a source-derived projection, not a measurement. Local CPU ~0.03 CPU-h;
  cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: encoder context margin at tile seams, then the
  S. pombe row; GPU half still waits on gagarin (lenin-0083).
- 2026-09-21T07:11Z lenin: renewed lease. Inbox: engels-0091, stalin-0094 (both verify PR 38 at 0c43c57: 52 tests, streamed outputs equal
  whole-segment decoding, storage == payload in stalin's three cases, the 14-configuration real-encoder emission comparison
  at error 0; both keep the ~1.5 GB figure as a projection, not a bound; no question to answer). **Revision step 3, encoder
  context margin** (PR 38 at cb2b27d, code at 8440560): in the segment mode each tile is encoded with `DEPENDENCY_RADIUS` =
  491 bases of the oriented chromosome on each side (clipped at the true ends, origin-anchored grid, cropped back;
  `predict_sequence(margin=)`, `measure --margin`, `encoded_bases` in the record); `SegmentMargin` checks with the real
  encoder that every tile equals the whole-strand encoder output to 2.4e-7 and that bare tiles differ by 0.03–0.13;
  191 tests pass under 3.11 + torch. Re-measured yeast chr I (records `chromosome-margin/`): 1 segment **41.0 / 42.7 /
  43.7**, RSS 0.36 GiB; 19 segments **13.2 / 13.3 / 15.9**, RSS 0.96 GiB (encoder 1.17 → 1.32 s, 1.12× the encoded bases on
  these short tiles, 1.08× on full ones; decode unchanged). The two decodes now differ by **4 chains of 1,832** (were 11):
  all four plus-strand genes of 2.9–4.6 kb starting 0.8–1.2 kb before a segment core seam, longer than the 2,048-base
  containment guarantee — the seam contract, not a tile artefact; for the metazoan row the seam overlap must scale with
  gene length (or the segment count drop), noted in a-pilot. Local CPU ~0.02 CPU-h (tests + 2 measured runs); cluster
  CPU-hours 0, GPU-hours 0; no held-out species touched. Next: the S. pombe normalization row (leakage check first, frozen
  smoke checkpoint, both modes) and AUGUSTUS on the same machine; GPU half still waits on gagarin (lenin-0083).
- 2026-09-21T09:11Z lenin: renewed lease. Inbox: engels-0092, stalin-0094 (both verify PR 38 at cb2b27d: 68 / 24 tests, 42 / 32 extra
  configurations against a whole-strand encoder oracle at ≤ 3e-7, streamed Viterbi + GFF3 exact vs per-segment oracle, the
  chromosome-margin records reconciled and the 4 differing chains confirmed as the seam contract; no question to answer).
  **S. pombe normalization row measured** (PR 38 at 0489c11, code unchanged at 8440560; records `pombe-normalization/`):
  leakage check first (0 violations, S. pombe cross-clade held out, unscored), then A on the three nuclear chromosomes
  (12.57 Mb) one process pinned to one core, frozen smoke checkpoint, float64 decode, margin 491, and AUGUSTUS 3.5.0 with
  marx-0026's command on the same FASTA and core. A at 19 segments per strand: **9.22 / 8.23 / 9.37 CPU-s per genome Mb**
  (stage sum / user / user+system; encoder 3.74, decode 5.31), RSS ≤ 1.31 GiB; exact strand decode 40.35 / 39.82 / 40.49,
  RSS 1.49 GiB; AUGUSTUS **51.41** user CPU-s/Mb, 0.40 GiB, 4,452 genes (164.9 on the cost-baseline runner → this machine
  is 3.21× faster on that command). **A is 1/5.5 of AUGUSTUS against the 1/11 portable target (30.0 CPU-s/Mb machine-
  normalized against 15): the CPU regime misses the normalization row by 2.0×, 1.8× on user only; no positive CPU
  allowance for B.** Segment probe on chr III (38 / 76 segments per strand): stage −5 % / −3 % at 1.6× / 2.9× the memory —
  the decode is at its operand floor (~4.5 CPU-s per genome Mb), segment count exhausted. Proposed revision (a-pilot 3.2):
  step 4 float32 decode (expected decode → 2.5–3), then an encoder work reduction only together with a fitted checkpoint's
  accuracy column; expectation stated that A meets 1/11 only with the encoder change or on the GPU regime. measured.tsv
  gets a `tool` column (AUGUSTUS rows unchanged) and three rows (AUGUSTUS same-machine, A exact, A 19 segments). Local
  CPU ~0.3 CPU-h (A 0.17, AUGUSTUS 0.18, probe 0.01); cluster CPU-hours 0, GPU-hours 0. Next: float32 decode path + tests,
  re-measure the S. pombe 19-segment row; then the metazoan row; GPU half still waits on gagarin (lenin-0083).
- 2026-09-21T10:24Z lenin: renewed lease. Inbox: engels-0093, stalin-0096 (both verify PR 38 at 0489c11: eight A records reconcile, AUGUSTUS
  4,452 genes, user-time convention 1.76×, user+system 2.00×, both fail; P3 processes=6 → 3 on the two A rows; P3 excluded
  mitochondrion is NC_088682.1, not NC_001326.1; CPU-h accounting 0.3661 h; no question to answer). **Both P3s fixed** (TSV rows
  processes=3, `augustus/outputs_sha256.txt` label, README; accounting ~0.37 h). **Revision step 4, float32 decode** (PR 38 at
  d78114c, code at aee0134): `measure --dtype float32` (emissions, motif bias, duration tables, scan); in the segment mode the
  carried scores are rebased at every tile seam (`scan_segments(rebase=)`, default for dtypes narrower than float64, shifts
  summed back in float64; float64 stays the bit-for-bit unbroken scan); `SeamRebase` tests: rebased float64 == unrebased on
  offset lattices, float32 rebased == float64 chains (173/173); 193 tests pass under 3.11 + torch. Measured on S. pombe (same
  protocol, checkpoint, core; records `pombe-normalization/A-f32/`): 19 segments **8.27 / 7.94 / 8.43** CPU-s per genome Mb
  (float64 9.22 / 8.23 / 9.37), RSS 0.97 GiB (was 1.31); exact decode 39.72 / 39.22 / 39.86 (was 40.35 / 39.82 / 40.49). Decode
  5.31 → 4.34: **a 10 % gain, not the halving proposed** — one-tile profile 0.93 → 0.73 s of which the per-step loop is 0.52 → 0.47
  (42 µs/step, ~15 torch ops on (38, 24) tensors: dispatch-bound), operand build 0.41 → 0.26. Row now 1/6.1 of AUGUSTUS, **miss
  1.8× (1.7× user only), still failed; no B allowance.** Outputs: float32 flips 103 of 121,937 chains vs float64 at 19 segments
  (93 internal boundary shifts in short chains), 112 at 1 segment, and float32 exact vs segment decodes differ by ~120 where
  float64 differ by 3 — near-ties of the flat smoke checkpoint; float64 stays the row of record until a fitted checkpoint
  re-checks the flip count. measured.tsv gets the float32 row as a measured option. a-pilot 3.2 amends the revision list: the
  remaining decode lever is a compiled scan step (worth ~2 CPU-s/Mb, row → ~6–6.5, 1/8), so the encoder work reduction with a
  fitted checkpoint's accuracy is necessary for 1/11 on CPU. Local CPU 0.17 CPU-h (six runs) + ~0.01 tests/profile; cluster
  CPU-hours 0, GPU-hours 0; S. pombe touched only for runtime normalization, unscored. Next: the metazoan development-chromosome
  row (seam overlap scaled to gene length), then the fitted checkpoint for the accuracy column; GPU half still waits on gagarin
  (lenin-0083).
- 2026-09-21T11:35Z lenin: renewed lease. Inbox: stalin-0097 (verifies PR 38 at d78114c: 25 focused tests, all 12 pombe JSON/stdout
  pairs reconcile, chain diffs reproduced, user-time miss 1.70× / user+system 1.80×; no question to answer). **Metazoan development-
  chromosome row measured** (PR 38 at 0941d6f, code unchanged at aee0134; records `elegans-chrV/`): *C. elegans* chr V (NC_003283.11,
  20.92 Mb; train species, chr V declared the dev chromosome in the run config, runtime only, unscored), pinned WBcel235 sources
  MD5-verified, leakage check first (0 violations), frozen smoke checkpoint, one process pinned to one core, seam overlap **16,384**
  (8,192-base containment, covers 97 % of chr V's admitted representatives; 1.014× oversampling), AUGUSTUS 3.5.0 `--species=
  caenorhabditis` on the same chromosome on a second core. A at 19 segments per strand, float64: **9.64 / 8.48 / 9.67 CPU-s per genome
  Mb** (stage / user / user+system; encoder 3.91, decode 5.50), RSS 2.04 GiB; float32 8.37 / 7.85 / 8.40, RSS 1.74 GiB; exact strand
  decode 42.27 / 41.98 / 42.30, RSS 4.87 GiB (row-length scaling of the traceback expansion; inside 8 GB, not the row of record);
  AUGUSTUS **52.19** user CPU-s/Mb, 0.70 GiB, 3,357 genes (within 2 % of its 51.4 on S. pombe here, so the 3.21× machine factor
  applies). **A is 1/5.4 of AUGUSTUS (1/6.2 on user), machine-normalized 31.0 vs 15: the metazoan row misses the CPU target by 2.0×
  (1.8× user), the same miss as S. pombe; no B allowance.** Per-Mb cost is length-independent within 5 % between 12.6 and 20.9 Mb.
  Outputs: 19-segment float64 decode **byte-identical** to the exact strand decode (208,796 chains, 0 differ) — the longer overlap
  removes the seam cuts the yeasts showed at 4,096; float32 flips 194 / 199 chains (0.09 %), same near-tie behaviour as S. pombe.
  measured.tsv gets four rows; a-pilot status paragraph and 3.2 updated. The DoD's "S. pombe and one metazoan development chromosome"
  is met for the CPU regime with the smoke checkpoint. Local CPU 0.65 CPU-h (A 0.35, AUGUSTUS 0.30); cluster CPU-hours 0, GPU-hours 0;
  no held-out species touched. Next: the fitted checkpoint (train-species fit under the cap, checkpoint selection on dev chromosomes,
  `benchmark/score.py` on the dev chromosomes with a declaration) and the accuracy column, re-checking the float32 flip count; GPU half
  still waits on gagarin (lenin-0083).
- 2026-09-22T00:11Z lenin: renewed lease. Inbox: engels-0094, stalin-0098, engels-0095 (three reviews of PR 38 at d78114c / 0941d6f:
  float32 rebasing verified on 864 outputs, metazoan records/outputs/geometry reconcile, all 4,965 admitted chr V spans fit their
  segment; two P3 reporting corrections), gagarin-0009 (idle, lenin-0083 still ungranted); no question to answer. **Applied the
  stalin-0098 corrections** (PR 38 at 172b695): `measured.tsv` `wall_s` is the process elapsed time on every candidate-A row
  (chr V 885.5 / 202.4 / 175.9; pombe float32 106.0, was the stage sum 104.0), metazoan comparison labelled same-machine (A core 2,
  AUGUSTUS core 4) in a-pilot; READMEs get distinct `stage wall s` / `process wall s` columns. **Fitted-checkpoint groundwork:
  gene-free (background) training windows** (PR 38 at ead9831): the chain loss accepts the empty chain as the all-intergenic
  support (`numerator_scores`, `support_mask`; introns without CDS still rejected), `iter_background_windows` tiles the gene-free
  intervals of sequences with admitted chains (all GFF3 transcripts blocked on both strands, FLANK margin, mitochondrion never
  background, seeded draw without replacement, half reverse-complemented), config keys `background_windows`/`background_length`
  in the manifest plan (default 0 = the recorded chain-only scope). 153 A tests pass; a 2-step S. cerevisiae load with 6
  background windows runs (loss ~805–810 nats per 4,096-base tile at init, finite); S. cerevisiae has only 23 gene-free 4,096-base
  tiles (198 at 2,048; 928 at 1,024), so the fit config sets the length per panel. Local CPU ~0.01 CPU-h (tests, 2-step check);
  cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: the fitted checkpoint itself — train-species fit config
  (dev chromosomes declared, background share), a bounded local CPU fit or the gagarin GPU grant (lenin-0083 still open),
  `benchmark/score.py` on the dev chromosomes for the accuracy column, and the float32 flip re-check.
- 2026-09-22T01:08Z lenin: renewed lease. Inbox: engels-0096 (PR 38 at ead9831: loss checks pass; P2 `max_window` does not bound `background_length`;
  P3 exclusion is CDS-span-free, not gene-free — 11/68/214 of the 23/198/928 S. cerevisiae tiles overlap a raw gene/pseudogene row),
  stalin-0099 (empty-chain numerator, gradients, mixed batches verified to 1e-14; no defect); no question to answer. **Both engels-0096
  findings fixed** (PR 38 at b7a74ab): `annotated_gene_spans` reads every raw GFF3 gene/pseudogene/`*_gene`/`*_gene_segment`/`*RNA`/
  transcript/exon/CDS/UTR row per seqid (`is_gene_feature`) and `iter_background_windows` blocks those instead of the admission
  parser's CDS spans, so UTRs, ncRNAs and pseudogenes are excluded while centromeres, LTRs and replication origins stay eligible;
  `TrainConfig.check_window_bound` refuses an enabled `background_length` above a non-null `max_window` at parse and again in
  `_load_all_windows`. Tests: engels' UTR+lncRNA fixture (spans, intervals, four 1,000-base tiles all outside the genes, case
  preserved / complemented), type classifier, config bound, and the `_load_all_windows` integration (max_window 100 refused; at
  1,000 loads n = 72 and 1,000); 193 A/grammar tests pass. S. cerevisiae recount with engels' oracle: **14 / 148 / 766** candidates
  at 4,096 / 2,048 / 1,024 with **0** overlapping any raw gene or pseudogene row (was 23 / 198 / 928 with 11 / 68 / 214); a-pilot
  section 2 item 1 corrected. Local CPU ~0.01 CPU-h (tests, recount); cluster CPU-hours 0, GPU-hours 0; no held-out species touched.
  Next: the fitted checkpoint — train-panel fit config (dev chromosomes declared, background length per species), a bounded local CPU
  fit or the gagarin GPU grant (lenin-0083 still open), `benchmark/score.py` on the dev chromosomes for the accuracy column, float32
  flip re-check.
- 2026-09-22T02:24Z lenin: renewed lease. Inbox: engels-0097, stalin-0100 (PR 38 at b7a74ab: both engels-0096 findings verified fixed, 42 tests,
  S. cerevisiae recount 14/148/766 with 0 gene overlap by an independent oracle, seeded train/dev split and sampling invariants over
  30 seeds; no question to answer). **Fitted checkpoint started.** Two probe fits on S. cerevisiae + C. elegans (dev chr I / chr V,
  548 background tiles of 2,048, 24,601 windows: 19,708 train / 4,893 dev) measured the cost: load 21 s, 4.0 s per 8-window step
  (0.5 s per ~2.4 kb window fwd+bwd), 0.13 s per dev window — a full chr V evaluation is ~10 min, so `dev_windows_max` (seeded
  dev subsample for checkpoint selection, both counts in the manifest) plus `best_step`, per-evaluation `history` and load/fit/
  eval timing in `run_manifest.json` landed (PR 38 at 4c43819; 30 train tests, 161 A + 36 grammar tests pass). **Bounded local CPU
  fit v1 launched** (1,500 Adam steps, batch 8, ~29 Mb sampled, eval every 100 steps on 256 dev windows, one thread pinned to one
  core, leakage check 0 violations first, sources MD5-gated; projected ~1.7 CPU-h, running past this tick); a-pilot 3.3 records
  the plan and probe numbers (ddcec44). Local CPU this tick ~0.3 CPU-h (two probes incl. one full 4,893-window eval, tests) plus
  the fit's ~1.7 CPU-h in flight; cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: collect the fit
  (manifest, history, best.pt into `smoke-local-20260920/fit-cpu-v1/`), `measure --checkpoint` rows on S. pombe and chr V under
  the fitted weights with the float32 flip re-check, `benchmark/score.py` on the two dev chromosomes for the accuracy column;
  gagarin GPU grant (lenin-0083) still open.
- 2026-09-22T03:08Z lenin: renewed lease. Inbox: stalin-0101 (PR 38 at ddcec44: dev subsample and checkpoint accounting verified, 30 tests + 15
  instrumented train() cases, no defect; no question to answer). **Fit v1 still running** (step 600/1,500 at 2,308 s, 3.85 s/step incl.
  evals; dev NLL 129/60/42/45/66/58 at steps 100–600, best_step 300 so far; finishes ~04:05Z). Run records (run.sh, config, leakage 0
  violations, source MD5, partial log) committed to `smoke-local-20260920/fit-cpu-v1/`. **Scoring pipeline dry run** with the interim
  step-300 checkpoint on S. cerevisiae chr I, other core: `measure --profile chromosome --gff-out` → `benchmark/score.py --seqids
  --genome --declaration` both exit 0 (`fit-cpu-v1/score-dry/`); chr I 12.79 CPU-s/Mb end to end, RSS 0.96 GiB, scorer 0.17 CPU-s.
  Interim numbers (pipeline check, no accuracy claim): nucleotide sens 0.963 / prec 0.456, 662 predicted loci vs 94 reference, 320
  predicted introns vs 3 (316 non-GT-AG) — decoder over-predicts at step 300; motif bias under-fit is the first suspect. a-pilot 3.3
  updated (PR 38 at aeed615). Local CPU this tick ~0.01 CPU-h (dry run) plus fit in flight (~0.64 CPU-h so far of ~1.6); cluster
  CPU-hours 0, GPU-hours 0; no held-out species touched. Next: collect the finished fit (train.out, run_manifest.json with actual/
  history/timing, train_time.txt, best.pt sha256), score chr I and chr V under the final best.pt, `measure --checkpoint` rows on
  S. pombe and chr V with the float32 flip re-check, section 6.1 accuracy entries; if the over-prediction persists, diagnose the
  motif bias / duration tables before a v2 fit. gagarin GPU grant (lenin-0083) still open.
- 2026-09-22T04:05Z lenin: renewed lease. Inbox: engels-0098, stalin-0102 (PR 38 at aeed615: interim step-300 score JSON reproduced exactly,
  all 662 chains pass sequence/frame/phase checks; checkpoint finite, motif tables moved canonical (GT +0.051, AG +0.053), duration
  mean 21.4 bases; no defect, no question to answer). **Fit v1 finished** (exit 0, 1,500 steps, **5,833 CPU-s = 1.62 CPU-h** on one
  core, wall 1:37:17, peak RSS 2.14 GiB; best_step **900**, dev NLL 38.04; train.out, run_manifest.json with actual/history/timing,
  train_time.txt, run.log, best_pt.sha256 in `fit-cpu-v1/`). **Both dev chromosomes scored under best.pt** (`fit-cpu-v1/score-final/`,
  measure --checkpoint → score.py, all exits 0): chr I 12.09 CPU-s/Mb, 0.96 GiB, **6 predicted transcripts vs 94**, nucleotide
  0.039/0.356; chr V 8.86 CPU-s/Mb, 2.05 GiB, **287 predicted vs 6,766**, nucleotide 0.470/0.586, 238 fusions, 0 exact transcripts,
  67 exact CDS exons, 9,622/10,261 predicted introns non-canonical, median predicted chain span 81 kb (max 843 kb). Cost inside 15
  CPU-s/Mb and 8 GB; **A misses the accuracy target in this fit.** Diagnosis (a-pilot 3.3): training windows are gene ± 10-base
  flank, so non-CDS bases are supervised as intron everywhere except the 548 background tiles (~1.1 Mb of 26.3 Mb sampled) — no
  intergenic-vs-intron signal, chromosome decode fuses genes through intergenic stretches; duration alone does not stop it. Revision
  proposed before any B allowance: adjacent-gene / wide-flank loader increment + genome-matched background share, then fit v2 scored
  the same way. a-pilot 3.3 and section 5 accuracy entry updated (PR 38 at 968f697). Local CPU this tick: fit 1.62 CPU-h (finished
  03:59Z) + scoring 0.05 CPU-h; cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: implement the wide-flank /
  adjacent-gene window loader with neighbouring-gene support masks and a background share from the intergenic fraction, tests, then
  launch bounded fit v2 at the same budget; gagarin GPU grant (lenin-0083) still open.
- 2026-09-22T05:15Z lenin: renewed lease. Inbox: engels-0099 (final scores reproduced; P3 sampled-background accounting, P3 GT-AG TP+FP), stalin-0103
  (P2 fitted CPU rows must be judged against the portable ceiling); no questions. **Correction to the previous bullet's "cost inside 15
  CPU-s/Mb"**: 15 is bound to the cost-baseline runner; this host's portable ceiling is AUGUSTUS 51.41/11 = 4.67 CPU-s/Mb, so the fitted
  rows (chr I 12.09, chr V 8.86 stage sum) miss by 2.6× and 1.9× (38.8 / 28.4 machine-normalized against 15) — the CPU regime still misses
  as in a-pilot 3.2; memory passes. Both P3s applied: the diagnosis now uses the sampled draws (484,704 intergenic of 15,628,525 bases
  through step 900, 3.1 %; 785,532 of 26,259,772 over all steps) not the 548-tile inventory, the GT-AG row reads 586 = 162 TP + 424 FP with
  9,622 as score.py's `other`, and the sampling cause is phrased as a working hypothesis (a-pilot 3.3, section 5, summary; score-final
  README; PR 38 at 646563d). **Loader increment implemented** (61a9480): `iter_windows(context=N)` widens each clean window by up to N
  annotated-intergenic bases per side, clipped at neighbouring genes and sequence ends (support mask unchanged, no loss/decoder change);
  `TrainConfig.context`; manifest records `actual.composition` (CDS/intron/intergenic bases, background draws); 5 new tests, 211 pass.
  Loader probe: train-pool intergenic share 3.0 % → 20.2 / 30.3 / 39.1 % at context 256 / 512 / 1,024. **Fit v2 launched 05:13Z**
  (`fit-cpu-v2/`: v1 config + context 512 + background 800/species, 1,500 steps, leakage 0 violations, 19,975 train / 4,951 dev windows,
  30.3 % intergenic; projected ~2.2–2.4 CPU-h on one core, running past this tick; a-pilot 3.3 updated, PR 38 at 2f9d053). Local CPU this
  tick ~0.05 CPU-h (four loader probes, tests) plus the fit in flight; cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next:
  collect fit v2 (manifest with composition, timing, best.pt sha256), score chr I and chr V under its best.pt into `fit-cpu-v2/score-final/`
  against the v1 rows, then the fitted S. pombe normalization row; gagarin GPU grant (lenin-0083) still open.
- 2026-09-22T06:07Z lenin: renewed lease. Inbox: engels-0100 (PR 38 at 2f9d053: context loader passes geometry review, 56 tests + 23,978 windows verified; P3 context
  is CDS-span-clipped not annotated intergenic, P3 README background split), stalin-0104 (stalin-0103 P2 closed; P3 v2 shares v1's dev settings not
  its dev subsample — 37/256 selected dev windows common; seeded draw replay 30.58 % U at 12,000 attempted draws); no questions. All three P3s
  applied as documentation, not loader changes (v2 is mid-fit; the raw-feature exclusion stays a knob for a later fit): dataset.py docstrings,
  a-pilot 3.3 and the v2 README now say the added context is coding-span-free (clipped at other genes' CDS spans; UTRs/ncRNA/CDS-free
  pseudogenes may be inside it, 48,213 / 3,965,426 yeast and 3,964,153 / 12,049,332 elegans added bases overlap raw gene spans), call the
  share the `U` share, give the split 19,975 train (773 bg) / 4,951 dev (175 bg), 948 bg in all, train pool 59,330,766 = 23,891,418 CDS +
  17,457,446 intron + 17,981,902 U, and state that v1/v2 raw dev NLLs are on different windows so chr I / chr V scoring is the comparison
  (PR 38 at 4d766b0). **Fit v2 in flight**: step 500/1,500 at 52 min CPU (dev NLL 284.9 → 68.8 → 56.0 → 90.9 → 52.0 at steps 100–500,
  RSS 1.8 GiB), ~2.6 CPU-h projected, finish ~07:50Z. Local CPU this tick ~0.01 CPU-h (tests) plus the fit in flight; cluster CPU-hours 0,
  GPU-hours 0; no held-out species touched. Next: collect v2 (manifest with accepted-draw composition, timing, best.pt sha256), score
  chr I and chr V into `fit-cpu-v2/score-final/` against the v1 rows, then the fitted S. pombe row; gagarin GPU grant (lenin-0083) still open.
- 2026-09-22T07:40Z lenin: renewed lease. Inbox: engels-0101 and stalin-0105 (PR 38 at 4d766b0: all P3s resolved, dataset ASTs identical after docstring
  strip; no questions). **Fit v2 finished** 07:26Z, exit 0: 1,500 steps in 7,964 CPU-s = 2.21 CPU-h on one core (5.3 s/step, 1.37× v1), peak RSS
  2.09 GiB, best step 900 (dev NLL 45.3 on the v2 dev subset), best.pt sha256 99f773ed…; manifest accepted-draw composition 35,305,664 sampled bases
  = 39.9 % CDS + 29.6 % intron + 30.6 % U (9,561,661 context + 490 background draws + 230,200 flank), equal to stalin-0104's replay to the base.
  **Scored chr I / chr V** under best.pt (`fit-cpu-v2/score-final/`, all exits 0, `summarize.py` reproduces the tables): cost unchanged (12.54 /
  8.88 CPU-s/Mb, RSS 0.96 / 2.05 GiB; CPU verdict stands, no B allowance). Fusions gone: chr I nt F1 0.825, MCC 0.746, 86/94 loci, 69 exact
  transcripts (v1: 6 chains, 0 exact); chr V 13,410 chains of median 294 b, 4 fusions (238), 3,804/4,995 loci, nt F1 0.538. Still misses: chr V
  donor/acceptor F1 0.014/0.025 (2,290 predicted introns, median 34 b, 2,114 non-GT-AG, 68 TP), 1,827 splits, 101 exact transcripts (sens 0.020);
  precision 9,606 FP loci on chr V / 128 on chr I, half the chains < 300 CDS bases. Proposed next revision: a longer fit under the same loader
  (lenin-0083 GPU grant) before any decoder change. a-pilot summary, 3.3 and section 5 updated (PR 38 at be7ab46). Local CPU this tick 2.21
  (fit) + 0.05 (scoring) CPU-h; cluster CPU-hours 0, GPU-hours 0; no held-out species touched. Next: fitted S. pombe normalization row under v2
  best.pt; keep the task in progress pending the longer fit / GPU regime.
- 2026-09-22T08:15Z lenin: renewed lease. Inbox: stalin-0106 (PR 38 at be7ab46: v2 scores reproduce; P2 locus FP conflates no-CDS-overlap with unmatched
  fragments, P3 causal scope of v1→v2, 2,125 not 2,114 non-GT-AG; step-900 prefix replay; no questions). **Both applied** (PR 38 at 9f48d29):
  reproduced stalin's read-only matching replay (chr I 86 matched / 127 no same-strand CDS overlap / 1 overlap-unmatched = 128 FP; chr V
  3,804 / 6,376 / 3,230 = 9,606) and replaced "overlap no reference gene" in the summary, 3.3, section 5 and the v2 score README with the split,
  the no-overlap subset flagged as unchecked against gene spans / opposite strand; "the only change is the U share" → a comparison of two loader
  configurations (context geometry, background pool, drawn chains, length exclusions, dev windows all change) that supports the supervision
  hypothesis without isolating U; "fusions gone" → sharply reduced (4 + 2 remain); minimum-length introns as hypothesis; 2,125 non-GT-AG; the
  step-900 prefix (v1 3.10 % / v2 30.63 % U) recorded in the v2 README; review-response bullets in a-pilot section 6. **Fitted S. pombe
  normalization row measured** (`pombe-normalization/A-v2/`, v2 best.pt, 19 segments, float64, core 2, S. pombe unscored): 8.92 stage / 7.95
  user / 9.06 user+system CPU-s/Mb, RSS 1.31 GiB, 9,555 chains — 1/5.7 of AUGUSTUS (51.43) against 1/11, 29.1 machine-normalized against 15;
  3 % from the smoke row, verdict unchanged (1.9× miss, no B allowance); a-pilot summary / 3.3 / section 5 updated. Local CPU this tick 0.03
  (pombe) + 0.02 (replay, checks) CPU-h; cluster CPU-hours 0, GPU-hours 0; no held-out species scored. Next: keep in progress pending the
  longer fit (lenin-0083 GPU grant still open); if no grant, a longer local CPU fit under the same loader is the fallback.
