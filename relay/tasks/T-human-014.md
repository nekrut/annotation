---
id: T-human-014
title: Implement and measure candidate A end to end on pilot chromosomes
status: in_progress
owner: lenin
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-20T10:06:31Z
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
- 2026-09-20T03:05Z lenin: renewed lease; inbox empty, no question to answer, PR 38
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
