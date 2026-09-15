---
id: T-human-013
title: Validate labels and grammar: admission audit and reference decoder for A
status: review
owner: marx
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-15T21:07:40Z
depends_on: [T-human-011]
touches: [model/, tests/]
pr: https://github.com/nekrut/annotation/pull/31
---

## Goal

Phase 4 milestone 1 of `docs/design/proposal.md` (section 8). Before any
fitting, make the label contract and the grammar real code.

- Run the full train-panel metadata and FASTA admission audit specified in
  proposal section 3.6 on the ten `train` species of `benchmark/panel.tsv`:
  representative selection, conflict masks, partial ends, phase and
  translation exceptions, ambiguous bases. Output one admitted-label
  manifest per species with counts and checksums.
- Implement candidate A's decoder as specified in sections 3.1 to 3.3: the
  coding-prefix states, the duration recurrence with exact delayed entry,
  checkpoints and traceback replay. Implement the slow reference recurrence
  first and the delayed-entry version second, and check them against each
  other.
- Pass every hand-checkable acceptance case in section 3.4, including
  partial edges, split codons, ambiguity and checkpoint replay.

Laptop only. No training, no cluster. Python 3.11, PyTorch for tensors,
tests under `tests/`. Keep the manifests under 5 MB each; commit fetch
scripts, not sequence data.

## Definition of done

- `model/grammar/` (or equivalent) with the reference and delayed-entry
  decoders, and `tests/` that run the section 3.4 cases in under a minute.
- `docs/design/admission-audit.md` with the per-species admitted-label
  counts, the fraction of representatives excluded and why, and the
  checksums of the manifests that later tasks train on.
- A finite legal numerator on every admitted sequence; any sequence where
  it is not is listed with the reason.
- Pull request from `work/T-human-013-<name>`, at least one `review` from
  another agent, coordinator merge.

## Log

- 2026-09-15 human: created.
- 2026-09-15 marx: claimed.
- 2026-09-15 marx: -> in_progress.
- 2026-09-15 marx (16:1xZ): claimed and started. Branch
  `work/T-human-013-marx` at 7e3e47a, draft PR #31
  (https://github.com/nekrut/annotation/pull/31). Added `pyproject.toml`
  (Python 3.11, standard library only for the grammar; torch and pytest
  optional) and `model/grammar/`: the state inventory of proposal 3.1
  (U, S(q), E(q), I(c,k), T(c,r)), genetic-code tables 1 and 6, the 11
  score channels, and the expanded reference recurrence of 3.2 as the
  specification oracle (forward and Viterbi over one transition generator,
  m-1 explicit mandatory intron states plus R tails, ambiguous bases
  branched under a uniform prior, traceback to CDS/intron segments with
  GFF3 phase and a sequence-uncertain flag). `tests/` pass the first five
  section 3.4 cases (phase 1 and 2 prefixes across an intron, split
  initiator, table 6 versus table 1 on TAA, m-1/m/m+1 lengths against
  pi(1-q) and pi q(1-q)) and brute-force path enumeration agrees with the
  partition and the Viterbi score on random small lattices; 8 tests, 0.02 s.
  Finding while testing: with zero intron emissions an intron pinned only
  by the CDS can legally shift one base and end on TGA instead of TAA, so
  the acceptance tests pin donor and acceptor scores. Next: edge partials
  with the partial-end prior, then the delayed-entry recurrence checked
  against this oracle exhaustively on tiny lattices, checkpoint and seam
  replay, reverse-complement agreement, then the section 3.6 admission
  audit on the ten train species.
- 2026-09-15 marx: -> in_progress.
- 2026-09-15 marx (17:0xZ): pushed 3450d07 on `work/T-human-013-marx` (PR #31).
  Fixed both interim findings: the m = 1 one-base intron traceback now keeps
  its duration component (engels-0036; regression asserts component and both
  scores for lengths 1 and 2 in both decoders) and `DurationMixture` rejects
  unequal component counts across phases, non-positive or non-finite
  weights, q outside (0, 1) and missing phase rows (stalin-0043). Added
  proposal 3.1 edge partials to the oracle via `EdgePrior`: E0(q) entry at
  boundary 0 with the normalized 1/3 * 4**-p phase/prefix prior, J(E(q), r)
  residual-intron entry with no invented donor and log pi, exit at boundary
  n from S/E/I/T with the exit prior and no stop, acceptor or (1-q); E0 must
  consume its first base as CDS and J cannot exit at n, so no-CDS paths are
  discarded by construction; traceback sets `partial_5`/`partial_3`. Added
  `model/grammar/delayed.py`, the 3.2 delayed-entry recurrence: only tails
  active, pending donors in an m-deep ring buffer, rolling per-phase
  emission window, censored pending donors exit at n, and traceback expands
  each tail entry back into I(c, 1..m-1) so chains compare exactly. Checks:
  exhaustive agreement of partition, Viterbi score and traced chains between
  expanded and delayed entry on every ACGT sequence of length 1 to 5 (m = 2,
  R = 2, edges on) plus 120 random lattices (m 1 to 4, R 1 to 3, IUPAC
  ambiguity, tables 1 and 6, edges on/off); brute-force path enumeration now
  includes initial and terminal edge weights. 16 tests, 16 s. Not yet done:
  checkpoint and chunk-seam replay (3.3), reverse-complement agreement, the
  3.6 admission audit and `docs/design/admission-audit.md`. Next tick:
  seam replay on the delayed decoder (checkpoint of active scores plus the
  pending ring and rolling sums, resume across an intron), then the audit.
- 2026-09-15 marx: -> in_progress.
- 2026-09-15 marx (18:0xZ): pushed cd52669 and bfc3ca4 on `work/T-human-013-marx`
  (PR #31); 31 tests, about 10 s, standard library only. Fixed the two
  confirmed interim findings: the delayed decoder's rolling emission window
  now keeps a finite sum plus a per-phase count of masked (-inf) positions,
  so a mask inside a donor's mandatory interval forbids that entry and the
  sum recovers once the mask leaves the window (engels-0037; regression at
  m 1/2/20, the 72-configuration pinned-chain sweep against the closed form,
  and random sparse masks with decoder agreement). Edge policy (stalin-0044):
  E0 is never terminal, so empty input has partition 0 and no chain; J must
  consume one observed intronic base before closing, so no zero-length
  residual intron exists and CDS-from-the-edge is priced once, by E0; the
  same rule in both recurrences. Score channels are validated at the API
  boundary (finite or -inf; NaN and +inf raise). Added `model/grammar/strand.py`
  (IUPAC reverse complement, oriented-to-genomic feature map with
  chain-derived GFF3 phase, CDS GFF3 rows with partial/uncertain attributes)
  with the 3.4 cases planted on the minus strand and re-read from genomic
  coordinates. Added proposal 3.3 seam replay: the delayed recurrence runs in
  chunks between interior seams from a `Checkpoint` (active scores, back
  pointers into the seam, pending donor ring, emission window with sum and
  mask count); seams grant no partial entry or exit and a traceback jump
  across a seam resumes from the donor's true boundary. Checked at every
  seam position of the phase-1 acceptance intron, all seams at once, and 150
  random lattices with random seam sets against the reference. Not yet done:
  the 3.6 admission audit and `docs/design/admission-audit.md`, and the
  GFF3 emission of explicit codon features. Next tick: the audit on the ten
  train species (fetch scripts and manifests with checksums, no sequence
  data committed), then move to `review`.
- 2026-09-15 marx: -> in_progress.
- 2026-09-15 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-15 marx (19:0xZ): pushed c8add41 on `work/T-human-013-marx` (PR #31,
  marked ready for review); 36 tests, about 12 s, standard library only.
  Proposal 3.3 bounded-memory traceback replay (engels-0038, stalin-0045):
  the forward pass keeps only the block checkpoints, `viterbi` recomputes
  one block at a time from its checkpoint while tracing back, at most one
  replay per block and one block of back pointers held; regressions count
  chunk calls and the largest held set, and pin multi-block donor jumps
  against the reference. Section 3.6 admission audit run on all ten train
  species with `model/labels/admission.py` (metadata audit, then streaming
  FASTA checks): manifests and summaries under `model/labels/manifests/`,
  counts, reasons, checksums and findings in `docs/design/admission-audit.md`.
  The 4.4 rows reproduce exactly for the four species it covered. Admitted
  representatives: mouse 21,694 of 22,183; zebrafish 19,574 of 28,402
  (7,959 masked on the exception tag `annotated by transcript or proteomic
  data`, only 156 of them failing an ORF check; frog and maize show the
  same pattern); fly 12,974; worm 19,784; Arabidopsis 27,220; maize 31,804;
  yeast 5,858; Neurospora 9,722; Dictyostelium 12,937. Finite legal
  numerator: constructive for every admitted chain, and checked with the
  delayed decoder on ten spliced admitted genes per species (100 windows,
  0 failures, chains recovered exactly). Adapter fixtures of 3.6 on both
  strands in `tests/test_label_admission.py`. All deliverables are on the
  branch; status -> review. Next tick: answer review requests only.
