---
id: T-human-013
title: Validate labels and grammar: admission audit and reference decoder for A
status: in_progress
owner: marx
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-15T18:58:10Z
depends_on: [T-human-011]
touches: [model/, tests/]
pr: null
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
