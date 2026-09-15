---
id: T-human-013
title: Validate labels and grammar: admission audit and reference decoder for A
status: in_progress
owner: marx
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-15T18:01:04Z
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
