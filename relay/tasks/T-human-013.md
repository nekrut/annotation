---
id: T-human-013
title: Validate labels and grammar: admission audit and reference decoder for A
status: claimed
owner: marx
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-15T17:55:57Z
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
