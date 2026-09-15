---
id: T-human-017
title: Reconsider candidate C with the encoder result: graph metering and edge head
status: open
owner: null
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: null
depends_on: [T-human-016]
touches: [model/, docs/design/]
pr: null
---

## Goal

Phase 4 milestone 5. Reconsider candidate C, the splice graph (proposal
section 5), in the light of the encoder comparison. Do not start before
T-human-016 has a `decision`.

- Meter true-path survival, vertex/edge/window counts, storage and output
  on train development sequences under the vertex and edge policy of
  section 5.1, before any fitting.
- Fit the edge head only if the metered resource limits hold under the
  Phase 4 cap. Retain the declared A fallback wherever the limits are
  exceeded.
- Compare against the winning arm of T-human-016 with the same decoder
  inputs and split.

## Definition of done

- `docs/design/c-admission.md` with the metered survival and resource
  numbers, and either the edge-head results in `docs/results/` or a clear
  statement that C is not admitted under the cap and why.
- Pull request, at least one `review`, coordinator merge.

## Log

- 2026-09-15 human: created.

