---
id: T-human-015
title: Fit and evaluate comparative support for candidate B over frozen A
status: open
owner: null
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: null
depends_on: [T-human-014]
touches: [model/, docs/design/]
pr: null
---

## Goal

Phase 4 milestone 3. Add candidate B's comparative residual over frozen A
(proposal sections 2 and 4) and fit its support policy.

- Freeze A's weights, decoder and support policy from T-human-014. Train
  only B's residual parameters, K=8 taxa from the alignment sources in
  `docs/data-sources.md`, patristic-distance bias as the geometry
  (Tree-RoPE stays a named comparison contingent on accessible code).
- Fit the seed threshold tau on reserved train chromosomes to the declared
  coverage criteria of section 4.2, meter the real token allocation, and
  choose the tile quota alpha from measured cost, not arithmetic. Freeze
  both before any held-out evaluation.
- Keep the CPU-capped and the gated and full GPU regimes as separate
  benchmark rows; a capped result never stands in for unrestricted B.
- If the support gate fails its recall criterion, report that and retain
  A plus the predeclared full-tile GPU experiment for diagnosis.

GPU work goes through gagarin under the Phase 4 training cap.

## Definition of done

- `model/b/` with the residual encoder, the streaming support rule, and a
  run manifest.
- `docs/design/b-support.md`: measured gate recall and token allocation per
  train species, the chosen tau and alpha with the measurements that chose
  them, and cost rows for each regime in `docs/cost-baseline/measured.tsv`.
- Actual GPU-hours recorded in the task log.
- Pull request, at least one `review`, coordinator merge.

## Log

- 2026-09-15 human: created.

