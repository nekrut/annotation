---
id: T-human-006
title: Synthesize the four reviews into one document and bibliography
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-10T05:24:23Z
depends_on: [T-human-002, T-human-003, T-human-004, T-human-005, T-human-012]
touches: [docs/review/, docs/refs/refs.bib, scripts/review/]
pr: https://github.com/nekrut/annotation/pull/16
---

## Goal

Merge the five independent reviews into one authoritative document the rest
of the project cites, without losing the disagreements.

## What to produce

- `docs/review/README.md`: the synthesized review, structured like the
  individual ones, with a single merged publications table and repository
  table.
- `docs/refs/refs.bib`: deduplicated bibliography; every entry has a DOI or
  preprint URL.
- `docs/review/disagreements.md`: every point where the four reviews
  reached different conclusions, what each said, and either a resolution
  with evidence or an open question for the coordinator.
- `docs/review/candidates.md`: the shortlist of ideas for the new model
  distilled from the four opinion sections, ranked, with the reasoning.

## Definition of done

Pull request merged after a `review` from at least one agent whose own
review was not the primary source of a section. A `note` to `all` summarizing
the top five findings.

## Log

- 2026-09-09 human: created.
- 2026-09-10 lenin: claimed.
- 2026-09-10 lenin: merged all five reviews. Wrote docs/review/README.md
  (merged publications and software tables, data sources, 14 consolidated
  failure modes, joint design position), docs/review/disagreements.md
  (factual conflicts resolved with evidence, 6 design disagreements, 4 open
  questions for the coordinator) and docs/review/candidates.md (8 candidates
  ranked by decision-relevance per unit compute). docs/refs/refs.bib: 111
  unique works deduplicated by DOI from 206 source entries. Added
  scripts/review/ (stdlib only) so every count is reproducible, and extended
  `touches` to cover it. Two verification passes: Crossref resolution of all
  111 DOIs (9 resolve to unrelated papers or nothing) and GitHub re-query of
  the 23 disputed repository rows (three slugs were called "Helixer"; one is
  a 0-star fork). -> review, PR #16.
  Next: merge the 44 engels/stalin annexes' own bibliography and repository
  additions, which this pass indexed but did not absorb (README.md §7).
