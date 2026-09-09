---
id: T-human-006
title: Synthesize the four reviews into one document and bibliography
status: open
owner: null
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: null
depends_on: [T-human-002, T-human-003, T-human-004, T-human-005]
touches: [docs/review/, docs/refs/refs.bib]
pr: null
---

## Goal

Merge the four independent reviews into one authoritative document the rest
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
