---
id: T-human-006
title: Synthesize the four reviews into one document and bibliography
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-10T07:10:00Z
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
- 2026-09-10 lenin: addressed both reviews of PR #16 (stalin
  20260910T033840Z-stalin-0027, request changes; marx 20260910T035517Z-marx-0027,
  accept with two edits). Head c66d753. (1) TSV serialization: fields are now
  whitespace-flattened and the file is read back and checked six-wide before it
  is written; csv.reader sees 111 records, was 123 of width {1,3,4,6}. Verdicts
  reproduce exactly (95 ok, 6 mismatch, 3 unresolved, 2 not-in-crossref, 1
  no-title, 4 no-doi), so refs.bib is byte-identical. merge_refs.py now refuses
  a malformed TSV instead of pinning a truncated title to a DO NOT CITE notice.
  Added scripts/review/test_verify_dois.py, 6 cases, incl. a newline in a
  registered title. (2) Alignment coverage: replaced the blanket "must be built
  with Cactus" in README §4 and candidates.md §2 with the actual limitation,
  citing dm6/multiz124way, ce11/multiz135way and sacCer3/multiz7way from the
  accepted T-human-008 inventory. (3) candidates.md §4 and disagreements.md §5.4
  now split representation symmetries (require invariance) from changes in
  evidence (measure degradation), as stalin's opinion actually asks. (4) The
  Methods 3.7 supervised-human-input caveat now sits beside the Tiberius human
  comparative metric, and candidates.md §3's one-GPU-day estimate is conditioned
  on prepared features plus a "what it cannot decide" paragraph. (5) marx's two
  requests: disagreements.md §5.1 cross-references candidate #5, and
  check_conflicts.py separates preprint/journal pairs (4 conflicts + 4 pairs,
  was 8). Recorded marx's correction to their own §6 CONTRAST figure and the one
  measured row (AUGUSTUS/S. pombe) in README §7. Task stays in `review`.
  Next: the 44 engels/stalin annexes' bibliography and repository additions,
  still the open gap in README §7 item 1.
- 2026-09-10 lenin: addressed the two follow-up reviews of PR #16 (engels
  20260910T042437Z-engels-0027, stalin 20260910T043806Z-stalin-0028; marx
  20260910T045310Z-marx-0028 accepts). Head 6dc286a. (1) Matched controls:
  candidates.md #4 now has six arms — KA/KS, codon-likelihood features, target
  DNA only, MSA/no-tree, MSA/tree-as-tokens, MSA/tree-as-metric — with the
  held-fixed list spelled out (informants, masks, candidate support, decoder,
  split, capacity), so (c)-(d) measures the alignment, (d)-(e)/(d)-(f) the tree
  given the alignment, and (e)-(f) the HyphAeon claim. (2) Scope: candidate #3,
  its summary-table row and disagreements.md 5.3 now present the Tiberius
  ab-initio-vs-ClaMSA reproduction as prioritization evidence about that
  pipeline on those mammalian regions — not a gate on #4, not a general verdict
  on comparative encoders or other clades, not grounds to revise the charter
  hypothesis — and carry the prepared-feature cost condition. (3) Helixer
  endpoint: the scoring qualification now sits in the README 2.4 row beside the
  F1 numbers and in 5.8 (GffCompare without --strict-match/-e matches identical
  intron chains regardless of terminal position; UTR stripping does not make it
  exact at coding termini), citing engels's endpoint audit. (4) Alignment
  reuse: README 4 and candidates.md #2 now separate a genome absent from the
  alignment (the real gap) from an included genome that is not the reference
  (HAL documents arbitrary-reference/subtree queries and MAF export; multiz MAF
  is reference-anchored, so it is an extraction cost), noting neither operation
  is evidence about held-out-species accuracy. (5) Wording: "may legitimately
  move" per stalin, plus the converse caveat that stable predictions on
  uninformative evidence prove nothing; disagreements.md 5.2 restates engels's
  position as order-of-operations rather than opposition to a parameter budget.
  Merged main into the branch so the cited review messages resolve at the PR
  head. Task stays in `review`.
  Next: the 44 engels/stalin annex bibliography and repository additions, still
  the open gap in README 7 item 1.
