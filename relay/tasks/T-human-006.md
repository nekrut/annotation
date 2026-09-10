---
id: T-human-006
title: Synthesize the four reviews into one document and bibliography
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-10T12:06:00Z
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
- 2026-09-10 lenin: annex integration, the gap README 7 item 1 has carried
  since the first pass. Head b648f75. New scripts/review/merge_annex_refs.py
  scans the 49 messages engels or stalin filed against T-human-003/T-human-004,
  extracts every DOI and subtracts the five artifact bibliographies: **15 works
  are cited only in the annexes** (10 engels, 7 stalin, 2 both), all 15 resolve
  in Crossref. The covers' own counts (8 engels, 12 stalin) overlap and five of
  stalin's twelve were already merged from other reviews, so 15 is the measured
  union, not 20. merge_refs.py now reads docs/review/annex-refs.bib as a sixth
  source; refs.bib is 111 -> 126 works, 110 crossref-ok, the same 9 bad DOIs
  (all trotsky's). check_conflicts.py scans the annex works too and surfaces two
  new preprint/journal pairs, PhyloCSF++ and ClaMSA - which is engels's own
  annex point that ClaMSA's recovered tables belong to the preprint. Two annex
  repository snapshots were in no base inventory: abacus-gene/paml (stalin 0016)
  and Jstacs/Jstacs, which hosts GeMoMa (engels 0018); both are transcribed into
  docs/review/annex-repos.tsv with their own pin and observation time and merged
  as "<agent> (annex)", repos.tsv 59 -> 61. That corrects section 3's claim that
  GeMoMa is not on GitHub at all. README 2.6 lists all 15 works with why the
  annex reached for them; 5.14 now cites Stentor's 15-base introns, ciliate
  context-dependent stop codons and the U12 database rather than pointing at
  "stalin's scope annexes"; 7 item 1 now says what is still not done - reading
  the 44 annexes as evidence rather than as a citation source - and carries
  stalin 0016's PAML finding (codeml overrides cleandata for pairwise runmodes,
  so a multispecies alignment deletes codons clean in the target/informant pair)
  to T-human-010. Also applied the two wording cleanups engels-0028 and
  stalin-0029 both asked for: candidates.md 4 no longer implies that passing an
  equality test shows an encoder ignores its comparative channel, and
  candidates.md 3 / disagreements.md 5.3 scope the ablation's attribution to the
  encoders, informants, species and split it evaluates. Task stays in `review`.
  Next: nothing further is planned on this task. It needs the coordinator to
  merge PR 16; reading all 44 annexes as primary evidence is a larger job than
  this task's remaining budget and should be its own task if the coordinator
  wants it.
- 2026-09-10 lenin: applied the three corrections both reviewers of the annex
  merge asked for (engels 20260910T062452Z-engels-0029, stalin
  20260910T064019Z-stalin-0030; both confirm the bibliography/repository merge
  reproduces byte for byte). Head 55e31da. (1) README 2.6 row for
  `flicek2006several` said "AUGUSTUS"; the paper describes MARS, a TWINSCAN
  extension that combines predictions from separate pairwise informants, as
  engels 0017 recorded. Renamed. (2) The U12 row and 5.14's last sentence
  summarized engels's annex 0021 as a warning against reading short gaps as
  U12 introns. Its actual point is that U12-type introns also occur with GT-AG
  termini, so zero recovered AT-AC introns does not establish zero U12
  recovery; both places now state that and say to keep terminal-dinucleotide
  strata separate from independently assigned U2/U12 labels. (3) 10 engels /
  7 stalin / 2 both are inclusive counts; grouping annex-dois.tsv by `agents`
  gives 8 + 5 + 2 (union 15, intersection 2), which I reproduced. Dropped
  "only" from the inclusive figures and gave the exclusive split. 2.6's
  rationale no longer claims each annex work was retrieved by exactly one
  agent (rogic2001evaluation and wei2006ests are shared) nor treats annex
  sourcing as a loss of blindness — the annexes were filed during the
  independent review tasks; the supported limitation is that these works are
  absent from the base bibliographies and their substantive integration is
  still incomplete. Also took engels's non-blocking point: 7 item 1 no longer
  reads as an instruction to T-human-010, which is `done` with a Nei-Gojobori
  counting estimator; the PAML `cleandata` warning is scoped to a future
  `codeml` adapter. Citation checker, conflict checker and the 6 DOI
  serialization tests re-run unchanged (95 cited of 126, 4 conflicts, 6
  preprint/journal pairs). Task stays in `review`.
  Next: nothing further planned. PR 16 needs a coordinator merge. Reading the
  44 annexes as primary evidence remains out of this task's budget and should
  be its own task if the coordinator wants it.
- 2026-09-10 lenin: no change to PR 16; head stays 55e31da. All three
  reviewers now accept it: engels 20260910T072225Z-engels-0030 and stalin
  20260910T073658Z-stalin-0031 close every item they raised against the annex
  merge, and marx 20260910T065543Z-marx-0030 (reviewing the earlier b648f75)
  asked for three things that are already in this head. marx's one remaining
  observation - README line 348, "ran on the fly 124-way", where "among four
  references" would be more precise - is explicitly not a change request, so I
  left the branch alone rather than reset three SHA-pinned acceptances for a
  cosmetic edit. Definition of done is met. Posted an alert to human
  (20260910T080618Z-lenin-0032): T-human-009 depends_on [T-human-006] and
  T-human-011 is two hops behind it, so this merge is the only thing gating
  Phase 2 and 3, and engels, stalin and marx have all reported owning nothing
  for several ticks. The alert asks for one of: merge, a list of what still
  has to change, or a decision on whether the disclosed section 7 item 1 scope
  limit (44 annexes merged as citations, not read as evidence) is accepted or
  becomes its own task. Lease renewed to 10:06Z.
  Next: act on the coordinator's answer. Nothing else planned on this task.
- 2026-09-10 lenin: no change to PR 16; head stays 55e31da, all three
  acceptances still pinned to it. Inbox this tick was marx
  20260910T082949Z-marx-0032, which re-verifies the same head against the
  checkers independently (126 refs entries, 95 cited, 4 conflicts, 6
  preprint/journal pairs, 61 repo rows, 6 DOI tests) and confirms its accept;
  it asks nothing of me. No `decision` has been posted since human-0010
  (02:12Z), so my merge alert 20260910T080618Z-lenin-0032 is still open and
  marx's proposal 20260910T032948Z-marx-0026 to relax T-human-009's
  `depends_on` is still unanswered. Verified the PR still applies to a moved
  `main`: `main` has advanced by relay-only commits since the last merge into
  the branch, and `git merge-tree` on the current bases reports no conflicting
  hunks, so the merge is clean whenever the coordinator takes it. Lease renewed
  to 11:06Z.
  Next: act on the coordinator's answer. Nothing else planned on this task.
- 2026-09-10 lenin: quiet tick, no change to PR 16; head stays 55e31da and all
  three acceptances (engels-0030, stalin-0031, marx-0032) remain pinned to it.
  Inbox empty. No `decision` since human-0010 (02:12Z), so my merge alert
  20260910T080618Z-lenin-0032 and marx's proposal 20260910T032948Z-marx-0026
  (relax T-human-009 `depends_on`) are both still open. Re-checked the branch
  against a moved `main` (4aae218): `git merge-tree --write-tree` returns clean,
  0 conflicts, 18 files / +4993 lines. Posted no note this tick - nothing
  happened that the last one did not already say. Lease renewed to 12:06Z.
  Next: act on the coordinator's answer. Nothing else planned on this task.
