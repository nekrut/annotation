---
id: T-human-013
title: Validate labels and grammar: admission audit and reference decoder for A
status: review
owner: marx
created_by: human
created: 2026-09-15T14:33:49Z
lease_until: 2026-09-17T00:53:18Z
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
- 2026-09-15 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-15 marx (20:0xZ): pushed cc0955a on `work/T-human-013-marx` (PR #31);
  44 tests, about 13 s, Python 3.11.15. Answered the change requests of
  engels-0039 and stalin-0046, all with both-strand regressions: (1) a declared
  partial end is admissible only where that end touches the sequence edge in
  transcriptional orientation, both when both are declared; an interior range
  declaration masks the chain as `partial_unlocated` and makes the junction at
  that row unknown; (2) every benchmark-accepted transcript is FASTA-audited,
  the auxiliary catalog reports annotated / retained / unknown separately and a
  site is retained only on a reliable observation (end not partial, coordinates
  in range, initiator or stop present on that transcript); (3) `transl_except`
  on the transcript row masks; (4) explicit `start_codon` / `stop_codon` GFF3
  features from the decoded chain, split across introns with the continuation
  phase, omitted on a censored end, terminal stop kept inside the CDS;
  (5) transcript identity is (seqid, strand, id) in every cache, component and
  membership set, invariant to FASTA order. All ten manifests, summaries and
  numerator checks regenerated from the same checksummed inputs (100 windows,
  0 failures): admitted counts unchanged in every species; fly loses the nine
  inherited `topology` / `no_initiator` / `internal_stop` / `phase_inconsistent`
  flags (539 to 530 topology, as stalin-0046 computed); interior partials no
  longer count as chain-end partials; 4 to 36 junctions per vertebrate species
  and about 30 starts/stops per species move to unknown.
  `docs/design/admission-audit.md` updated with the new tables, checksums and a
  paragraph on what changed. Next tick: answer further review requests only.
- 2026-09-15 marx: -> in_progress.
- 2026-09-15 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-15 marx (21:0xZ): pushed de6d599 on `work/T-human-013-marx` (PR #31);
  45 tests, about 15 s, Python 3.11.15. Answered the one open change request
  (engels-0040, reproduced by stalin-0047; lenin-0049 concurred with requests
  already answered in cc0955a): the FASTA audit read the terminal stop off the
  last complete codon, so a complete stop followed by leftover bases passed the
  end check and the annotated end stayed a retained auxiliary stop. The stop
  must now be the last three observed bases in frame; leftover bases fail the
  end check (`no_stop` with `seq_frame_length`) and every complete codon, an
  earlier stop included, is `internal_stop`. Regressions on both strands,
  split and unsplit, non-stop remainder, in-frame control and a 3'-partial
  edge control (still admitted, stop unknown). Re-fetched all ten train
  genomes (every GFF and FASTA md5 matches the documented inputs) and
  regenerated all ten manifests: nine decompress byte-identical and keep
  their committed files and checksums; one representative changes in the
  whole panel, zebrafish rna-NM_001100045.1 (gene-b3gnt5b, already masked),
  which gains `no_stop` and `internal_stop` and moves its end from retained
  to unknown (stops 33,263/93 to 33,262/94; `internal_stop` 153 to 154,
  `no_stop` 15 to 16). No admitted set changed. Zebrafish numerator check
  rerun on the new manifest: same ten spliced genes, identical scores, 0
  failures. `docs/design/admission-audit.md` updated (rule text, both tables,
  the 168/131 prose to 154/112, checksum row, a paragraph on the change).
  Next tick: answer further review requests only; once a reviewer accepts
  de6d599, nothing more is pushed to the branch.
- 2026-09-15 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-15 marx (21:5xZ): lease renewed, status stays `review`. All three
  reviewers accepted head de6d599 (lenin-0050, engels-0041, stalin-0048) and no
  change request is open, so nothing further is pushed to
  `work/T-human-013-marx`; PR #31 waits on the coordinator's merge. Read
  gagarin-0007 (tiberius-ten finished for T-human-009; egapx-ciona failed on a
  missing `-lc` cache directory) and will push the one-line sbatch fix on
  `work/T-human-009-marx` only if the coordinator answers gagarin-0008 with
  option 1. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; otherwise review or prepare T-human-014.
- 2026-09-15 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-15 marx (22:5xZ): lease renewed, status stays `review`. No new
  messages since note 0057 and no open change request; PR #31 head de6d599
  (accepted by lenin-0050, engels-0041, stalin-0048) still merges cleanly onto
  main (47 files, 5,075 insertions, no conflicts). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) is still unanswered, so `work/T-human-009-marx` is untouched.
  Next tick: nothing on this branch unless a reviewer or the coordinator
  asks; after the merge, review or prepare T-human-014.
- 2026-09-15 marx: -> in_progress.
- 2026-09-15 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-15 marx (23:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057 and no open change request. PR #31 head
  de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is unchanged,
  still open, and still merges cleanly onto origin/main. Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (00:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057 and no open change request. PR #31 head
  de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is unchanged,
  still open, and still merges cleanly onto origin/main (checked with
  `git merge-tree`). Nothing pushed to `work/T-human-013-marx`; no note
  posted. gagarin-0008 (egapx-ciona cache directory) still awaits the
  coordinator, so `work/T-human-009-marx` is untouched. Next tick: nothing
  on this branch unless a reviewer or the coordinator asks; after the merge,
  review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (01:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057 and no open change request. PR #31 head
  de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is unchanged,
  still open, and still merges cleanly onto origin/main (`git merge-tree`
  after unshallowing the clone; main has moved only under `relay/` since
  the PR base 986722e). Nothing pushed to `work/T-human-013-marx`; no note
  posted. gagarin-0008 (egapx-ciona cache directory) still awaits the
  coordinator, so `work/T-human-009-marx` is untouched. Next tick: nothing
  on this branch unless a reviewer or the coordinator asks; after the merge,
  review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (02:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057 and no open change request. PR #31 head
  de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is unchanged,
  still open, and still merges cleanly onto origin/main (`git merge-tree`
  on the unshallowed clone; main has moved only under `relay/` since the
  PR base 986722e). Nothing pushed to `work/T-human-013-marx`; no note
  posted. gagarin-0008 (egapx-ciona cache directory) still awaits the
  coordinator, so `work/T-human-009-marx` is untouched. Next tick: nothing
  on this branch unless a reviewer or the coordinator asks; after the merge,
  review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (03:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057 and no open change request. PR #31 head
  de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is unchanged,
  still open, and still merges cleanly onto origin/main at 209298c
  (`git merge-tree`; main has moved only under `relay/` since the PR base
  986722e). Nothing pushed to `work/T-human-013-marx`; no note posted.
  gagarin-0008 (egapx-ciona cache directory) still awaits the coordinator,
  so `work/T-human-009-marx` is untouched. Next tick: nothing on this
  branch unless a reviewer or the coordinator asks; after the merge, review
  or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (04:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057 and no open change request. PR #31 head
  de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is unchanged,
  still open (last PR update 2026-09-15T21:01Z), and still merges cleanly
  onto origin/main at 3bb7773 (`git merge-tree` after unshallowing the
  clone; main has moved only under `relay/` since the PR base 986722e).
  Nothing pushed to `work/T-human-013-marx`; no note posted. gagarin-0008
  (egapx-ciona cache directory) still awaits the coordinator, so
  `work/T-human-009-marx` is untouched. Next tick: nothing on this branch
  unless a reviewer or the coordinator asks; after the merge, review or
  prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (05:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review or comment on PR #31,
  and no open change request. PR #31 head de6d599 (accepted by lenin-0050,
  engels-0041, stalin-0048) is unchanged, still open (last PR update
  2026-09-15T21:01Z), and still merges cleanly onto origin/main at 265f7f9
  (`git merge-tree` after unshallowing; main has moved only under `relay/`
  since the PR base 986722e). Nothing pushed to `work/T-human-013-marx`;
  no note posted. gagarin-0008 (egapx-ciona cache directory) still awaits
  the coordinator, so `work/T-human-009-marx` is untouched. Next tick:
  nothing on this branch unless a reviewer or the coordinator asks; after
  the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (06:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  f688092 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (07:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  0f18599 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (08:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  795a5ba (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (09:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  3287601 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (10:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review, thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  55a3c08 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (11:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review, thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  646b919 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (12:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review, thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  1e76cb0 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (13:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review, thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  6b26e59 (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (14:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message since marx-0057, no GitHub review, thread or comment on
  PR #31, and no open change request. PR #31 head de6d599 (accepted by
  lenin-0050, engels-0041, stalin-0048) is unchanged, still open (last PR
  update 2026-09-15T21:01Z), and still merges cleanly onto origin/main at
  4f8f3ae (`git merge-tree` after unshallowing; main has moved only under
  `relay/` since the PR base 986722e). Nothing pushed to
  `work/T-human-013-marx`; no note posted. gagarin-0008 (egapx-ciona cache
  directory) still awaits the coordinator, so `work/T-human-009-marx` is
  untouched. Next tick: nothing on this branch unless a reviewer or the
  coordinator asks; after the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (15:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since marx-0057, no GitHub review,
  thread or comment on PR #31, and no open change request. lenin-0051
  (alert to human, task T-human-014) asks the coordinator to merge PR #31
  since it blocks all of Phase 4; noted, nothing for marx to do. PR #31
  head de6d599 (accepted by lenin-0050, engels-0041, stalin-0048) is
  unchanged, still open (last PR update 2026-09-15T21:01Z), and still
  merges cleanly onto origin/main at e599af3 (`git merge-tree` after
  unshallowing; main has moved only under `relay/` since the PR base
  986722e). Nothing pushed to `work/T-human-013-marx`; no note posted.
  gagarin-0008 (egapx-ciona cache directory) still awaits the coordinator,
  so `work/T-human-009-marx` is untouched. Next tick: nothing on this
  branch unless a reviewer or the coordinator asks; after the merge, review
  or prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (16:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since marx-0057, no GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), and no
  open change request. lenin-0051 (alert to human) still awaits the
  coordinator's merge. PR #31 head de6d599 (accepted by lenin-0050,
  engels-0041, stalin-0048) is unchanged, still open, and still merges
  cleanly onto origin/main at 37ec9e2 (`git merge-tree` after unshallowing;
  main has moved only under `relay/` since the PR base 986722e). Nothing
  pushed to `work/T-human-013-marx`; no note posted. gagarin-0008
  (egapx-ciona cache directory) still awaits the coordinator, so
  `work/T-human-009-marx` is untouched. Next tick: nothing on this branch
  unless a reviewer or the coordinator asks; after the merge, review or
  prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (17:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since marx-0057, no GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), and no
  open change request. lenin-0051 (alert to human) still awaits the
  coordinator's merge. PR #31 head de6d599 (accepted by lenin-0050,
  engels-0041, stalin-0048) is unchanged, still open, and still merges
  cleanly onto origin/main at 89c7c0e (`git merge-tree` after unshallowing;
  main has moved only under `relay/` since the PR base 986722e). Nothing
  pushed to `work/T-human-013-marx`; no note posted. gagarin-0008
  (egapx-ciona cache directory) still awaits the coordinator, so
  `work/T-human-009-marx` is untouched. Next tick: nothing on this branch
  unless a reviewer or the coordinator asks; after the merge, review or
  prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (18:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since marx-0057, no GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), and no
  open change request. lenin-0051 (alert to human) still awaits the
  coordinator's merge. PR #31 head de6d599 (accepted by lenin-0050,
  engels-0041, stalin-0048) is unchanged, still open, and still merges
  cleanly onto origin/main at 6466fc2 (`git merge-tree` after unshallowing;
  main has moved only under `relay/` since the PR base 986722e). Nothing
  pushed to `work/T-human-013-marx`; no note posted. gagarin-0008
  (egapx-ciona cache directory) still awaits the coordinator, so
  `work/T-human-009-marx` is untouched. Next tick: nothing on this branch
  unless a reviewer or the coordinator asks; after the merge, review or
  prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (19:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since marx-0057, no GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), and no
  open change request. Since last tick main gained only the other agents'
  heartbeat commits (engels, stalin, lenin, gagarin all report inbox empty
  and waiting on this merge). lenin-0051 (alert to human) still awaits the
  coordinator's merge. PR #31 head de6d599 (accepted by lenin-0050,
  engels-0041, stalin-0048) is unchanged, still open, and still merges
  cleanly onto origin/main at 659c9e1 (`git merge-tree` after unshallowing;
  main has moved only under `relay/` since the PR base 986722e). Nothing
  pushed to `work/T-human-013-marx`; no note posted. gagarin-0008
  (egapx-ciona cache directory) still awaits the coordinator, so
  `work/T-human-009-marx` is untouched. Next tick: nothing on this branch
  unless a reviewer or the coordinator asks; after the merge, review or
  prepare T-human-014.
- 2026-09-16 marx (20:5xZ): lease renewed, status stays `review`. Inbox had
  one message, lenin-0052, a re-review of de6d599 confirming lenin-0049's
  request-changes (transcript-id collision) is resolved and accepting; it
  is a `review`, not a question, so no answer is owed. No GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), PR still
  open, head de6d599 unchanged, and it still merges cleanly onto
  origin/main (`git merge-tree`; zero non-relay commits on main since the
  PR base 986722e). All reviewers (lenin-0050/0052, engels-0041,
  stalin-0048) accept the current head and no change request is open, so
  nothing pushed to `work/T-human-013-marx`; no note posted. lenin-0051
  (alert to human) and gagarin-0008 (egapx-ciona cache directory) still
  await the coordinator; `work/T-human-009-marx` untouched. Next tick:
  nothing on this branch unless a reviewer or the coordinator asks; after
  the merge, review or prepare T-human-014.
- 2026-09-16 marx: -> in_progress.
- 2026-09-16 marx (21:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since lenin-0052, no GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), PR still
  open, head de6d599 unchanged, and it still merges cleanly onto
  origin/main at 4c3fd83 (`git merge-tree` after unshallowing; zero
  non-relay commits on main since the PR base 986722e). All reviewers
  (lenin-0050/0052, engels-0041, stalin-0048) accept the current head and
  no change request is open, so nothing pushed to `work/T-human-013-marx`;
  no note posted. lenin-0051 (alert to human) and gagarin-0008
  (egapx-ciona cache directory) still await the coordinator;
  `work/T-human-009-marx` untouched. Next tick: nothing on this branch
  unless a reviewer or the coordinator asks; after the merge, review or
  prepare T-human-014.
- 2026-09-16 marx: -> review (https://github.com/nekrut/annotation/pull/31).
- 2026-09-16 marx (22:5xZ): lease renewed, status stays `review`. Inbox
  empty; no message addressed to marx since lenin-0052, no GitHub review,
  thread or comment on PR #31 (last PR update 2026-09-15T21:01Z), PR still
  open, head de6d599 unchanged, and it still merges cleanly onto
  origin/main at 13c2581 (`git merge-tree` after unshallowing: merge base
  is the PR base 986722e, zero conflict hunks, zero non-relay commits on
  main since). All reviewers (lenin-0050/0052, engels-0041, stalin-0048)
  accept the current head and no change request is open, so nothing pushed
  to `work/T-human-013-marx`; no note posted. lenin-0051 (alert to human)
  and gagarin-0008 (egapx-ciona cache directory) still await the
  coordinator; `work/T-human-009-marx` untouched. Next tick: nothing on
  this branch unless a reviewer or the coordinator asks; after the merge,
  review or prepare T-human-014.
