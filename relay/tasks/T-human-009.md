---
id: T-human-009
title: Cost baseline of existing tools and a compute budget for ours
status: review
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-12T14:52:25Z
depends_on: [T-human-006]
touches: [docs/cost-baseline.md]
pr: https://github.com/nekrut/annotation/pull/20
---

## Goal

Put numbers on "current tools are slow and expensive", and set the budget
the new model must fit in.

## What to produce

- `docs/cost-baseline.md`: for EGAPx, Tiberius, BRAKER3, AUGUSTUS, Helixer,
  and GeneMark-ETP at minimum: documented or measured wall clock, CPU or GPU
  hours, peak memory, and failure rate per genome, normalized per megabase.
  Use the Galaxy EGAPx analysis in nekrut/scalingPaper (1,409 jobs, 290
  genomes, about 416,000 CPU-hours, 45% failure, 24% CPU efficiency) as one
  primary source and cite others.
- Where you can, run one tool on one small benchmark genome on a laptop and
  record the actual numbers. If a tool cannot run within the laptop
  constraint, say so; that is itself a finding. Measuring EGAPx or Tiberius
  end to end will need the cluster: post an `alert` with an estimate, as
  the charter describes, and run there once granted.
- A target budget: what "fast" must mean for the new model in seconds per
  megabase on one CPU core and on one consumer GPU, and the parameter and
  memory ceiling, with justification.

## Definition of done

Pull request merged after one `review`. Every number has a source.

## Log

- 2026-09-09 human: created.
- 2026-09-10 marx: claimed.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: claimed once T-human-006 went `done` (the draft had been
  built on `work/T-human-009-marx` over the previous ticks, see marx-0026,
  0031, 0032). This tick: merged `main` into the branch (clean), rewrote the
  status line, cited the merged synthesis instead of the open PR, added the
  Helixer consumer-GPU memory note from `benchmark/validation/README.md`
  (default batch 32 exhausts 16 GB at the vertebrate window; the validation
  declarations for Tiberius and Helixer record no time or memory, checked),
  posted the cluster `alert` 20260910T115307Z-marx-0033 for EGAPx on
  *C. intestinalis* (32 CPU, 256 GB, 70 to 100 CPU-h) and Tiberius on the
  panel (one GPU, about 3 GPU-h), and opened PR #20
  (https://github.com/nekrut/annotation/pull/20, head 86a9091, 7 files,
  +532). Status -> review. Definition of done needs one `review` from
  another agent and a merge. Next: answer review comments; run the two
  cluster rows only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: lenin's `review` (lenin-0033) received, verdict accept
  after six fixes; all six numbers re-derived here and confirmed (panel
  13.959 Gb; covered Tiberius species 10.602 Gb; `cds_fraction_pct` sums
  145,477 human transcripts). Pushed 01985d9 on `work/T-human-009-marx`
  (PR #20): panel 14 Gb and 1.9 h; Tiberius row lists the eleven species,
  10.60 Gb, 5.6 GPU-h at 1.9 GPU-s/Mb; §5.3 conclusion 2 now cites a
  measured genomic CDS union (new `docs/cost-baseline/cds_union.py` and
  `cds_union.tsv`: human 1.34%, mouse 1.38%, maize 1.96%, yeast 72.5%,
  md5-checked RefSeq GFFs) instead of the isoform-summed column; one-core
  peak stated with two FMA ports and conclusion 1 shown to hold at 9e10;
  GPU ceiling normalized to tensor throughput (4090 165 vs A100 312 TFLOPS
  BF16, sources cited); CPU ceiling bound to the §3.1 runner with the
  AUGUSTUS ratio as the portable form. §6 records lenin's two suggested rows
  (gigabase AUGUSTUS: 35 to 100 CPU-h, outside one tick; timed Tiberius and
  Helixer reruns: need the RTX 5080 box) as not measured, with reasons. The
  corrected cluster estimate is restated in marx's note this tick because
  marx-0033 cannot be edited. Status stays `review`; the merge by `human`
  is what remains. Next: any further review comments; cluster rows only on
  a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: lenin's re-review (lenin-0034) received, verdict accept,
  one non-blocking finding: the human CDS union counted GRCh38.p14's alt
  and patch scaffolds, which are not in the panel's `genome_bp`. Landed it
  now rather than deferring, since the PR is unmerged and the fix is small:
  `cds_union.py` fetches the NCBI assembly report, keeps only
  primary-assembly nuclear seqids (roles assembled-molecule, unlocalized-
  and unplaced-scaffold; unit not non-nuclear) and checks their lengths sum
  to `genome_bp` (exact for all four species); `--all` restores the old
  count. Re-derived: human primary union 36,126,376 bp = 1.17% (lenin's
  36,137,717 includes chrM; the 11,341 bp difference is chrM CDS); the
  unfiltered run still gives 41,405,486, so primary plus excluded
  (5,279,110) reproduces the old table exactly. Isoform inflation on the
  primary assembly is 7.6x (273.37/36.13), not the 8.3x in lenin-0034.
  Other rows: yeast 72.35% (chrM dropped), mouse 1.38%, maize 1.95%.
  `cds_union.tsv` regenerated with two new columns (`sequences`,
  `excluded_cds_union_bp`); §3.4 explains the filter; §5.3 conclusion 2
  reads 1.2% human, 8.8% isoform sum, 7.6x, candidate stage 1.7 to 4.3x.
  Pushed as one commit on `work/T-human-009-marx` (PR #20). Status stays
  `review`; nothing else in lenin-0034 asked for a change; the merge by
  `human` is what remains. Next: further review comments if any; cluster
  rows only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: lenin's third review (lenin-0035, accept, nothing to
  change) received. Used the tick on the row §6 named next: AUGUSTUS 3.5.0
  (`human` set) on human chromosome 21, NC_000021.9, 46.71 Mb (40.09
  non-N), eight `--predictionStart/--predictionEnd` windows of 5.84 Mb,
  four in parallel under GNU time on the §3.1 runner. 1,154 CPU-s, 347 s
  wall, 1.65 GB peak, 308 loci against 238 in RefSeq; nucleotide F1 0.693,
  locus F1 0.575, transcript F1 0.110 (`benchmark/score.py --seqids`).
  25 CPU-s/Mb is 2.3 to 6.8x under the five small genomes, so the budget's
  flat-in-genome-size assumption failed in the GHMM's favour; a
  `--softmasking=0` control on one window (105 vs 169 CPU-s, 140 vs 42
  genes) rules masking out as the cause; read overhead per process is
  0.5 CPU-s. §5.2's CPU ceiling restated: 15 CPU-s/Mb is 1/11 of AUGUSTUS
  on *S. pombe* but only 1.6x under it on chr21, so on mammals the model
  buys accuracy, not CPU. The whole-genome FASTA transfer was cut off by
  the proxy after chr21 (checksum unconfirmed; chr21 extracted at its full
  length; GFF md5 ok), recorded in §3.2. Pushed 781df39 on
  `work/T-human-009-marx` (PR #20): `run_augustus_windows.sh`,
  `augustus-Homo_sapiens-chr21.yaml`, `measured.tsv` row, §1.1, §3.2 to
  3.4, §5.1 to 5.3, §6. Status stays `review`; merge by `human` remains.
  Next: one maize chromosome by the same script if a tick is free before
  the merge; cluster rows only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: engels's `review` (engels-0031) received: three requested
  corrections, all landed on `work/T-human-009-marx` (PR #20, head f9ddbfa,
  three commits 3a0bfbf, b5033fa, f9ddbfa). (1) `run_augustus_windows.sh`
  now runs under `set -euo pipefail`, checks every scheduled window for a
  GNU time exit status of 0 before combining, exits 1 with no
  `augustus.gff3` and no success summary otherwise, treats a zero-gene
  window as success, records the FASTA md5, and writes `windows.tsv` (exit
  status, CPU, RSS, genes per window). Tested: a shim that fails window 1
  of 4 (wrapper exit 1, no combined GFF, windows.tsv shows status 3), a
  2 kb all-N sequence (2/2 windows ok, genes=0), yeast chrI in 4 windows
  (82 genes, 4/4 ok). (2) New `window_cuts.py` reproduces the seven cut
  coordinates and counts reference loci (scorer's filter, 238 on chr21)
  whose merged CDS span contains one: 2, not "at most seven" (CXADR at
  17,516,244|245, RUNX1 at 35,032,488|489; 4 within 50 kb). The document
  now gives the coordinates, the filter and the count, says the near-cut
  effect is not measured, and frames the intron-length explanation as a
  hypothesis for T-human-011. (3) §1.1, §3.3, §5.2 and §6: candidate
  count, masking (scoped to window 1), parameter set versus sequence, and
  the AUGUSTUS spread transferring to the new model are labelled
  hypotheses with the observations kept separate; every whole-genome
  21 CPU-h figure says it assumes uniform per-Mb cost; the idle-threads
  claim is removed. Beyond the asks: reran chr21 under the checked runner
  on a 2.10 GHz instance of the same class (the 2.80 GHz original was
  not available): all 8 windows exit 0, 1,225.5 CPU-s, 359 s wall,
  1.69 GB, and predictions identical to the row gene for gene (308 genes,
  F1 0.69286 / 0.57509 / 0.10989), so the reported run had no failed
  window and the accuracy is reproducible; `windows-Homo_sapiens-chr21.tsv`
  committed, FASTA md5 in the declaration, row unchanged. Status stays
  `review`; merge by `human` remains. Next: further review comments if
  any; cluster rows only on a `decision`.
- 2026-09-10 marx: -> in_progress.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: engels's second `review` (engels-0032, two requested
  edits, one non-blocking edge case) and lenin's fourth (lenin-0036,
  accept at f9ddbfa) received. All three engels items landed as one commit
  febf828 on `work/T-human-009-marx` (PR #20). (1) §3.3: the sentence
  "windowing accounts for little of that" replaced by the measured fact
  (two retained reference-locus CDS spans intersect a cut) and the
  statement that windowing's contribution to the error is unquantified.
  (2) §3.2: the replicate paragraph now says it completed every window
  and reproduced the original predictions and accuracy, that its user CPU
  was 6.16% higher (1,225.45 vs 1,154.3 CPU-s), and that the two runs
  neither characterize timing variability nor recover the original
  workers' exit statuses. (3) `run_augustus_windows.sh` and
  `window_cuts.py` reject a window count whose ceil(length/K) step leaves
  an empty final window (length 10, K=6), with the same rule in both, so
  runner and cut counter always agree; stub-augustus tests: length 10/K=6
  exits 1 before any window runs, length 400/K=4 completes 4/4, chr21's
  seven cuts unchanged. No numbers changed. Status stays `review`; engels's
  acceptance and the merge by `human` remain. Next: nothing further unless
  a reviewer asks; cluster rows only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: lenin's fifth review (lenin-0037) and engels's third
  (engels-0033) received, both accept at febf828 with no change requests
  (lenin's `cuts(length, 0)` exception-ordering remark is explicitly not a
  request; engels's acceptance requests no revision). Checked PR #20: open,
  head febf828, mergeable state clean, base 8700d27. Per review discipline
  nothing was pushed to `work/T-human-009-marx`, so the coordinator merges
  the head both reviewers have seen. Lease renewed; status stays `review`.
  Next: nothing until the merge or a new review comment; the cosmetic
  `window_cuts.py` items fold in only if the branch is touched again for
  another reason; cluster rows (marx-0033) only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: inbox empty; no new review message and no GitHub-side
  review or comment on PR #20 since lenin-0037 and engels-0033 (both accept
  at febf828). Checked PR #20 this tick: open, head febf828, base 8700d27,
  mergeable state clean, 13 commits, no required checks. Per review
  discipline nothing was pushed to `work/T-human-009-marx`. lenin's alert
  lenin-0038 to `human` already states that the merge is the only blocker
  for T-human-011, so no note this tick. Lease renewed; status stays
  `review`. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-10 marx: -> in_progress.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: inbox held one item, lenin's `proposal` lenin-0039 to
  `human` (let T-human-011 start against the accepted §5 budget at
  febf828); it asks nothing of marx and its reading of §5 is consistent
  with the branch, so no answer. No new `review` message and no
  GitHub-side review or comment on PR #20 since lenin-0037 and engels-0033
  (both accept at febf828). Checked PR #20 this tick: open, head febf828,
  base 8700d27, mergeable state clean, 13 commits. Per review discipline
  nothing was pushed to `work/T-human-009-marx`. Lease renewed; status
  stays `review`; no note. Next: nothing until the merge or a new review
  comment; cluster rows (marx-0033) only on a `decision`.
- 2026-09-10 marx: -> in_progress.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick: open,
  head febf828, base 8700d27, mergeable state clean, 13 commits. Per
  review discipline nothing was pushed to `work/T-human-009-marx`. Lease
  renewed; status stays `review`; no note. Next: nothing until the merge
  or a new review comment; cluster rows (marx-0033) only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (21:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 16:53Z. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (22:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 16:53Z. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-10 marx: -> in_progress.
- 2026-09-10 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-10 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (23:51Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 16:53Z. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (00:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (01:55Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (02:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (03:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (04:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (05:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (06:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (07:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (08:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (09:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (10:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (11:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (12:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (13:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (14:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (15:51Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (16:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (17:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (18:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (19:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> in_progress.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (20:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (21:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (22:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-11 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-11 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (23:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> in_progress.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (00:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (01:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (02:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (03:53Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> in_progress.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (04:52Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Per review discipline nothing was pushed
  to `work/T-human-009-marx`. Lease renewed; status stays `review`; no
  note. Next: nothing until the merge or a new review comment; cluster
  rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (05:53Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Only relay ticks (lenin holding on
  T-human-011, marx lease renewal) landed on `main` since the last tick.
  Per review discipline nothing was pushed to `work/T-human-009-marx`.
  Lease renewed; status stays `review`; no note. Next: nothing until the
  merge or a new review comment; cluster rows (marx-0033) only on a
  `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. No new `review` message and no GitHub-side
  review, review thread or comment on PR #20 since lenin-0037 and
  engels-0033 (both accept at febf828). Checked PR #20 this tick (06:55Z):
  open, head febf828, base 8700d27, mergeable state clean, 13 commits,
  last updated 2026-09-10 16:53Z. Only relay ticks (lenin holding on
  T-human-011, marx lease renewal) landed on `main` since the last tick.
  Per review discipline nothing was pushed to `work/T-human-009-marx`.
  Lease renewed; status stays `review`; no note. Next: nothing until the
  merge or a new review comment; cluster rows (marx-0033) only on a
  `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. Checked PR #20 this tick (07:51Z): open,
  head febf828, base 8700d27, mergeable state clean, 13 commits, last
  updated 2026-09-10 16:53Z; no new review, thread or comment on GitHub and
  no new `review` message in the relay since lenin-0037 and engels-0033
  (both accept at febf828). Only relay ticks (lenin holding on T-human-011,
  marx lease renewal) landed on `main` since the last tick. Per review
  discipline nothing was pushed to `work/T-human-009-marx`. Lease renewed;
  status stays `review`; no note. Next: nothing until the merge or a new
  review comment; cluster rows (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> in_progress.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. Checked PR #20 this tick (08:51Z): open,
  head febf828, base 8700d27, mergeable state clean, 13 commits, last
  updated 2026-09-10 16:53Z; zero GitHub reviews and zero comments on the
  PR, and no new `review` message in the relay since lenin-0037 and
  engels-0033 (both accept at febf828). Only relay ticks (lenin, engels and
  stalin holding on T-human-011, marx lease renewal) landed on `main` since
  the last tick. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. Checked PR #20 this tick (09:52Z): open,
  head febf828, base 8700d27, mergeable state clean, 13 commits, last
  updated 2026-09-10 16:53Z; zero GitHub reviews and zero comments on the
  PR, and no new `review` message in the relay since lenin-0037 and
  engels-0033 (both accept at febf828). Only relay holding ticks (lenin,
  engels, stalin waiting on T-human-011, marx lease renewal) landed on
  `main` since the last tick. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. Checked PR #20 this tick (10:53Z): open,
  head febf828, base 8700d27, mergeable state clean, 13 commits, last
  updated 2026-09-10 16:53Z; zero GitHub reviews and zero comments on the
  PR, and no new `review` message in the relay since lenin-0037 and
  engels-0033 (both accept at febf828). Only relay holding ticks (lenin,
  engels, stalin waiting on T-human-011, marx lease renewal) landed on
  `main` since the last tick. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. Checked PR #20 this tick (11:52Z): open,
  head febf828, base 8700d27, mergeable state clean, 13 commits, last
  updated 2026-09-10 16:53Z; zero GitHub reviews and zero comments on the
  PR, and no new `review` message in the relay since lenin-0037 and
  engels-0033 (both accept at febf828). Only relay holding ticks (lenin,
  engels, stalin waiting on T-human-011, marx lease renewal) landed on
  `main` since the last tick. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
- 2026-09-12 marx: -> review (https://github.com/nekrut/annotation/pull/20).
- 2026-09-12 marx: inbox empty. Checked PR #20 this tick (12:53Z): open,
  head febf828, base 8700d27, mergeable state clean, 13 commits, last
  updated 2026-09-10 16:53Z; zero GitHub reviews and zero comments on the
  PR, and no new `review` message in the relay since lenin-0037 and
  engels-0033 (both accept at febf828). Only relay holding ticks (lenin,
  engels, stalin waiting on T-human-011, marx lease renewal) landed on
  `main` since the last tick. Per review discipline nothing was pushed to
  `work/T-human-009-marx`. Lease renewed; status stays `review`; no note.
  Next: nothing until the merge or a new review comment; cluster rows
  (marx-0033) only on a `decision`.
