---
id: T-human-002
title: Independent review of gene prediction literature and software (slot 1 of 4)
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T07:15:00Z
depends_on: []
touches: [relay/artifacts/T-human-002/]
pr: null
---

## Goal

One independent, careful review of the state of eukaryotic gene prediction:
the published methods and the software that is actually alive. Four
identical tasks exist (T-human-002 to T-human-005), one per agent. Claim any
open one; hold at most one; do not read the other agents' review artifacts
until yours is in `review`.

## What to produce

Write `relay/artifacts/T-human-002/review.md` with these sections, plus
`relay/artifacts/T-human-002/refs.bib` and `relay/artifacts/T-human-002/repos.tsv`.

1. **Search log.** Queries run, sources used (PMC OA, Europe PMC, bioRxiv,
   arXiv, Semantic Scholar/OpenAlex, GitHub search), dates, counts.
2. **Publications table.** For each method: name, year, DOI, approach class
   (ab initio HMM/GHMM, comparative or alignment-based, evidence-based
   pipeline, deep learning, hybrid), inputs required, clades it was trained
   or evaluated on, reported accuracy with the benchmark and metric used,
   runtime and hardware, cross-species generalization evidence, code link.
   Cover at least: GENSCAN, AUGUSTUS and BRAKER, GeneMark family, SNAP,
   Gnomon and EGAPx, MAKER, N-SCAN and TWINSCAN, CONTRAST, Helixer,
   Tiberius, and any 2023 to 2026 deep learning work you find. Include the
   KA/KS approach (PMC155263) and other simple comparative signals.
3. **Repository inventory** (`repos.tsv`): URL, last commit date, commits in
   the last 12 months, open issues, stars, language, licence, whether the
   README install worked for you in a fresh environment, notes. Start from
   Gaius-Augustus/Tiberius and ncbi/egapx and search outward.
4. **Data sources noted in passing.** Anything the papers used that we
   might reuse: alignments, reference annotations, RNA-seq compendia.
5. **Failure modes of current tools.** With citations: cross-clade
   generalization, intron length sensitivity, GC isochores, resource use.
6. **Your opinion, one page.** What you would build, what inputs it should
   take, what geometry or inductive bias it should encode, what the biggest
   risk is. Be specific and be willing to be wrong.

## Definition of done

All three artifact files exist, every row in the publications table has a
DOI or preprint link, `refs.bib` parses, the review is under 6,000 words
excluding tables, and the task is in `review` with a `note` to `all`
announcing it. The coordinator marks it `done`.

## Log

- 2026-09-09 human: created.
- 2026-09-09 lenin: claimed.
- 2026-09-09 lenin: run 1. Literature pass via Europe PMC REST (20 broad + 25
  targeted title queries + 19 DOI lookups) and a GitHub API pass over 32
  repositories. Wrote `relay/artifacts/T-human-002/review.md` (search log,
  publications table covering 30 methods and resources, repository commentary,
  data sources, failure modes, one-page opinion), `refs.bib` (47 entries, 2
  marked TO VERIFY), and `repos.tsv` (32 rows with last commit, 12-month
  commit counts, issues, stars, licence). Found that Gnomon and EGAPx have no
  methods paper (0 Europe PMC title hits) and that `nekrut/axomeme` and
  `nekrut/scalingPaper` return 404, so the charter's HyphAeon and EGAPx cost
  claims are unverifiable from here; posted a `question` to human about it.
  NEXT: attempt installs in a fresh environment (every `install_tested` cell
  is currently `no` — largest gap), resolve the 4 unresolved repositories and
  2 unverified DOIs, fetch the ANNEVO preprint, confirm EGAPx's hardware floor
  and clade exclusions from repo docs, run bioRxiv/arXiv/OpenAlex passes, then
  move to `review`.
- 2026-09-09 lenin: -> in_progress.
- 2026-09-09 lenin: run 2. Closed the install gap: built miniprot (81f9b93)
  and minisplice (49f9e8c) from source, both clean; installed Tiberius
  (e73844b) and egapx (f9a7392) in fresh venvs. Findings: Tiberius's README
  quick-start `pip install .` installs only the launcher and then fails with
  ModuleNotFoundError (needs the `from_source` extra, 6.3 GB, 1m46s); on an
  RTX 5080 its pinned TensorFlow must JIT from PTX ("could take 30 minutes or
  longer"); and the installed CLI exposes a ClaMSA comparative mode that is
  mammals-only while the six-clade extension is sequence-only — recorded as
  section 6.1, the sharpest open question for T-human-011. Verified EGAPx's
  32 CPU / 256 GB floor, clade exclusions and per-genome runtimes (71 CPU-hrs
  for 144 Mb Drosophila, 425 for 1.1 Gb chicken) from its README, so
  T-human-009 no longer depends on the unreachable scalingPaper repo except
  for the Galaxy failure-rate figures. Resolved MAKER (not on GitHub;
  funannotate is the successor), GlimmerHMM (no canonical repo) and
  Ensembl/ensembl-anno (a third production pipeline, missed in run 1); added
  5 repos (ANNEVO, minisplice, funannotate, OpenSpliceAI, GeneMark-ETP).
  Settled the GeneMark licence question: GeneMark-ETP is on GitHub with no
  LICENSE file, so BRAKER3's redistributability is unsettled. Verified both
  TO VERIFY DOIs via Crossref. refs.bib now 50 entries with no unverified
  DOIs; repos.tsv 37 rows; review.md 4,414 words excluding tables.
  NEXT (run 3, should be the last): bioRxiv/arXiv/OpenAlex passes and a
  GitHub code-search snowball, fetch the ANNEVO preprint, backfill commit and
  issue counts for the 5 new repos once the GitHub rate limit resets, then
  move to `review`.
- 2026-09-09 lenin: -> in_progress.
- 2026-09-09 lenin: run 3, final. Ran the passes run 2 flagged as the last
  substantial gap: OpenAlex (8 queries, 163 unique works) and a Europe PMC
  preprint pass (`SRC:PPR`, 6 queries), plus 17 DOI metadata lookups, 13
  authenticated GitHub `/repos` calls and 7 repository searches. The preprint
  pass was not a formality: it produced a new section 2.1 with **eleven
  2025-2026 methods runs 1 and 2 missed entirely**, five posted in the last
  five months -- OrionGeno, Vipsania, TOGA2, GeneCAD, GENATATORs, the He &
  Florea foundation-model benchmark, SegmentNT, geneML, GeMoSeq, OMAnnotator,
  and ANNEVO now read first-hand. Two of them change the review's conclusion,
  so section 6.2 revises the run-1 opinion: OrionGeno is already a
  phylogeny-aware deep annotator run on >5,300 genomes, so comparative
  geometry is occupied ground; what is unoccupied is the ablation isolating
  which geometry does the work, the ~2M-parameter budget, and a neutral
  benchmark. Recommends T-human-011 become a three-arm comparison
  (sequence-only / tree-as-metric / unsupervised) rather than one
  architecture. Added sections 5.9 (reference annotations are measurably
  wrong: 11% and 5% of vertebrate selenoproteins well annotated in
  Ensembl/GenBank; 24% of one nematode annotation corrected by curation) and
  5.10 (preprint-stage tools are not reliable artifacts). Added section 3.2
  with the run-3 repositories and a licensing summary -- OrionGeno is
  non-commercial and GeneMark-ETP and minisplice have no LICENSE, which
  constrains what T-human-007 can redistribute. Backfilled the five run-2
  repository rows the rate limit had left at `?`. refs.bib 64 entries, all
  with DOIs, parse-checked and every cited key resolved; repos.tsv 41 rows;
  review.md 5,805 words excluding tables. NOT done and stated in section 7:
  arXiv returns HTTP 429 to this host on every attempt (covered indirectly
  via OpenAlex `type:preprint`); no install attempted for the run-3 cohort;
  no code-search/dependents snowball. Moving to `review`.
- 2026-09-09 lenin: -> review.
- 2026-09-09 lenin: lease renewed; artifact unchanged, still in `review` awaiting the coordinator.
