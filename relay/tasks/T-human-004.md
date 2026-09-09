---
id: T-human-004
title: Independent review of gene prediction literature and software (slot 3 of 4)
status: in_progress
owner: stalin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T14:35:38Z
depends_on: []
touches: [relay/artifacts/T-human-004/]
pr: null
---

## Goal

One independent, careful review of the state of eukaryotic gene prediction:
the published methods and the software that is actually alive. Four
identical tasks exist (T-human-002 to T-human-005), one per agent. Claim any
open one; hold at most one; do not read the other agents' review artifacts
until yours is in `review`.

## What to produce

Write `relay/artifacts/T-human-004/review.md` with these sections, plus
`relay/artifacts/T-human-004/refs.bib` and `relay/artifacts/T-human-004/repos.tsv`.

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
- 2026-09-09 stalin: claimed.
- 2026-09-09 stalin: -> in_progress.
- 2026-09-09 stalin: Completed the first bounded independent research pass
  in `relay/artifacts/T-human-004/`: draft review with 21 method rows,
  26-reference bibliography, search/access records, and 13-repository inventory
  with pinned commits and dated activity counts. BibTeX parsed without warnings.
  SNAP built and ran both README examples in a fresh temporary checkout;
  TWINSCAN's documented build failed on this host, with diagnostics preserved.
  No other review artifacts were read. Next: finish primary-text and supplement
  extraction, audit benchmark splits/metrics and input costs, and extend bounded
  installation checks. Keep in_progress; the review is not ready for submission.

- 2026-09-09 stalin: Renewed lease and processed all three unread messages.
  Checked bioRxiv version metadata and pinned OrionGeno/Vipsania repositories;
  audited Vipsania train/test sequence exposure and recomputed selected clade
  means from published scores. Findings, provenance, and remaining limitations
  are in [the tick addendum](../messages/20260909T023936Z-stalin-0002.md).
  Artifact edits remain excluded by AGENTS.md; retain in_progress.
  Next: full-text/supplement access and endpoint/split verification.

- 2026-09-09T03:41:25Z stalin: Renewed lease and processed four unread broadcasts; no
  addressed question required an answer. Audited PhyloCSF evaluation units,
  checked ClaMSA exported rate matrices and logistic-head arithmetic in memory,
  and traced PhyloCSF++ smoothing parameters to coding-label inputs. Findings,
  a pinned repository inventory addition, and checksums are in
  [this tick's addendum](../messages/20260909T034125Z-stalin-0003.md).
  No other review artifact was read; no installation or accuracy run was claimed.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: unresolved primary-text/supplement evidence and historical method gaps.

- 2026-09-09T04:40:13Z stalin: Renewed lease and processed all three unread broadcasts; no
  addressed question required an answer. Checked SNAP transfer metrics and
  evaluation filters, AUGUSTUS manuscript tables and denominator/version caveats,
  and GeneMark-ES 2005/2008 adaptation and intron-model evidence. Findings and
  source checksums are in [this tick addendum](../messages/20260909T044013Z-stalin-0004.md).
  No other review artifact was read; no model or installation run was claimed.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: N-SCAN/AUGUSTUS-CGP primary evidence, journal reconciliation, and
  remaining input-cost checks.

- 2026-09-09T05:40:06Z stalin: Renewed lease and processed all three unread broadcasts; no
  addressed question required an answer. Audited CONTRAST training/validation,
  timing scope, splice-class coverage, and an inconsistent informant count;
  checked N-SCAN_EST evidence preparation and pinned N-SCAN adaptation/input
  documentation. Findings and source checksums are in
  [this tick addendum](../messages/20260909T053948Z-stalin-0005.md).
  No other review artifact was read; no installation or model run was claimed.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: outstanding evidence-pipeline primary-text/supplement checks and
  comparative-model version reconciliation; original N-SCAN journal tables
  remain inaccessible through the attempted open routes.

- 2026-09-09T06:41:53Z stalin: Renewed lease and processed all three unread broadcasts; no
  addressed question required an answer. Recovered AUGUSTUS-CGP journal
  metrics and clade-level runtime scope, audited supplied annotation hints
  and pinned candidate/training documentation, and reconciled all 72 MAKER2
  accuracy cells with supplementary Sn/Sp pairs within rounding. Findings,
  source checksums and access gaps are in
  [this tick addendum](../messages/20260909T064135Z-stalin-0006.md).
  No other review artifact was read; no installation or predictor run was claimed.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: CGP supplement/training provenance and remaining Gnomon/EGAPx or
  modern-model primary evidence.

- 2026-09-09T07:40:36Z stalin: Renewed lease and processed all three unread broadcasts; no
  addressed question required an answer. Audited EGAPx conditional HMM retraining,
  taxonomy/genome-size intron limits, documented costs and EGAP agreement endpoints,
  and RefSeq evidence/reference provenance. Traced active training and parameter
  call paths and checked three intron-limit examples by arithmetic. Findings,
  source checksums and limitations are in
  [this tick addendum](../messages/20260909T074015Z-stalin-0007.md).
  No other review artifact was read; no installation, training or predictor run
  was claimed. Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: unresolved CGP supplement/training provenance and modern-model checkpoint
  and input audits; artifact integration remains pending.

- 2026-09-09T08:40:30Z stalin: Renewed lease and processed both unread broadcasts; no
  addressed question required an answer. Parsed all current Tiberius checkpoint
  configurations, identified folded training-species lists and a mammalian
  split-metadata discrepancy, and recovered the original supplement's ClaMSA
  label provenance and preprocessing costs. Findings, checksums and limitations
  are in [this tick addendum](../messages/20260909T084030Z-stalin-0008.md).
  No other review artifact was read; no installation, training or prediction ran.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: remaining primary-text gaps and checkpoint provenance reconciliation;
  artifact integration remains pending.

- 2026-09-09T09:39:36Z stalin: Renewed lease and processed both unread broadcasts; neither
  required an answer. Audited geneML's enforced intron limits, input-dependent
  score filtering and shuffled-chunk validation; reproduced the numeric score
  override's type mismatch with isolated scalar functions. Findings, checksums,
  full-text access gaps and the unavailable NumPy check are in
  [this tick addendum](../messages/20260909T093936Z-stalin-0009.md).
  No other review artifact was read; no installation, training or prediction ran.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: remaining modern-model endpoint and checkpoint-provenance gaps;
  artifact integration remains pending.

- 2026-09-09T10:40:21Z stalin: Renewed lease and processed all three unread broadcasts; none
  required an answer. Audited Helixer's native structure-matching rules against
  GffCompare 0.12.8 source, parsed checkpoint-specific species membership,
  checked automatic model selection in isolation, and recovered runtime hardware
  and stage scope. Findings, checksums and limitations are in
  [this tick addendum](../messages/20260909T104021Z-stalin-0010.md).
  No other review artifact was read; no installation, training or prediction ran.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: remaining modern-model primary-text/provenance gaps; artifact integration
  remains pending.

- 2026-09-09T11:39:59Z stalin: Renewed lease and processed both unread broadcasts; neither
  required an answer. Audited GENATATOR's crop-conditioned recovery metric,
  default CDS heuristic and stage-model pinning, parsed the dataset species
  table, and reproduced splice-filter strand asymmetry on synthetic motifs.
  Findings, checksums and limitations are in
  [this tick addendum](../messages/20260909T113959Z-stalin-0011.md).
  No other review artifact was read; no installation, training or prediction ran.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: remaining primary-text and checkpoint-provenance gaps; artifact
  integration remains pending.

- 2026-09-09T12:39:11Z stalin: Renewed lease and processed both unread broadcasts; neither
  required an answer. Audited ANNEVO preprint versus current evaluation rules,
  training-label exclusions and active checkpoint selection; verified the decoder's
  short/single-exon score gate and an annotation skip with seven isolated outcomes.
  Findings, checksums and access limits are in
  [this tick addendum](../messages/20260909T123911Z-stalin-0012.md).
  No other review artifact was read; no installation, training or prediction ran.
  Keep in_progress; artifact edits remain excluded by AGENTS.md.
  Next: journal/supplement and checkpoint-membership gaps; artifact integration
  remains pending.
