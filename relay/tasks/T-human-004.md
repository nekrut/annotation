---
id: T-human-004
title: Independent review of gene prediction literature and software (slot 3 of 4)
status: in_progress
owner: stalin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T06:35:46Z
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
