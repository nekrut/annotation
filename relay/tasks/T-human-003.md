---
id: T-human-003
title: Independent review of gene prediction literature and software (slot 2 of 4)
status: in_progress
owner: engels
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T03:23:03Z
depends_on: []
touches: [relay/artifacts/T-human-003/]
pr: null
---

## Goal

One independent, careful review of the state of eukaryotic gene prediction:
the published methods and the software that is actually alive. Four
identical tasks exist (T-human-002 to T-human-005), one per agent. Claim any
open one; hold at most one; do not read the other agents' review artifacts
until yours is in `review`.

## What to produce

Write `relay/artifacts/T-human-003/review.md` with these sections, plus
`relay/artifacts/T-human-003/refs.bib` and `relay/artifacts/T-human-003/repos.tsv`.

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
- 2026-09-09 engels: claimed.
- 2026-09-09 engels: -> in_progress.
- 2026-09-09 engels (01:51 UTC): completed the bounded first-run research pass begun at 01:23 UTC. Created the independent [review draft](../artifacts/T-human-003/review.md), bibliography, repository inventory, source/search logs and reproducible inspection helpers. Read core open papers and separated exact-structure metrics, base-level scores, simulated classifier errors, and preprocessing costs. BibTeX parsing and artifact checks passed; Tiberius launcher installation and configuration listing passed in a fresh venv, with inference untested. No other agent review artifacts were read. Next: finish historical-method and supplementary-table extraction, unresolved code/license and installation checks, and the Gnomon/EGAPx citation requirement. Task remains in_progress.
