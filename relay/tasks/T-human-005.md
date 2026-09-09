---
id: T-human-005
title: Independent review of gene prediction literature and software (slot 4 of 4)
status: review
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T13:52:35Z
depends_on: []
touches: [relay/artifacts/T-human-005/]
pr: null
---

## Goal

One independent, careful review of the state of eukaryotic gene prediction:
the published methods and the software that is actually alive. Five
identical tasks exist (T-human-002 to T-human-005 and T-human-012), one per agent. Claim any
open one; hold at most one; do not read the other agents' review artifacts
until yours is in `review`.

## What to produce

Write `relay/artifacts/T-human-005/review.md` with these sections, plus
`relay/artifacts/T-human-005/refs.bib` and `relay/artifacts/T-human-005/repos.tsv`.

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
- 2026-09-09 marx: claimed.
- 2026-09-09 marx (tick 1): claimed. Ran 87 Europe PMC queries plus abstract retrieval for 54 DOIs (OpenAlex and Semantic Scholar rate-limited; GitHub API blocked from this runner, so repository activity measured with anonymous shallow clones). Wrote `review.md` draft (search log, 38-row publications table, repo summary, outlines for sections 4 to 6), `refs.bib` (54 entries, braces balanced), `repos.tsv` (30 repos with last commit, 12-month commit count, licence, HEAD hash; stars and issues `n/a`). Next: fetch OA full text (Tiberius 2024, BRAKER3, GeneMark-ETP, Helixer 2021, CONTRAST, G3PO, fitild) for exact accuracy and runtime numbers; bioRxiv full text of the four 2026 preprints for code links and species panels; add aligner and evaluator rows; fresh-venv install tests for Tiberius, Helixer, egapx, BRAKER, EviAnn; then write sections 4 to 6.
- 2026-09-09 marx: -> in_progress.
- 2026-09-09 marx (tick 2): Europe PMC full-text XML, OpenAlex, Semantic Scholar and raw.githubusercontent.com became reachable (GitHub API still 403; bioRxiv JATS 429). Read 17 OA full texts (Tiberius 2024, BRAKER3, GeneMark-ETP, Helixer 2021/2026, CONTRAST, G3PO, fitild, SNAP, GeneMark-EP+, BRAKER2, sensor-NN 2023, GeneMark-ES, GALBA, PhyloCSF, EGASP, MAKER2) and replaced abstract-level placeholders in the publications table with exact F1, Sn/Sp, runtime and hardware figures; added rows for ClaMSA, EVM, Mikado, aligners, evaluation tools, scaling studies (44 rows). Filled stars and open issues for all 30 repos via shields.io; fresh-venv install tests: Tiberius installs on Python 3.12 only (package without entry point; needs container image), Vipsania on PyPI for 3.12 only (TF 2.19 + CUDA wheels), Helixer installs on 3.11 but Helixer.py fails on two undeclared dependencies (pyyaml, scikit-learn). Wrote sections 4 to 6 in full. refs.bib now 61 entries, braces balanced. Review is 2,713 words excluding tables. Next: read ClaMSA and AUGUSTUS-cgp full text, retry bioRxiv full text for the four 2026 preprints (species panels, code links), test the Helixer container route if possible, final proofread, then move to review with a note to all.
- 2026-09-09 marx: -> in_progress.
- 2026-09-09 marx (tick 3): bioRxiv HTML/PDF/JATS (Cloudflare 429 on all seven attempts), PeerJ and OUP (403), NCBI efetch (publisher-restricted XML) and Europe PMC (no full text for the five PPR records) all refused the remaining full texts, so ClaMSA and AUGUSTUS-cgp are covered from journal abstracts plus repository documentation (clamsa README, Augustus cgp-mode docs), and the row notes say so. Read the v2 abstracts of the Tiberius multi-clade (2026-07-29) and OrionGeno (2026-08-24) preprints via the bioRxiv API and updated those rows. Found and shallow-cloned xjtu-omics/ANNEVO (v2.3.3, 158 stars, custom non-commercial licence, six lineage models, README runtime table) and Gaius-Augustus/clamsa (dormant since 2024-07); repos.tsv now 32 rows. Proofread sections 3 to 6; review is 2,977 words excluding tables; refs.bib 61 entries, braces balanced. Moved to review with a note to all. Not done: no install test for container-only tools (Tiberius, Helixer images, egapx) because the runner has no container runtime; no full text for the five 2026 preprints.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx (tick 4): lease renewed; still in `review`, awaiting the coordinator. No edits to the artifact.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx (tick 5): lease renewed; still in `review`, awaiting the coordinator. No edits to the artifact.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx (tick 6): lease renewed; still in `review`, awaiting the coordinator. No edits to the artifact.
- 2026-09-09 marx: -> review.
- 2026-09-09 marx (tick 7): lease renewed; still in `review`, awaiting the coordinator. No edits to the artifact.
