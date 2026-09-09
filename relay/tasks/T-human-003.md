---
id: T-human-003
title: Independent review of gene prediction literature and software (slot 2 of 4)
status: in_progress
owner: engels
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T10:20:36Z
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
- 2026-09-09 engels: -> in_progress.

- 2026-09-09 engels (02:08 UTC): second bounded tick. Renewed the lease and processed the full inbox; no questions required answers. Expanded classical/comparative results and costs, added the original GeneMark-ES and directly relevant RefSeq/Gnomon/EGAPx citations, and added three repository snapshots with further license inspection. SNAP compiled and both README examples ran in a fresh source directory with the host compiler; reproducible commands and hashes are retained. Artifact verification passed (29 BibTeX entries, 24 DOI-bearing publication rows, 18 repository rows), as recorded in [verification.json](../artifacts/T-human-003/verification.json). No other agent review artifacts were read. Next: remaining historical/full-text and supplementary extraction, modern-model split/metric reconciliation, and remaining software checks. Task remains in_progress.
- 2026-09-09 engels: -> in_progress.

- 2026-09-09 engels (02:24 UTC): third bounded tick. Renewed the lease and processed the full unread inbox; no questions required answers. Recovered original TWINSCAN benchmark/split details, MAKER Table 1 metrics and historical runtime (flagging a prose/table sensitivity-specificity reversal), and PhyloCSF++ speed/agreement results. Citations, PDF checksums, search queries and access failures are retained in [the research addendum](../messages/20260909T022456Z-engels-0003.md). Current operator write limits restrict this tick to a new message and owned task/roster edits; existing review artifacts were not edited. No other agent review artifacts were read. Next: remaining original/full-text and modern-model split/metric checks, and artifact integration when authorized. Task remains in_progress.
- 2026-09-09 engels (03:27 UTC): fourth bounded tick. Renewed the lease and processed all four unread messages; no addressed question required an answer. Audited Helixer's published GffCompare command and terminal-boundary semantics, recovered its matched-training mammal comparison and hardware, and reconciled Tiberius's longest-CDS policy/small-model ablation with SegmentNT's supervised holdout and released pretraining corpus. Findings, limitations, source pins, and checksums are in [the research addendum](../messages/20260909T032715Z-engels-0004.md). No other agent review artifact was read. Operator write limits were respected: only a new message and owned task/roster edits. Next: remaining full-text and scorer/version checks; artifact integration remains pending. Task stays in_progress.
- 2026-09-09 engels (04:26 UTC): fifth bounded tick. Renewed the lease and processed all four unread messages; no addressed question required an answer. Recovered the open ClaMSA preprint and supplement, recording cross-clade AUCs, candidate construction and version differences from the journal abstract. Audited ANNEVO's preprint, published supplement and pinned repository for clade-specific training, architecture changes, scorer boundary tolerance and runtime comparability. Findings, source pins, access limits and checksums are in [the research addendum](../messages/20260909T042636Z-engels-0005.md). No other agent review artifact was read. Only the new message and owned task/roster were edited under the operator's write limits. Next: remaining journal/version and evaluation checks; artifact integration remains pending. Task stays in_progress.

- 2026-09-09 engels (05:25 UTC): sixth bounded tick. Renewed the lease and processed both unread broadcasts; neither required an answer. Audited the final AUGUSTUS-CGP paper for geometric structure, clade-level CPU/memory accounting and exact-versus-overlap exon scoring; checked pinned current documentation for candidate filters, separate intron mechanisms, tree units and reference-trained feature scores. Recovered PhyloCSF's length-calibration improvement and composition-only behavior when informants are absent. Citations, access limits and checksums are in [the research addendum](../messages/20260909T052453Z-engels-0006.md). No other agent review artifact was read. Only a new message and owned task/roster files were changed under the operator's write limits. Next: unresolved original-paper/supplementary and version checks; artifact integration remains pending. Task stays in_progress.

- 2026-09-09 engels (06:26 UTC): seventh bounded tick. Renewed the lease and processed all three unread broadcasts; no question required an answer. Recovered GeneMark-EP+ supplementary gene accuracy across the six evaluated species and checked 24 selected values against pinned repository tables using decimal rounding. Audited zebrafish reference filtering, flagged the paper/repository annotation-release mismatch, and checked GeneMark-ETP reference denominators, unsupported-candidate pruning and GC-dependent training. Citations, checksums and access details are in [the research addendum](../messages/20260909T062545Z-engels-0007.md). No other agent review artifact was read. Only a new message and owned task/roster files were changed under the operator write limits. Next: remaining modern-model/original-paper checks and annotation provenance; artifact integration remains pending. Task stays in_progress.

- 2026-09-09 engels (07:28 UTC): eighth bounded tick. Renewed the lease and processed all three unread broadcasts; no question required an answer. Audited geneML's open preprint, supplements and pinned source: reconciled its 761/752 training manifests, checked 45 accuracy means and five runtime means, recovered strict zero-tolerance scoring, and confirmed the default hard intron ceiling. Also recovered the original Helixer's base-level metrics and combined development/test reporting. Findings, discrepancies, source pins and hashes are in [the research addendum](../messages/20260909T072742Z-engels-0008.md). No other agent review artifact was read. Only a new message and owned task/roster files were changed under the operator write limits. Next: remaining historical-source and checkpoint/input provenance gaps; artifact integration remains pending. Task stays in_progress.

- 2026-09-09 engels: first push raced with lenin and was rejected. Rebased cleanly, read the newly arrived lenin-0009 broadcast (no question to answer), and advanced last_seen to the fetched main commit b29282a719d2e7caa664577c8c802bce4a0dea39. The contract, charter and owned task were unchanged remotely.

- 2026-09-09 engels (08:26 UTC): ninth bounded tick. Renewed the lease and processed the one unread broadcast; no question required an answer. Audited original KA/KS input construction, simulated negatives and an apparent likelihood-ratio sign inconsistency; checked historical HMR195 selection/averaging rules and GENSCAN training overlap, independent test results and GC-conditioned duration priors. Findings, arithmetic checks, source hashes and access limitations are in [the research addendum](../messages/20260909T082602Z-engels-0009.md). No other agent review artifact was read. Only a new message and owned task/roster files were changed under the operator write limits. Next: remaining checkpoint/input provenance and historical-source gaps; artifact integration remains pending. Task stays in_progress.

- 2026-09-09 engels: first push raced with another relay update and was rejected. Rebased cleanly and read the newly arrived lenin-0010 broadcast; no question required an answer. The protocol, charter and owned task were unchanged remotely. Advanced last_seen to fetched main 762f901f17ef2e24cb98044816a0e6717d4c56f4 before retrying.
