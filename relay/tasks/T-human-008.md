---
id: T-human-008
title: Inventory usable data: alignments, conservation, expression, annotations
status: review
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T13:52:35Z
depends_on: []
touches: [docs/data-sources.md, scripts/data/]
pr: https://github.com/nekrut/annotation/pull/8
---

## Goal

Know exactly what data the new model can be trained on, how to get it, and
how to cut it into training examples.

## What to produce

- `docs/data-sources.md`: for UCSC (multiz and Cactus whole-genome
  alignments and their trees, phyloP and phastCons, RNA-seq, regulation,
  variation tracks), Ensembl Compara, NCBI RefSeq, Zoonomia, and any others
  you find: what assemblies are covered, formats (MAF, HAL, bigWig, GFF3),
  sizes, access methods (rsync, bigDataUrl, APIs), licences and rate limits.
- A section on non-coding alignment: how well introns, UTRs, and
  intergenic sequence are aligned in each source and at what evolutionary
  distances alignment breaks down.
- `scripts/data/fetch_window.py`: given an assembly, a locus, and a
  flanking size, fetch the multiple alignment block with its Newick tree
  and the reference annotation for that window, in a format a model can
  consume. Demonstrate on one human locus and one non-mammal locus.
- A short assessment of which sources are suitable for training versus
  held-out evaluation given the benchmark's leakage rules.

## Definition of done

Pull request merged after one `review`. The window script runs in a fresh
clone for both demonstrated loci without committing data over 5 MB.

## Log

- 2026-09-09 human: created.
- 2026-09-09 trotsky: claimed.
- 2026-09-09 marx: claimed (previous lease by trotsky expired).
- 2026-09-09 marx (tick 1 on this task): claimed after trotsky's lease expired with no logged work or branch (note 20260909T055335Z-marx-0005). Probed UCSC API/hgdownload, Ensembl REST/FTP, NCBI, Zoonomia and GenArk from the runner; wrote `docs/data-sources.md` (inventory with measured sizes, access, licences, rate limits; measured non-coding alignability on three windows; training-versus-held-out assessment under the T-007 leakage rules: 12 of 20 panel species have no public multiple alignment, 7 of the 10 held-out ones) and `scripts/data/fetch_window.py` (stdlib; UCSC bigMaf via API, wigMaf via HTTP Range reads on the uncompressed /gbdb MAF, Ensembl EPO with ancestral rows; MAF + Newick + FASTA + annotation + conservation + manifest with per-species coverage by CDS/UTR/intron/intergenic). Demonstrated on human HBB (470-way and Cactus 241-way), fly Adh (124-way), chicken GAPDH (Ensembl sauropsids), mouse Gapdh (Ensembl mammals, mouse as query) and a GenArk-only Plasmodium window; table with checksums in `scripts/data/README.md`. Branch `work/T-human-008-marx` pushed; no PR yet. Next: aggregate the coverage measurement over the benchmark gene sample, draft the window-cutting convention for training examples, re-read the Ensembl licence text, then open the PR and move to review with a note.
- 2026-09-09 marx: -> in_progress (work/T-human-008-marx).
- 2026-09-09 marx: -> in_progress.
- 2026-09-09 marx: -> review (https://github.com/nekrut/annotation/pull/8).
- 2026-09-09 marx (tick 2 on this task): closed open items 1, 2 and 5 of section 9. Aggregated non-coding alignability over a fixed-seed gene sample (12 fly windows on dm6 124-way, 10 human on hg38 100-way, 122.6 MB in 258 requests) binned by patristic distance with new `scripts/data/coverage_by_distance.py`: introns/UTR/intergenic are mostly aligned below 0.5 substitutions per site and essentially unaligned past 1 in both clades, CDS persists (0.73 fly, 0.48 human beyond 2); 110 of 123 fly and 34 of 99 human informants are past that horizon. Wrote the window-cutting convention as `scripts/data/cut_windows.py` (stdlib `.npz`, gap/unaligned states, insertion channel, frame, boundaries, `--drop-species`, reverse complement, 12-check self-test) and section 6.3. Re-read and quoted the Ensembl licence. `hgdownload.soe.ucsc.edu` reset every connection for several minutes mid-run; `fetch_window.py` now rotates mirrors and takes `--download-host`. Opened PR #8, moved to review. Next: answer review comments; remaining section 9 items (mouse/worm aggregation, 470-way clade counts) are optional follow-ups.
- 2026-09-09 marx: -> review (https://github.com/nekrut/annotation/pull/8).
- 2026-09-09 marx (tick 3 on this task): lease renewed; PR #8 has no reviews, comments or CI yet, so no changes to the branch this tick. Spent the tick reviewing T-human-007 (PR #5), which with this task gates T-human-010 and T-human-011: posted `review` 20260909T075645Z-marx-0007 (accept; pseudogene CDS scored as truth, greedy-versus-optimal wording, and a proposed split of the section 3.2 alignment rule into reference-anchored versus jointly inferred alignments, citing `docs/data-sources.md` section 7). If lenin adopts the split, section 7 here and the `--drop-species` note in 6.3 need one sentence each on `alignment_rows_dropped` in the declaration. Next: answer any review of PR #8; optional section 9 items unchanged.
- 2026-09-09 marx: -> review (https://github.com/nekrut/annotation/pull/8).
- 2026-09-09 marx (tick 4 on this task): lease renewed; PR #8 still has no reviews or comments. lenin's note 20260909T082157Z-lenin-0010 adopted the reference-anchored versus jointly inferred split into `docs/benchmark.md` 3.2 and added the required `alignment_rows_dropped` key, so the follow-up flagged last tick is done: pushed one commit to `work/T-human-008-marx` (PR #8) that restates the split in section 7 as the benchmark now has it, names the rows to drop for the two affected panel alignments (`apiMel4` in dm6 124-way, `hg38` in mm39 35-way, the latter checked against `mm39.35way.nh`), marks the jointly inferred fish and mammal EPO sets as rebuild-only, says in 6.3 and the `cut_windows.py` docstring that the sidecar's `dropped_species` list is the value of that key, and records lenin's answer to section 9 item 3 (with/without-informants column left to T-human-011). Self-test still 12/12. Next: answer any review of PR #8; optional section 9 items unchanged.
- 2026-09-09 marx: -> review (https://github.com/nekrut/annotation/pull/8).
- 2026-09-09 marx (tick 5 on this task): lease renewed; answered lenin's review 20260909T090919Z-lenin-0011 on PR #8 with two commits on `work/T-human-008-marx` (520a92b, 0e554d1). All four items done: alignment class from an explicit track table with a strict default and `--reference-anchored` override (sidecar records which decided); EPO ancestral rows read from the manifest's `ancestral_sequences`, tied to clades through the block trees' internal node names, and dropped with any `--drop-species` member of their clade (`dropped_ancestors`); `label_counts` counted through the mask so padded `.rc` sidecars agree; EPO row order from the species set (`fetch_window.py` now records `species_set_members`), ancestors appended sorted, `--drop-ancestors`, and 6.3 rewritten to say which part of `K` is fixed. Also from the review's smaller list: `--transcript-types` defaulting to the benchmark truth rule, a minus-strand reference block counter, and the `.rc` last-base insertion documented. Found and fixed on the way: Ensembl tree leaves (`species_region_start_end[strand]`) never matched the distance lookup, so `dist` was NaN for every EPO informant. Self-test 12 -> 27 checks including lenin's fixtures. Verified live on the chicken GAPDH sauropsid window (14 requests, 93 KB; Ensembl REST returned transient 500s twice before succeeding): K=16, 5 ancestors dropped with zebra finch, forward and rc label counts identical, real distances. Posted note replying to the review. Task stays in review awaiting a second look. Next: answer further review; optional section 9 items unchanged.
- 2026-09-09 marx: -> review (https://github.com/nekrut/annotation/pull/8).
- 2026-09-09 marx (tick 6 on this task): lease renewed; PR #8 still has no GitHub reviews or comments. Acted on lenin's handback in note 20260909T101258Z-lenin-0012 (T-human-007 keeps zebrafish on GRCz12ab; both Ensembl fish EPO sets contain the held-out fugu, so the fish alignment is rebuild-only on either assembly): pushed e2bf1aa to `work/T-human-008-marx` counting zebrafish with the species that lack a usable training alignment (thirteen of twenty) in section 1, the section 7 table row and first consequence, section 8 (with GRCz12ab as the assembly to build a fugu-free fish clade on), and section 9 item 6 marked answered. PR description updated to match. Self-test 27/27. Nothing else in the review thread is open; the branch is ready for the coordinator to merge. Next: answer further review if any; optional section 9 items 1, 4 and 5 unchanged.
- 2026-09-09 marx: -> review (https://github.com/nekrut/annotation/pull/8).
- 2026-09-09 marx (tick 7 on this task): lease renewed; PR #8 still has no GitHub reviews or comments. Closed section 9 item 5 (the optional 470-way clade count) with one commit 219243c on `work/T-human-008-marx`: new stdlib `scripts/data/tree_composition.py` (Newick parser, patristic distance from the reference, per-clade leaves/species/min/median/max and counts within 0.5 and beyond 1.0 subst/site; validated by reproducing the 6.1 human 100-way histogram exactly) and `scripts/data/hg38.470way.orders.tsv` (470 leaves to order/family via GBIF, 35 KB). Findings folded into the document: the 470-way is mammal-only with 431 distinct species for 470 leaves (35 duplicated as assemblies or subspecies, about 8% of rows); 413 of 469 informants lie within 0.5 of human and only platypus and echidna beyond 1, so section 6.2 item 2 had it backwards (it predicted the 470-way would be sparser than the 100-way) and now says the informant axis is dense and near-redundant on the mammal tree, sparse on vertebrate-wide and fly trees. Also corrected section 2.1: the 470-way does have raw per-chromosome MAF under `goldenPath/.../maf/` (499 files, 5.9 TB, chr1 471 GB, no gzipped copy), so rsync is not a laptop path for it (section 9 item 7 amended). PR description updated. Self-test 27/27. Next: answer further review if any; optional section 9 items 1 and 4 unchanged.
