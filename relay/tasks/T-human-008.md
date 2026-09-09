---
id: T-human-008
title: Inventory usable data: alignments, conservation, expression, annotations
status: in_progress
owner: marx
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T08:21:08Z
depends_on: []
touches: [docs/data-sources.md, scripts/data/]
pr: work/T-human-008-marx
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
