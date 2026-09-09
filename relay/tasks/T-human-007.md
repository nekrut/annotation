---
id: T-human-007
title: Design the benchmark: species panel, references, metrics, leakage rules
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T07:15:00Z
depends_on: []
touches: [docs/benchmark.md, benchmark/]
pr: https://github.com/nekrut/annotation/pull/5
---

## Goal

A benchmark that a model can be judged on and that rewards cross-clade
generalization rather than memorization.

## What to produce

- `docs/benchmark.md`: a species panel of 12 to 20 genomes covering
  mammals, birds, fish, insects, nematodes, plants, fungi, and at least one
  protist, chosen to span intron length distributions (median and tail),
  genome size, GC content, and gene density. Include the numbers that
  justify each pick. State which reference annotation is truth for each
  (RefSeq, Ensembl, MANE, community curated) and its known weaknesses.
- Metrics: nucleotide, exon, and gene level sensitivity and specificity;
  splice-site accuracy; BUSCO completeness of predicted proteomes; wall
  clock and peak memory per megabase. Define them precisely.
- Held-out rules: which species are never used in training, how to detect
  leakage through alignments that include held-out species, and a minimum
  phylogenetic distance rule.
- `benchmark/`: a manifest (`panel.tsv`) with assembly accessions, annotation
  versions, and checksums, plus a fetch script that does not commit data.

## Definition of done

Pull request merged after one `review`. The fetch script runs for at least
one species in a fresh clone.

## Log

- 2026-09-09 human: created.
- 2026-09-09 lenin: claimed.
- 2026-09-09 lenin: run 1. Built the panel. 20 genomes across 8 clades, chosen
  by NCBI Datasets metadata and then justified with statistics recomputed
  locally rather than quoted: `benchmark/annotation_stats.py` streams each
  reference GFF3 and derives isoform-de-duplicated intron and exon
  distributions. Spread achieved: median intron 56 bp (S. pombe) to 1,781 bp
  (human), 32x; intron p99 407 bp to 107,644 bp, 264x; GC 19.5% (P.
  falciparum) to 48.5% (N. crassa); gene density 6.5 to 499 per Mb; coding
  fraction 3.6% (maize) to 73.1% (S. cerevisiae). Deliberately oversamples the
  tail EGAPx excludes (3 fungi, 3 protists, 1 nematode) and includes
  Tetrahymena for nuclear genetic code 6, which makes every hard-coded stop
  prior wrong and turns that into a diagnostic rather than a nuisance.
  Three splits: train (10), heldout (6), heldout_paired (4). The minimum
  distance rule (MRCA of a held-out species with the nearest training species
  at or above CLASS) is enforced by `benchmark/leakage_check.py` against a
  committed NCBI Taxonomy cache, so it is checked and not asserted; the four
  paired species break it on purpose to isolate a regime shift from a clade
  shift. docs/benchmark.md defines four leakage channels (label, sequence,
  alignment, evidence), a required per-run declaration block, and the metrics
  to implementable precision. Verified in a fresh clone: fetch downloads and
  double-checksums S. cerevisiae, stats and leakage check both run with
  nothing installed. PR #5.
  NOT done, listed as section 7 of the doc: no scorer (`benchmark/score.py`);
  no high-confidence subset to bound annotation error; the distance rule is
  taxonomic rank rather than substitutions per site; no RNA-seq accessions
  pinned; BUSCO/OMArk lineage datasets not pinned; `--what fasta` untested at
  3 Gb scale.
  NEXT: address review feedback, then write the scorer.
- 2026-09-09 lenin: -> review (https://github.com/nekrut/annotation/pull/5).
- 2026-09-09 lenin: run 2. Wrote the scorer, which run 1 listed as open item 1.
  `benchmark/score.py` (stdlib, 1 file) implements section 4.1 nucleotide with
  MCC, 4.2 exon with the four-way type stratification and the
  overlap-but-no-boundary class, 4.3 donors and acceptors stratified by
  intron-length decile and, with `--genome`, by dinucleotide class and local
  GC, 4.4 transcript exact match under the isoform rule plus locus matching
  with fusion and split counts, and 4.5 start/stop codons. It refuses to run
  without the section 3.3 declaration and records its SHA-256.
  `benchmark/report.py` assembles the 4.8 table and the two aggregates and
  refuses to average an incomplete panel without `--partial`.
  Real data contradicted two definitions in the doc, both now amended:
  intron-length deciles cannot come from `panel.tsv` (it carries only
  p10/median/p90/p99) so the scorer computes them from the reference; and a
  CDS gap under 20 bp is not a splice junction, because 47 of the 343 CDS gaps
  in the S. cerevisiae reference are Ty programmed-frameshift 1 bp gaps, which
  would have put a 14% floor of non-splice-sites into the yeast donor counts.
  A third came from the locus metric: a fusion now requires the two annotated
  genes not to overlap each other, since yeast has 91 same-strand
  CDS-overlapping gene pairs and without that condition a *perfect*
  prediction scores 137 fusions and 132 splits.
  Verified: `--self-test` checks 30 counts over three fixture pairs plus a
  self-comparison that must be exactly 1.0 everywhere and a check of the
  streaming FASTA window reader; on the real yeast reference (17 seqs, 12.16
  Mb, 6,027 transcripts) it runs in 1.1 s, scores 1.0/MCC 1.0 against itself
  with 0 fusions and 0 splits, and on a degraded copy (10% transcripts
  dropped, 10% of CDS 3' ends shifted 3 bp) returns nucleotide F1 0.947, exon
  0.856, donor 0.892, transcript 0.852, locus 0.949, recovering 270 GT-AG, 8
  GC-AG and 18 other donors from the FASTA. `fetch.py --what fasta` also
  exercised for the first time (3.8 MB, checksum verified), closing part of
  run 1's open item 7.
  NOT done, now open items 1 and 2 of the doc: no BUSCO/OMArk or cost columns
  (external; T-human-009 owns the cost half), and the scorer has been run on
  one species only -- the primary-assembly filter for alt loci and patches is
  a RefSeq convention untested on the human GFF3 that actually has them, and
  untested on Ensembl input.
  Task stays in `review`; PR #5 updated. NEXT: run the scorer on human and one
  Ensembl-annotated species to test the primary-assembly filter, and address
  review feedback when it arrives.
- 2026-09-09 lenin: run 3. Closed run 2's open item 2: ran `benchmark/score.py`
  on real references beyond yeast for the first time -- human, Arabidopsis,
  C. elegans, maize, Nematostella, Plasmodium, Tetrahymena, and the Ensembl
  (not RefSeq) C. elegans annotation. It found three defects, all now fixed and
  all invisible on yeast.
  (1) Organelles were being scored. `genome=mitochondrion` and
  `genome=chloroplast` regions passed the filter, so the human, yeast,
  Arabidopsis and maize mitochondria and the Arabidopsis and maize chloroplasts
  were scored against the nuclear genetic code; human chrM is translation
  table 2, so every start/stop check on it was wrong, and organelles are out of
  charter scope anyway.
  (2) The alt-locus rule discarded every unplaced scaffold in the panel. It
  keyed on `genome=genomic` plus `chromosome=`, but all 680 human, all 675
  maize and all 32 Nematostella genomic regions carry a `chromosome=` (unplaced
  ones say `chromosome=Unknown`). The real discriminator is `map=`: a
  cytogenetic band means alt locus or patch, `map=unlocalized` or no `map=`
  means scaffold. Human goes from 25 scored sequences to 102, +13.3 Mb and
  +179 transcripts.
  (3) A *perfect* human prediction scored 82 fusions and 82 splits. Run 2's
  pairwise "do the two annotated genes overlap each other" guard does not
  survive overlap chains: 82 human loci have a gene B overlapping both A and C
  while A and C are disjoint. Fusions are now counted across connected
  components of the reference's own same-strand CDS-overlap graph, which makes
  an identity run score exactly 0 fusions and 0 splits on all eight species and
  subsumes the yeast fix.
  Added, because two of the three were silent: every result now carries
  `sequence_selection` (kept, and dropped grouped by reason) and
  `predicted_transcripts_not_scored` /
  `predicted_sequences_absent_from_reference`, with a stderr warning when more
  than half the prediction is off-panel -- scoring the Ensembl C. elegans
  annotation against the RefSeq one returns F1 0.0 everywhere purely because
  Ensembl calls chromosome I `I` and RefSeq calls it `NC_003279.8`.
  Verified: `--self-test` now checks 35 counts including two sequence-selection
  fixtures (RefSeq-shaped and Ensembl-shaped); eight identity runs at F1 1.0 /
  MCC 1.0 with 0 fusions and 0 splits; human 37 s and 0.93 GB peak on one
  laptop core over 3.10 Gb and 132,030 transcripts; degraded human (10%
  transcripts dropped, 10% of CDS 3' ends shifted 3 bp) returns nucleotide F1
  0.984, exon 0.821, donor 0.891, transcript 0.366, locus 0.986 with 0 fusions
  and 7 splits; `report.py` unchanged and still assembles the table.
  docs/benchmark.md sections 4, 4.4, 6 and 7 rewritten to match.
  Task stays in `review`; PR #5 updated.
  NOT done: 12 panel species still unscored; `--genome` never exercised at
  vertebrate scale, so the section 4.3 dinucleotide and GC stratifications are
  yeast-only; the no-`region`-features path is a name heuristic, now open
  item 5. NEXT: address review feedback on PR #5 when it arrives.
