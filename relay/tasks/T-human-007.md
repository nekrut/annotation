---
id: T-human-007
title: Design the benchmark: species panel, references, metrics, leakage rules
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T15:06:03Z
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
- 2026-09-09 lenin: run 4. Closed run 3's two open items and found one defect by
  doing so. (a) Identity runs now cover all 20 panel references, not 8: every
  species scores F1 1.0 and MCC 1.0 on every metric against itself with 0
  fusions and 0 splits, 1.0 s / 43 MB (P. falciparum, 23 Mb) to 52 s / 335 MB
  (Z. mays, 2.18 Gb), human 38 s / 0.92 GB over 3.10 Gb; all 20 take the RefSeq
  `region` path. (b) `--genome` ran for the first time beyond yeast and at
  vertebrate scale: T. rubripes (384 Mb, 230,048 introns, 32 s / 0.69 GB),
  D. melanogaster, C. elegans, S. pombe, P. falciparum, S. cerevisiae, so the
  4.3 dinucleotide and local-GC strata are now exercised across four kingdoms.
  Non-GT-AG donor fractions recovered from the genomes: 1.57 / 0.98 / 0.97 /
  0.16 / 0.18 / 3.91%.
  DEFECT (b) found: `WindowFetcher` silently dropped windows near the end of a
  sequence and everything queued behind them. It pops a window when its *end*
  has been read, but the queue is sorted by window *start*, so a local-GC
  window centred on a splice site within GC_WINDOW/2 (100 bp) of a sequence end
  can never satisfy that test; it stalls at the head and hides every later
  window on that sequence. The fetcher returns None, which `dinuc_class`
  reports as an `unknown` class and `gc_bin` drops from the denominator -- a
  silent loss, not an error. One T. rubripes intron ending 52 bp from the end
  of a 43 kb scaffold was counted `unknown` rather than the GT-AG it is. Fixed
  by flushing the remaining queue against what was read at each record end,
  clipping the window to the sequence; the buffer is never trimmed past the
  head window's start, so nothing is lost by then. Impact was 1 intron in
  230,048 here and 0 on the other five, because a stall only reaches windows
  behind it, but on a fragmented assembly or with a larger GC_WINDOW it would
  take out the tail of every affected scaffold.
  Verified: regression fixture added (two-record FASTA, overrun on the first
  record and on the last) which fails on the old code with four wrong values;
  self-test now 41 checks, all pass; all six --genome runs repeated after the
  fix, `unknown` count 0 everywhere and the T. rubripes GT-AG total up by
  exactly 1. `fetch.py --what fasta` now exercised on 6 of 20 genomes, largest
  391 Mb.
  Also documented, as open item 10: the local-GC stratification degenerates on
  AT-rich genomes by design -- fixed absolute GC bands keep the column
  comparable across species, but 98.9% of P. falciparum donors land in one bin,
  so that column is to be read across species, not within one.
  docs/benchmark.md sections 6 and 7 rewritten to match. Task stays in `review`;
  PR #5 updated.
  NOT done, now open item 2: degraded-copy runs still cover only 2 of 20
  species, and no real predictor output has ever been scored -- an identity run
  exercises every code path but pins only the fixed points. NEXT: score a real
  AUGUSTUS or Helixer GFF3 on one species to exercise the section 4 conventions
  (stop codon in/out of CDS, Parent shapes, no `region` features) that no
  RefSeq-vs-RefSeq run can reach; address review feedback on PR #5 when it
  arrives.
- 2026-09-09 lenin: run 5. Closed run 4's open item 2 in the half that mattered:
  scored real predictor output for the first time. Installed AUGUSTUS 3.5.0
  (bioconda `augustus-3.5.0-pl5321h5653ebf_10`, no root) and ran it ab initio
  on S. cerevisiae (47 s wall on 17 cores, 237 MB peak) and on S. pombe twice,
  once with its own parameter set and once with the S. cerevisiae one. Its
  GFF3 is the shape no RefSeq-vs-RefSeq run can reach: `transcript` not `mRNA`,
  no `exon`, no `region`, no UTR, stop codon outside the CDS. Three defects,
  all fixed, all invisible on a reference.
  (1) `--stop-outside-cds` WAS A NO-OP. It was recorded in the result and never
  applied, so the 3 bp were never restored. On S. cerevisiae that is exon F1
  0.027, terminal-exon and single-exon F1 0.000, stop-codon F1 0.000 -- while
  nucleotide F1 0.959 and locus F1 0.919, i.e. the metrics a 3 bp 3' shift does
  not move stayed healthy and the run reads as a real result. The flag now
  unions `stop_codon` features into the CDS chain (right for a stop split
  across an intron, and right for a 3'-partial gene, which has none and must
  not be extended), falls back to a counted blind 3 bp extension when there are
  no such features, and -- the part that matters -- the convention is now
  DETECTED from the prediction and reported as `stop_codon_convention_detected`
  with a loud stderr warning when it disagrees with the flag. Corrected, the
  same prediction scores exon F1 0.757 and stop-codon F1 0.912.
  (2) A transcript id reused across sequences was silently welded into one
  chain. AUGUSTUS restarts gene numbering at `g1` in every invocation, so the
  obvious per-chromosome parallel run concatenates to one `g1.t1` per
  chromosome; the loader keyed on the id alone and turned 5,154 transcripts
  into 663 cross-chromosome chimaeras scoring nucleotide F1 0.198 without a
  word. Now counted in `predicted_conflicting_transcript_ids` (466 on the naive
  concatenation) and warned about.
  (3) `report.py` dropped a duplicate species silently. Two results for one
  species -- a run and its ablation, which is exactly what the two S. pombe
  runs are -- kept the last and printed "1 of 20 species scored". It now names
  both files and exits 2.
  THE RESULT THAT MATTERS: swapping AUGUSTUS's parameter set between two
  ascomycete yeasts of nearly the same size and GC costs 62% of exon F1
  (0.773 -> 0.294) and 80% of donor F1 (0.853 -> 0.174) while nucleotide F1
  falls only 9% (0.955 -> 0.867). That is the argument for section 4.8
  reporting a vector rather than a headline, and it is a floor for the
  charter's clade-independence claim: beat a two-yeast parameter swap before
  a mammal-to-fungus claim means anything.
  Verified: `--self-test` 47 -> 69 checks (the "41" in the previous doc
  revision undercounted; the same counting method gives 47 for that code), all
  pass, including an AUGUSTUS-shaped fixture that must score exactly 1.0 with
  the flag and lose exactly its terminal/single exons and stop codons without
  it, a stop_codon-stripped variant, and an id-collision fixture. Identity runs
  on S. cerevisiae, S. pombe and human still 1.0 / MCC 1.0 on every metric with
  0 fusions and 0 splits. Declarations and full JSON results committed under
  `benchmark/validation/` (36 KB, no genome or prediction data).
  docs/benchmark.md: section 4.5 now specifies the stop-codon convention and
  its detection, section 6 gains the AUGUSTUS subsection and the ablation
  table, open items 2 and 5 rewritten, new open item 10.
  Task stays in `review`; PR #5 updated.
  NOT done: real predictor output still covers one tool and two species --
  Helixer and Tiberius emit no `stop_codon` features at all, so the blind 3 bp
  extension is fixture-only, and no evidence-based pipeline (UTRs, alternative
  isoforms) has been scored. Degraded-copy runs still 2 of 20. AUGUSTUS's
  yeast donor F1 0.394 (571 predicted introns against 281 scored reference
  ones in a 95.3% single-exon genome) is recorded as open item 10 rather than
  read as settled. NEXT: Helixer or Tiberius on one vertebrate; address review
  feedback on PR #5 when it arrives.
- 2026-09-09 lenin: run 6. Addressed marx's review of PR #5 (message
  20260909T075645Z-marx-0007, verdict accept). All six items are closed; the
  two the review called defects were real and are fixed in code.
  (1) PSEUDOGENE AND GENE-FRAGMENT CDS ROWS WERE SCORED AS TRUTH. A reference
  CDS row is not automatically a protein-coding gene: RefSeq annotates
  pseudogenes with real CDS blocks carrying `pseudo=true` under a
  `gene_biotype=pseudogene` parent, and immunoglobulin/TCR segments
  (`V_segment`, `D_segment`, `J_segment`, `C_region`) with no start or stop
  codon of their own. Both are now dropped from truth and counted by reason in
  a new `reference_transcript_selection`, and dropped from the *prediction* by
  the same test so an identity run is still exactly 1.0;
  `--score-all-transcripts` restores the old behaviour for auditing. The size
  of the correction is species-specific and larger than the review's yeast
  example suggested: C. elegans has 1,958 pseudogene transcripts with CDS rows,
  6.4% of its 30,548, and human 588 on scored sequences (198 pseudogene, 390
  gene fragments) of 132,030 -- while Arabidopsis, maize, Tetrahymena,
  P. falciparum and Nematostella have none at all, because those references
  annotate pseudogenes without CDS. Finding the human gene fragments needed a
  second fix: the loader only followed `mRNA`/`transcript` parents, so a CDS
  under a `V_gene_segment` had no path to the gene row that says what it is.
  (2) PARTIAL CDS ENDS ENTERED THE 4.5 DENOMINATORS. RefSeq marks an
  incomplete CDS end with `start_range=`/`end_range=` on the CDS row, and
  which *biological* end that is depends on the strand. A 5'-partial gene has
  no annotated start codon, so charging a false negative for missing it is
  unavoidable by any predictor. Both ends are now excluded independently. The
  half the review did not ask for but that the metric needs: dropping only the
  reference position moves the charge to the predictor, whose own start
  somewhere inside that gene becomes a false positive, so a predicted position
  inside the span of a reference transcript partial at that end is dropped
  too, unless it coincides with a surviving reference position. Read the CDS
  row and not the mRNA: 788 of S. pombe's 5,166 mRNA rows carry
  `partial=true` (incomplete UTRs) against 6 of its CDS rows.
  (3) `leakage_check.py` now prints every training species tied at the deepest
  MRCA rank, truncated with a count. Chicken ties with X. tropicalis as well
  as mouse and the frog is the closer relative in time; at Eukaryota all ten
  training species tie. Section 3.1 requotes the new output.
  (4) Section 2.2 now states that `cds_fraction_pct`'s denominator is the
  Datasets API's `genome_bp`, not the GFF's `##sequence-region` total, which
  also counts organelles (S. cerevisiae 73.1 against 72.6, the 85,779 bp
  mitochondrion); `annotation_stats.py` emits
  `cds_fraction_of_sequence_region_pct` so both are visible.
  (5) Section 3.2's alignment rule is split, accepting marx's argument:
  jointly inferred alignments (Cactus, Ensembl EPO) must be rebuilt without
  the held-out species; reference-anchored ones (multiz) may have the row
  dropped at cut time, because that is a filter and not a realignment. The
  affected panel alignments are dm6 124-way (contains apiMel4) and mm39
  35-way (contains human). Section 3.3 gains a required
  `alignment_rows_dropped` key and score.py enforces it.
  (6) Section 4.4 now says the matching is greedy, with its tie-break
  (descending shared CDS bases, then reference id, then predicted id), since
  that is what the code does and "optimal one-to-one" would make a second
  implementation of the spec disagree with this one. Typo in 4.6 fixed.
  Verified: self-test 69 -> 86 checks, all pass, including a fixture built
  from the four RefSeq shapes (pseudo CDS, V segment, plus-strand
  `start_range=`, and a minus-strand `start_range=` that is a *3'* end) which
  must score 1.0 with the filter on and lose exactly two loci with
  `--score-all-transcripts`. Nine identity runs re-verified at F1 1.0 / MCC
  1.0 on every metric with 0 fusions and 0 splits (human, C. elegans,
  Arabidopsis, maize, Tetrahymena, P. falciparum, Nematostella and both
  yeasts); the Ensembl-vs-RefSeq naming-mismatch warning still fires. The
  three AUGUSTUS results were regenerated and are committed: every F1 moved by
  at most 0.002, so the run-5 finding stands -- the two-yeast parameter swap
  costs exon F1 0.774 -> 0.296 and donor F1 0.854 -> 0.175 while nucleotide F1
  falls 0.955 -> 0.868.
  Task stays in `review`; PR #5 updated (head 9325342).
  NOT done, now open items 12 and 13: which non-`protein_coding` biotypes to
  exclude is a judgement call that costs 390 human transcripts and nothing
  elsewhere, and the partial-CDS rule has never met a predictor that itself
  emits truncated genes at contig ends -- `predicted_partial_5prime/3prime`
  are reported for when it does. Marx's remaining note (a with/without
  informants column in report.py, since 7 of 10 held-out species have no
  public alignment) is left for the T-human-011 owner as the review suggests.
  NEXT: Helixer or Tiberius on one vertebrate, which is the input shape no
  AUGUSTUS run reaches (no `stop_codon` features at all, so the blind 3 bp
  extension is still fixture-only).
- 2026-09-09 lenin: lease renewed; task stays in `review` (PR #5 unchanged,
  head 9325342). No work-product change this tick. Spent the tick reviewing
  marx's T-human-008 PR #8, which needed a non-owner review before the
  coordinator can merge and which is the document `docs/benchmark.md`
  section 3.2 was written against: review message 20260909T090919Z-lenin-0011,
  approve after four fixes. Three are in `scripts/data/cut_windows.py`, each
  reproduced with a runnable fixture: (1) `dropped_rows_only` decides
  reference-anchored versus jointly inferred by substring-matching "cactus"
  in the track name, so a jointly-inferred UCSC track not named that (the
  document itself names `hprc90way`) falls to the permissive branch of
  section 3.2 and passes `benchmark/score.py` with row-dropping alone;
  (2) the ancestral-row test is `startswith("ancestor")` but Ensembl EPO
  rows are `Ggal-Mgal[2]`, so `sidecar["ancestral"]` is always empty and
  `--drop-species` cannot reach an ancestor inferred with a held-out
  species; (3) `label_counts` slices `[:b-a]` on reverse-complement
  examples, where the padding has moved to the front, so counts are wrong
  for any window shorter than `L` (verified 3/9 intergenic/CDS forward
  against 11/5 reverse on a 20 bp window at L=30; the tensors are correct,
  only the sidecar lies). The fourth is documentation: 6.3's "K is fixed per
  track" holds for multiz but not for EPO, where per-block trees send row
  order to first-appearance and two windows are not stackable.
  Also flagged that `paint_labels` ignores the transcript type the fetcher
  already carries, so pseudogenes paint labels at loci section 4 removes
  from the truth set.
  ACCEPTED for me to resolve: marx's section 9 item 6, whether the panel
  should keep zebrafish on GRCz12ab (current RefSeq annotation, no
  alignment anywhere) or carry GRCz11 (Ensembl fish EPO) with a declared
  liftover. NEXT: decide that, and Helixer or Tiberius on one vertebrate,
  still the input shape no AUGUSTUS run reaches.
- 2026-09-09 lenin: run 7. Answered the one item run 6 accepted: marx's
  `docs/data-sources.md` §7 item 6, GRCz11 versus GRCz12ab for zebrafish.
  DECIDED: keep GRCz12ab. New `docs/benchmark.md` §2.4 carries the reasoning,
  §2.3's caveat row and `panel.tsv` `notes` point at it.
  The argument, all checked today rather than assumed. (1) Moving to GRCz11
  buys exactly one thing, the Ensembl fish EPO, and that alignment is
  unusable here: both fish species sets, EPO (32 species) and EPO-extended
  (65), contain `takifugu_rubripes`
  (`rest.ensembl.org/info/compara/species_sets/{EPO,EPO_EXTENDED}`), and fugu
  is the `heldout_paired` species paired with zebrafish. EPO is jointly
  inferred, so §3.2 channel 3 already requires a rebuild with fugu removed;
  a dropped row does not qualify. Zebrafish is rebuild-only for the
  comparative arm on either assembly, so the move buys nothing. The
  contamination also sits exactly on the axis the pair exists to test, a
  6.8x median-intron shift inside one subclass.
  (2) Nothing else is on either assembly. UCSC `danRer11` has no
  multiz/phastCons/phyloP at all -- 40 tracks from `list/tracks`, none
  conservation, and no `*way` directory on hgdownload, only vsHg38/vsMm39/
  vsMm10/vsGalGal6 chains. GRCz12ab has no GenArk hub as of today.
  (3) Since a rebuild is needed either way, rebuild on the better assembly:
  GRCz12ab is 25 contigs, N50 59.43 Mb, zero gap bases; GRCz11 is 1,917
  scaffolds, 19,725 contigs, N50 1.42 Mb, 4,689,282 gap bases (NCBI Datasets
  `dataset_report`, both accessions). Zebrafish has the panel's third-longest
  introns (p99 62,699 bp, max 1,090,140 bp), and GRCz11 puts 4.7 Mb of N
  through them. GRCz11's annotation is also two years older
  (RS_2024_08, 27,158 coding genes, against RS_2026_07 and 28,415).
  (4) The liftover option is real and is rejected in §2.4 with the chain it
  would have used: UCSC ships `danRer11ToGCA_052040795.1.over.chain.gz`
  (9.4 MB, 2026-05-08) and `GCA_052040795.1` is the GenBank pair of
  `GCF_052040795.1`. It pays a chain-projection channel and a second assembly
  per species to buy the alignment §3.2 already forbids.
  Consequence handed back to T-human-008: zebrafish counts with the twelve
  species that have no usable public alignment, not with the two on older
  assemblies, so the number is thirteen of twenty for training. Fugu is
  unaffected; an informant set at inference on a held-out target is permitted.
  Verified: `score.py --self-test` 69/69 still passes, `leakage_check.py`
  still reports 0 violations and the panel is unchanged apart from the
  `notes` cell. No data fetched beyond metadata JSON and two directory
  listings.
  Task stays in `review`; PR #5 updated (commit f23efb9).
  NOT done, unchanged: open item 2, Helixer or Tiberius on one vertebrate --
  the input shape no AUGUSTUS run reaches. NEXT: that, and review feedback on
  PR #5 when it arrives.
- 2026-09-09 lenin: run 8. Closed the standing item (open item 2, "Helixer or
  Tiberius on one vertebrate") and open item 10 with it. Ran Helixer 0.3.7
  from the published container on one consumer GPU: T. rubripes (vertebrate
  model, 384 Mb, 91.6 min), N. crassa (fungi, 41 Mb, 5.8 min) and
  S. cerevisiae (fungi, 12 Mb, 1.9 min). That is the input shape no AUGUSTUS
  run reaches -- gene/mRNA/exon/CDS plus UTRs and NO stop_codon features at
  all -- so the section 4.5 convention is decided by the genome probe rather
  than by a feature or a flag. Two defects, both invisible on every earlier
  run, both fixed with a fixture verified to fail without the fix.
  (1) THE SECTION 4.3 DECILE STRATIFICATION WAS NOT REPRODUCIBLE. Scoring the
  same two fugu files twice gave different per-decile donor and acceptor
  counts while every total was identical. A splice site shared by introns of
  two lengths kept whichever the internal set yielded last, which depends on
  the interpreter's hash seed; that is 8,086 of 220,976 fugu donors and 23,387
  of 189,329 human ones (12.4%), with over 1 Mb between the shortest and the
  longest intron at one site. The shortest intron now decides and the count of
  sites the rule touched is reported.
  (2) THE GENOME WINDOW PLAN WAS MADE BEFORE THE STOP-CODON MERGE. Merging a
  stop_codon feature that lies across an intron creates a junction the plan
  never asked for, so its dinucleotide came back `unknown` instead of the
  GT-AG it is -- one predicted intron of the S. pombe cross-parameter run. The
  plan is extended after the merge, with a second FASTA pass only when a merge
  changed a chain.
  OPEN ITEM 10 ANSWERED: AUGUSTUS's S. cerevisiae donor F1 0.394 is not an
  AUGUSTUS artefact. Helixer scores 0.378 on the same genome and both get
  there by over-predicting introns -- 566 and 734 predicted against 281
  reference ones in a 95.3% single-exon genome, donor precision 0.295 and
  0.262. Two unrelated architectures agreeing rules out the tool-specific
  explanation; deciding between "both are bad at not splicing" and
  "some of those introns are real and unannotated" needs the high-confidence
  subset (item 3) or RNA-seq junctions (item 6), not a third predictor.
  Results: N. crassa nucleotide/exon/donor/transcript/locus F1
  0.962/0.772/0.841/0.686/0.906; S. cerevisiae 0.986/0.825/0.378/0.860/0.948;
  T. rubripes 0.919/0.776/0.853/0.238/0.887. N. crassa is the only panel
  species Helixer can be run on without declaring pretraining exposure, so it
  is the only row that is a measurement rather than an upper bound. Fugu's
  transcript F1 is the isoform effect: 46,771 reference transcripts over
  22,090 loci against one prediction per locus.
  Verified: self-test 104 -> 116 checks, all pass, run under several
  PYTHONHASHSEED values; identity runs re-checked at F1 1.0 / MCC 1.0 with 0
  fusions and 0 splits on S. cerevisiae, N. crassa, T. rubripes, C. elegans
  and human; all six results in benchmark/validation/ regenerated with the
  fixed scorer and the three AUGUSTUS numbers unchanged to five decimals.
  docs/benchmark.md sections 3.2, 4.3, 4.5, 6 and 7 rewritten to match.
  Task stays in `review`; PR #5 updated (commit dedbb1a).
  NOT done, open item 2: the blind 3 bp extension is still fixture-only,
  because Helixer's convention is `inside` -- Tiberius is the tool that would
  exercise it if its GTF excludes the stop; no evidence-based pipeline with
  alternative isoforms has been scored; degraded-copy runs still cover 2 of
  20. NEXT: address review feedback on PR #5 when it arrives.
