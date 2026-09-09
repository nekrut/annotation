---
id: T-human-007
title: Design the benchmark: species panel, references, metrics, leakage rules
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T21:30:00Z
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
- 2026-09-09 lenin: run 9. PR #5 still has no review, and nothing open is
  claimable (T-human-006 waits on T-human-003/004, T-human-009 on T-human-006,
  T-human-010 and T-human-011 on this task reaching `done`), so the tick went
  to the largest unexplained number in the run-8 results: fugu's transcript F1
  0.238. It was not, as run 8 assumed, only the isoform-density effect.
  DEFECT FOUND (the fifth from a real run, and the first that needs a deeply
  annotated reference to see): AN EXACT TRANSCRIPT MATCH COULD LOSE THE
  WITHIN-LOCUS PAIRING TO AN ISOFORM THAT MERELY CONTAINS IT. Section 4.4
  paired predicted transcripts to annotated isoforms in descending order of
  shared CDS bases with the reference id as tie-break. An isoform that
  contains the prediction's whole CDS shares *exactly* as many bases with it
  as the isoform the prediction equals, so the two tie and the accession
  string decides. Fugu `rab44` is the case: the prediction is
  `rna-XM_029826178.1` base for base, `rna-XM_011613896.2` is that chain plus
  51 bases at one end, both share 10,749, `011` sorts before `029`, and the
  exact hit was scored as a miss. 343 of fugu's 5,907 exact matches (5.8%) and
  34 of N. crassa's 6,923 were lost this way, each charged as a false positive
  and a false negative. Exactness now outranks overlap; the ordering is in
  4.4 as a specification, not just in the code.
  Results: T. rubripes transcript F1 0.23782 -> 0.25248 (tp 5,564 -> 5,907),
  N. crassa 0.68609 -> 0.68947 (tp 6,889 -> 6,923). Nothing else moved in any
  of the six runs: the three AUGUSTUS runs and Helixer on S. cerevisiae are
  unchanged, as is every non-transcript column of the other two, because a
  single-isoform reference cannot exhibit the bug. So fugu's transcript F1 is
  still mostly the isoform-density effect; 0.014 of it was this.
  ALSO FIXED, a documentation defect the work surfaced: docs/benchmark.md 6
  claimed the self-test checks 116 values while it ran 104. The self-test now
  prints the count it ran (111 with the new fixture) so the sentence cannot
  drift again.
  ANSWERED, for marx's T-human-008 question in 20260909T140006Z-marx-0013:
  the benchmark designates no representative isoform per locus and 4.4 now
  says why. Picking one would pick a winner the reference does not, and the
  isoform rule already removes the penalty deep annotation would carry. A
  window cutter needing one target per locus is making a training decision,
  not a scoring one. The constraint that does bind: against a
  one-prediction-per-locus model the transcript column is bounded by how close
  the chosen target is to *some* annotated isoform.
  Verified: new fixture (two isoforms, one containing the other, prediction
  equals the shorter) fails under the old order and passes under the new;
  self-test 104 -> 111 checks, all pass under PYTHONHASHSEED 0/1/2/42/1337;
  identity runs re-checked at F1 1.0 and MCC 1.0 with 0 fusions and 0 splits
  on S. cerevisiae, N. crassa, T. rubripes, C. elegans and human (human
  131,442 transcript tp), with shared-donor counts 0/73/8,086/2,120/23,387
  unchanged; all six benchmark/validation/ JSONs regenerated.
  docs/benchmark.md 4.4, 6 and 7 item 2 and benchmark/validation/README.md
  updated. Task stays in `review`; PR #5 updated (commit b89f938).
  NOT done, unchanged: open item 2's remaining gaps -- the blind 3 bp
  extension is still fixture-only, no prediction with more than one isoform
  per locus has been scored (so the false-positive branch of the 4.4 pairing
  is fixture-only too), no evidence-based pipeline, degraded-copy runs still
  cover 2 of 20. NEXT: address review feedback on PR #5 when it arrives.
- 2026-09-09T15:30Z (lenin): lease renewed to 17:30Z; task stays in `review`,
  PR #5 still has no review from another agent. Bounded work went to the
  largest remaining gap in open item 2: no submission had ever carried more
  than one isoform per locus, so the 4.4 branch that charges an unmatched
  extra prediction as a false positive had only ever run on fixtures.
  SCORED: GENCODE 50 (Ensembl 116, primary-assembly GFF3, sha256 272f9972...)
  against the human RefSeq reference, sequence names rewritten from the
  GRCh38.p14 assembly report and nothing else changed. 370,476 scored chains
  against the reference's 131,442 -- 2.8 isoforms for every one the reference
  has. 64 s, 1.6 GB, identical under PYTHONHASHSEED 0 and 1. Not an accuracy
  measurement (both files are human curation of the same genome); it is a real
  multi-isoform input of the shape an evidence-based pipeline produces.
  Result: nucleotide F1 0.933, exon 0.698, donor 0.908, acceptor 0.826,
  transcript 0.290 (precision 0.190, 300,036 unmatched predicted chains),
  locus 0.968, start 0.637, stop 0.360.
  FIXED, two defects it found, both invisible to a RefSeq-shaped submission:
  (1) the section 4 biotype filter read only `gene_biotype`; GENCODE writes
  `gene_type` and Ensembl writes `biotype`. Because the filter is applied to
  the prediction as well as the reference, reading one spelling dropped the
  reference's 390 immunoglobulin and T-cell receptor transcripts as
  unanswerable while keeping GENCODE's 421 and charging them as false
  positives -- penalising a submission for answering a question section 4
  forbids the reference to ask. All three spellings are now read; locus F1
  0.959 -> 0.968 (false-positive loci 1,133 -> 723), nucleotide 0.931 ->
  0.933. (2) `predicted_transcripts_not_scored` pooled the filter's drops with
  predictions on sequences the reference does not have, which is the thing its
  warning tells the submitter to go and fix; it read 434 where 13 transcripts
  are actually off-panel. It now counts sequence exclusion only.
  Verified: self-test 111 -> 118 checks, all passing under PYTHONHASHSEED
  0/1/42; the new fixture fails on the pre-fix code with locus fp 2 and
  transcript fp 4. The six earlier validation runs all have
  `predicted_transcript_selection.dropped == 0`, so none of them moves.
  The 21% excess of distinct CDS acceptors over donors in GENCODE (257,153 vs
  212,078, against RefSeq's 193,359 vs 188,913) was recounted independently
  from the GFF3 outside the scorer and agrees exactly, so acceptor F1 0.826
  below donor F1 0.908 is annotation depth, not a scoring artefact.
  docs/benchmark.md 4, 6 and 7 item 2 and benchmark/validation/README.md
  updated; PR #5 updated (commit f4ea4fc).
  NOT done, unchanged: the blind 3 bp extension is still fixture-only, no
  prediction whose sequence set genuinely diverges from the reference's, no
  evidence-based pipeline run end to end, no predictor emitting partial genes
  at contig ends, degraded-copy runs still cover 2 of 20.
  NEXT: address review feedback on PR #5 when it arrives.
- 2026-09-09T17:20Z (lenin): run 10. Lease renewed to 19:20Z; task stays in
  `review`, PR #5 still has no review from another agent. The tick went to the
  last untouched shape in open item 2: no submission had ever declared its own
  partial genes, so all of section 4.5's predicted-partial handling was
  fixture-only -- and the fixtures were built around the *reference*
  convention.
  RAN: AUGUSTUS 3.5.0 `--genemodel=partial` twice. *T. thermophila*
  (`--species=tetrahymena`, one invocation per scaffold over its 1,158
  scaffolds, genetic code 6) and *A. mellifera* (`--species=honeybee1`, 177
  sequences, 36.4 min on 22 cores, 570 MB). Both are in AUGUSTUS's own
  training set, so neither is a measurement.
  Results: T. thermophila nucleotide F1 0.420, exon 0.328, donor 0.392,
  transcript 0.208, locus 0.536, start 0.310, stop 0.434; A. mellifera 0.891 /
  0.693 / 0.811 / 0.225 / 0.793 / 0.402 / 0.618.
  DEFECT 1 (the sixth from a real run): COUNTING A REUSED TRANSCRIPT ID WAS
  NOT ENOUGH; IT HAD TO BE RESOLVED. Run 2 added
  `predicted_conflicting_transcript_ids` and a warning but still welded the
  colliding chains. On 1,158 scaffolds that is not an edge case: 214 ids
  collide and the file drops to 218 scored chains from 10,355, nucleotide F1
  0.020 against 0.420. A CDS chain lies on one sequence and one strand by
  construction, so ids are now keyed by (sequence, strand). Control: the same
  AUGUSTUS run repeated with `--uniqueGeneId=true` scores identically block
  for block, which is what the new fixture asserts instead of a count.
  DEFECT 2: `predicted_partial_5prime`/`_3prime` COULD NOT FIRE. They read
  `start_range=`/`end_range=`, a reference convention no predictor writes, so
  they were structurally zero on every prediction ever scored here. A
  GTF-lineage predictor states a truncated end by omitting the
  `start_codon`/`stop_codon` feature instead. T. thermophila states 250 chains
  with no start codon and 69 with no stop; reading that removes 208
  false-positive starts and 22 false-positive stops. `predicted_partial_source`
  now says which convention the counts came from.
  CROSS-CHECKED, because the inference is a guess unless something independent
  agrees: GENCODE 50 writes codon features *and* tags incomplete ends
  `cds_start_NF`/`cds_end_NF`. Omissions 13,203 and 19,535; tags 13,226 and
  19,858; agreement 99.7% and 99.2% over 370,910 transcripts. That is why the
  GENCODE row's codon columns moved -- start F1 0.637 -> 0.750, stop 0.360 ->
  0.408, every other column identical to five decimals.
  DEFECT 3: A 3'-PARTIAL GENE COULD BE HANDED A STOP CODON IT NEVER PREDICTED.
  The bare-GTF case and the truncated-gene case both reach the 4.5 3 bp
  extension as a chain with no `stop_codon` feature and need opposite
  treatment. A. mellifera is the run that exercises it (`honeybee1` excludes
  the stop, `tetrahymena` includes it): 13 chains skipped, worth 8
  false-positive stops, 2 false-positive terminal exons, 27 false-positive
  bases. The fixture shows the worse form, where the invented 3 bp land on the
  reference's real stop and score a *true positive*.
  ALSO DOCUMENTED: "AUGUSTUS excludes the stop codon" is not a property of
  AUGUSTUS. `stopCodonExcludedFromCDS` is a line in each *species* parameter
  file; of the 166 shipped with 3.5.0, 44 set it true and 122 false. Section
  4.5 now carries the one-line census that produces those counts.
  Verified: self-test 118 -> 135 checks under PYTHONHASHSEED 0/1/42/1337, the
  partial fixture confirmed to fail on the pre-fix code (stop tp 3 of 3 where
  2 is right); identity runs F1 1.0 and MCC 1.0 with 0 fusions and 0 splits on
  S. cerevisiae, N. crassa, C. elegans, T. rubripes, T. thermophila,
  A. mellifera and human; all eight benchmark/validation/ JSONs regenerated,
  and the six earlier ones are unchanged in every value apart from the two new
  codon fields.
  docs/benchmark.md 4.5, 6, 7 items 2 and 13 and benchmark/validation/README.md
  updated; PR #5 updated (commit 025ef2d).
  NOT done, unchanged: the blind 3 bp extension is still fixture-only (Helixer
  is `inside` and both AUGUSTUS conventions ship the feature), no prediction
  whose sequence set genuinely diverges from the reference's, no
  evidence-based pipeline end to end, degraded-copy runs still cover 2 of 20.
  Open item 13 is closed.
  NEXT: address review feedback on PR #5 when it arrives.
- 2026-09-09T18:55Z (lenin): lease renewed to 20:15Z. Run 9: the control run,
  panel-wide, and a defect it found in itself.
  ADDED `benchmark/degrade.py` -- a seeded, known-perturbation copy of a
  reference: delete 10% of transcripts, move the downstream boundary in
  transcription order of 10% of CDS segments 3 bp further downstream. An
  identity run makes every metric 1.0 by construction, so it cannot tell a
  correct scorer from one dropping the same thing from both sides; a degraded
  copy can. Degraded-copy coverage 2 of 20 species -> 20 of 20. Results in
  `benchmark/validation/degraded/` (one scored JSON and one degrade summary
  per species) plus its README and the §3.3 declaration they were run under.
  THE PERTURBATION IS ASYMMETRIC ON PURPOSE, AND THAT IS THE TEST. Which
  coordinate a "downstream" shift moves depends on the strand -- `end` on +,
  `start` on - -- so it lands on stop codons and donors and never on start
  codons or acceptors. Every species reproduces it: donor F1 0.811-0.870
  against acceptor F1 0.934-0.989.
  LOCUS PRECISION AND START-CODON PRECISION ARE EXACTLY 1.0 ON ALL 20, over
  775,253 scored reference transcripts: zero false-positive loci, zero
  false-positive start codons. The earlier 2-species run reported human start
  F1 0.912 *with* false positives, which was correct for the ad-hoc script
  behind it -- it moved `end` regardless of strand, so on half the annotation
  it perturbed starts while the text said stops. Hence a reviewable file.
  FUSION AND SPLIT ARE NOT INVARIANT under a perturbation that is neither.
  Six species report 1-2 fusions, six report 1-8 splits. Ablating the two
  perturbations at one seed separates them exactly: every split comes from the
  deletions (human 8, A. mellifera 2, X. tropicalis 1), every fusion from the
  3 bp shifts (human 1, A. mellifera 1). Both are §4.4 as specified, but a
  low-single-digit count is inside the noise a *correct* submission produces.
  Report it, do not rank on it.
  DEFECT IN degrade.py, FOUND BY THE ABLATION AND NOT BY THE COMBINED RUN: the
  shift was drawn only for CDS lines surviving the deletion pass, making the
  shift stream a function of the drop rate, so the ablation was not a
  decomposition of the combined run. T. rubripes showed it -- a combined
  fusion neither single perturbation could produce, which is impossible since
  combined CDS blocks are a subset of shift-only's. Fixed; the self-test now
  asserts the invariant over four seeds and fails on the pre-fix code.
  Verified: degrade.py --self-test 27 checks and score.py --self-test 135
  checks under PYTHONHASHSEED 0/1/42/1337; all four ablated species reproduce
  the combined fusion and split counts exactly after the fix; leakage_check
  still 0 violations. Scoring cost 0.68 s / 44 MB (S. pombe) to 52 s / 1.04 GB
  (D. rerio); human 43 s / 1.06 GB.
  docs/benchmark.md §6 rewritten with the 20-row table and the four findings,
  §7 item 2 updated; benchmark/README.md and validation/README.md updated.
  PR #5 updated (commit 834a2fa).
  NOT done, unchanged: no BUSCO/OMArk or cost columns; the blind 3 bp
  extension is still fixture-only; no prediction whose sequence set genuinely
  diverges from the reference's; no evidence-based pipeline end to end. Item
  2's coverage half is closed; what is left is real predictor output.
  NEXT: address review feedback on PR #5 when it arrives.
- 2026-09-09T19:25Z (lenin): lease renewed to 21:30Z. Run 10: the 3 Gb genome,
  and what it found about the section 4.3 stratifications.
  ITEM 9 CLOSED. `fetch.py --what fasta` had never run above 391 Mb. Human:
  `GCF_000001405.40_GRCh38.p14_genomic.fna.gz`, 972,898,531 bytes for 3.1 Gb,
  MD5-verified against the NCBI manifest in 19.6 s. What is left there is
  eleven more species, not scale.
  `score.py --genome` AT HUMAN SCALE for the first time -- the case the
  streaming window reader was written for. Identity 3,101,538,863 scored bp
  over 102 sequences in 67 s / 1.35 GB; GENCODE 50 (347,110 predicted introns
  against the reference's 218,446) in 88 s / 2.12 GB, against 43 s / 1.06 GB
  without it. NOT ONE WINDOW CAME BACK UNSERVED: no `unknown` dinucleotide
  class and no dropped GC bin in either run, over the FASTA's 705 records.
  Every GENCODE metric is identical with and without `--genome`; what the
  genome adds is the two strata and a second reading of the stop-codon
  convention, which agrees with the file's own `stop_codon` features (349,754
  chains inside, 579 outside, 18,253 neither). Human reference splice census:
  215,956 GT-AG, 1,939 GC-AG, 219 AT-AC, 332 other, so 1.14% non-GT-AG.
  THE FINDING: THE THREE 4.3 STRATIFICATIONS HAVE THREE DENOMINATORS, AND ONE
  WAS SINGLE-SIDED. `by_dinucleotide` is per intron, the deciles and local GC
  are per site, and 4.3 read as if all three were per donor and per acceptor.
  Invisible on S. cerevisiae (281 introns over 281 donors); on human 218,446
  introns over 188,913 donors and 193,359 acceptors, so the dinucleotide row
  sums to 15.6% more than the donor total printed beside it. Per-intron is the
  only defensible unit -- a class is the pair -- and the sites that sit in two
  classes at once are now counted: 106 donors and 644 acceptors in the human
  reference, 129 and 2,515 in GENCODE, 185 and 260 in T. rubripes, 47 and 88
  in A. mellifera. The 6:1 acceptor-to-donor ratio is the mechanism: a shared
  donor usually keeps its class because both introns end AG, a shared acceptor
  changes it whenever one intron starts GC. `splice.strata_units` now states
  each table's unit, `by_local_gc` gains an acceptor table (it was donors
  only; the old table is `by_local_gc.donor` value for value), and
  `reference/predicted_sites_multiple_dinuc_classes` report the counts.
  ALSO: human is the counter-case for open item 11 -- its 188,913 donors
  spread 7,913 / 50,492 / 43,960 / 45,405 / 41,143 across the five fixed GC
  bands, so the bands that degenerate on P. falciparum carry signal on a
  vertebrate.
  Verified: self-test 135 -> 145 checks under PYTHONHASHSEED 0/1/42/1337, the
  new fixture confirmed to fail on the pre-change code (`gc acceptor sites:
  got 0, want 2`); degrade.py 27 checks; leakage_check 0 violations; all nine
  validation JSONs regenerated from the original prediction files with NO
  EXISTING FIELD MOVING in any of them.
  docs/benchmark.md 4.3, 6, 7 items 9 and 11, benchmark/README.md and
  validation/README.md updated; PR #5 updated (commit 795e0f8).
  NOT done, unchanged: no BUSCO/OMArk or cost columns; the blind 3 bp
  extension is still fixture-only; no prediction whose sequence set genuinely
  diverges from the reference's; no evidence-based pipeline end to end.
  NEXT: address review feedback on PR #5 when it arrives.
