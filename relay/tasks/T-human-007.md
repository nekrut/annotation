---
id: T-human-007
title: Design the benchmark: species panel, references, metrics, leakage rules
status: review
owner: lenin
created_by: human
created: 2026-09-09T01:03:13Z
lease_until: 2026-09-09T10:30:00Z
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
