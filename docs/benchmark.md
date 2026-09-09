# Benchmark design: species panel, references, metrics, and leakage rules

Task [T-human-007](../relay/tasks/T-human-007.md). Owner: `lenin`.
Status: first draft, for review.

Every number in this document was either read from the NCBI Datasets API on
2026-09-09 or computed from the reference annotation of the assembly named in
the same row by [`benchmark/annotation_stats.py`](../benchmark/annotation_stats.py).
Nothing here is quoted from a paper's own benchmark. The manifest that
carries all of it is [`benchmark/panel.tsv`](../benchmark/panel.tsv).

## 1. What this benchmark is for, and what it is not

The charter asks for a model that is accurate *across clades* at a small
fraction of current compute. That is a generalization claim and a cost claim,
and neither is measurable on the benchmarks the field currently uses. The
review under `relay/artifacts/T-human-002/review.md` (§5.8) found that no two
methods in its publications table share an evaluation set: GENSCAN reports
exact exons on 1997 vertebrate sets, BRAKER3 transcript-level F1 on 11
species, Tiberius gene-level F1 on human, and the 2025 overview uses G3PO.
Comparing published numbers across those is meaningless.

So this benchmark has exactly three jobs.

1. **Make cross-clade generalization measurable.** One panel, one metric
   definition, one set of held-out species that no training run may touch.
2. **Make cost comparable.** Wall clock and peak memory measured on the same
   hardware over the same genomes, normalized per megabase.
3. **Refuse a single headline number.** He and Florea show that every
   sequence foundation model peaks on the exon class best represented in its
   training data and "decreases drastically" on the rest
   ([10.64898/2026.02.22.707219](https://doi.org/10.64898/2026.02.22.707219)),
   and GENATATORs argues that standard per-token and per-sequence scores
   "fail to capture the challenges of real-world gene annotation"
   ([10.64898/2026.06.17.732686](https://doi.org/10.64898/2026.06.17.732686)).
   Both point the same way: stratified reporting is mandatory here (§4.8).

It is **not** a leaderboard for picking a production annotator, and it is not
a training set. A model may train on the ten `train` species; it may not
train on the panel.

## 2. The panel

Twenty genomes. The selection rule was: cover the eukaryotic clades the
charter names, and inside that constraint maximize the spread of the two
variables that break existing tools — intron length distribution and GC
composition — while keeping every genome small enough that a full-genome
prediction run fits on a laptop or one consumer GPU.

### Table 1. The panel

| species | common name | clade | assembly | size | GC% | coding genes | genes/Mb | CDS% | split |
|---|---|---|---|---|---|---|---|---|---|
| *Gallus gallus* | chicken | bird | `GCF_016699485.2` | 1.05 Gb | 42 | 18,023 | 17.1 | 14.4 | `heldout` |
| *Nematostella vectensis* | starlet sea anemone | cnidarian | `GCF_932526225.1` | 269 Mb | 40.5 | 19,231 | 71.4 | 25.5 | `heldout` |
| *Schizosaccharomyces pombe* | fission yeast | fungus-ascomycete | `GCF_000002945.2` | 13 Mb | 36 | 5,124 | 407.6 | 57.3 | `heldout` |
| *Plasmodium falciparum* | malaria parasite | protist-apicomplexan | `GCF_000002765.6` | 23 Mb | 19.5 | 5,282 | 226.8 | 52.8 | `heldout` |
| *Tetrahymena thermophila* | ciliate | protist-ciliate | `GCF_000189635.1` | 103 Mb | 22.5 | 26,996 | 262.1 | 49.4 | `heldout` |
| *Ciona intestinalis* | vase tunicate | tunicate | `GCF_018327825.1` | 140 Mb | 36 | 15,365 | 109.5 | 24.8 | `heldout` |
| *Xenopus tropicalis* | western clawed frog | amphibian | `GCF_000004195.4` | 1.45 Gb | 40.5 | 21,826 | 15.0 | 6.7 | `train` |
| *Danio rerio* | zebrafish | fish | `GCF_052040795.1` | 1.45 Gb | 37 | 28,415 | 19.6 | 14.9 | `train` |
| *Saccharomyces cerevisiae* | baker's yeast | fungus-ascomycete | `GCF_000146045.2` | 12 Mb | 38.5 | 6,021 | 498.8 | 73.1 | `train` |
| *Neurospora crassa* | red bread mould | fungus-ascomycete | `GCF_000182925.2` | 41 Mb | 48.5 | 9,757 | 237.8 | 41.3 | `train` |
| *Drosophila melanogaster* | fruit fly | insect | `GCF_000001215.4` | 144 Mb | 42 | 13,962 | 97.2 | 42.6 | `train` |
| *Mus musculus* | house mouse | mammal | `GCF_000001635.27` | 2.73 Gb | 42 | 22,198 | 8.1 | 7.3 | `train` |
| *Caenorhabditis elegans* | roundworm | nematode | `GCF_000002985.6` | 100 Mb | 35.5 | 19,983 | 199.3 | 42.8 | `train` |
| *Arabidopsis thaliana* | thale cress | plant-dicot | `GCF_000001735.4` | 119 Mb | 36 | 27,562 | 231.3 | 52.6 | `train` |
| *Zea mays* | maize | plant-monocot | `GCF_902167145.1` | 2.18 Gb | 47 | 34,313 | 15.7 | 3.6 | `train` |
| *Dictyostelium discoideum* | slime mould | protist-amoebozoan | `GCF_000004695.1` | 34 Mb | 22.5 | 13,289 | 389.1 | 62.0 | `train` |
| *Takifugu rubripes* | fugu | fish | `GCF_901000725.3` | 384 Mb | 45.5 | 22,090 | 57.5 | 27.1 | `heldout_paired` |
| *Apis mellifera* | honey bee | insect | `GCF_003254395.2` | 225 Mb | 32.5 | 9,935 | 44.1 | 25.5 | `heldout_paired` |
| *Homo sapiens* | human | mammal | `GCF_000001405.40` | 3.10 Gb | 41 | 20,076 | 6.5 | 9.6 | `heldout_paired` |
| *Oryza sativa* | asian rice | plant-monocot | `GCF_034140825.1` | 386 Mb | 43.5 | 29,427 | 76.3 | 15.8 | `heldout_paired` |

### Table 2. Intron and exon structure

| species | unique introns | intron p10 | median | p90 | p99 | max | exons/tx | single-exon tx % |
|---|---|---|---|---|---|---|---|---|
| *Schizosaccharomyces pombe* | 5,253 | 41 | **56** | 162 | 407 | 2,526 | 2.03 | 51.7 |
| *Caenorhabditis elegans* | 114,350 | 45 | **65** | 809 | 3,747 | 100,912 | 6.83 | 4.2 |
| *Tetrahymena thermophila* | 94,194 | 52 | **72** | 282 | 873 | 6,975 | 4.49 | 31.8 |
| *Neurospora crassa* | 17,862 | 56 | **76** | 225 | 629 | 1,884 | 2.85 | 16.9 |
| *Arabidopsis thaliana* | 133,552 | 78 | **101** | 355 | 922 | 146,982 | 5.99 | 18.2 |
| *Drosophila melanogaster* | 58,679 | 57 | **102** | 3,621 | 27,362 | 268,107 | 5.96 | 9.2 |
| *Dictyostelium discoideum* | 16,808 | 73 | **104** | 220 | 767 | 2,842 | 2.26 | 32.0 |
| *Plasmodium falciparum* | 8,445 | 96 | **142** | 248 | 693 | 2,425 | 2.63 | 45.3 |
| *Zea mays* | 166,366 | 81 | **146** | 1,448 | 11,431 | 398,243 | 6.82 | 11.1 |
| *Saccharomyces cerevisiae* | 296 | 75 | **147** | 511 | 1,571 | 2,483 | 1.05 | 95.3 |
| *Apis mellifera* | 77,561 | 72 | **158** | 3,746 | 59,490 | 596,047 | 10.36 | 1.2 |
| *Oryza sativa* | 139,097 | 82 | **161** | 1,065 | 4,232 | 228,480 | 6.38 | 13.3 |
| *Takifugu rubripes* | 250,741 | 76 | **170** | 1,729 | 13,498 | 1,065,045 | 14.36 | 1.3 |
| *Ciona intestinalis* | 133,190 | 78 | **345** | 897 | 6,287 | 198,088 | 10.72 | 8.1 |
| *Nematostella vectensis* | 164,908 | 115 | **459** | 1,968 | 10,150 | 99,820 | 11.44 | 7.3 |
| *Gallus gallus* | 218,705 | 123 | **889** | 8,164 | 55,595 | 983,230 | 14.71 | 1.2 |
| *Xenopus tropicalis* | 221,818 | 153 | **1149** | 7,487 | 46,894 | 1,192,258 | 13.77 | 3.3 |
| *Danio rerio* | 302,826 | 89 | **1150** | 7,610 | 62,699 | 1,090,140 | 15.08 | 1.2 |
| *Mus musculus* | 247,313 | 167 | **1528** | 12,823 | 88,830 | 1,067,346 | 13.59 | 2.8 |
| *Homo sapiens* | 304,911 | 176 | **1781** | 16,983 | 107,644 | 1,160,411 | 13.87 | 1.2 |

### Table 3. Reference annotation

| species | annotation | provider | released | transcripts (mRNA) |
|---|---|---|---|---|
| *Neurospora crassa* | Annotation submitted by Broad Institute | Broad Institute | 2023-04-03 | 10,812 |
| *Drosophila melanogaster* | FlyBase Release 6.54 | FlyBase | 2023-12-26 | 30,802 |
| *Homo sapiens* | GCF_000001405.40-RS_2025_08 | NCBI RefSeq | 2025-08-01 | 145,477 |
| *Mus musculus* | GCF_000001635.27-RS_2024_02 | NCBI RefSeq | 2024-02-01 | 97,339 |
| *Gallus gallus* | NCBI Gallus gallus Annotation Release 106 | NCBI RefSeq | 2021-12-23 | 68,684 |
| *Danio rerio* | GCF_052040795.1-RS_2026_07 | NCBI RefSeq | 2026-07-20 | 94,006 |
| *Takifugu rubripes* | GCF_901000725.3-RS_2026_03 | NCBI RefSeq | 2026-03-02 | 46,771 |
| *Xenopus tropicalis* | NCBI Xenopus tropicalis Annotation Release 104 | NCBI RefSeq | 2019-12-06 | 45,100 |
| *Ciona intestinalis* | GCF_018327825.1-RS_2026_01 | NCBI RefSeq | 2026-01-16 | 19,614 |
| *Nematostella vectensis* | NCBI Nematostella vectensis Annotation Release 101 | NCBI RefSeq | 2022-06-10 | 32,370 |
| *Apis mellifera* | NCBI Apis mellifera Annotation Release 104 | NCBI RefSeq | 2018-09-13 | 23,471 |
| *Oryza sativa* | GCF_034140825.1-RS_2024_06 | NCBI RefSeq | 2024-06-21 | 42,933 |
| *Zea mays* | GCF_902167145.1-RS_2025_02 | NCBI RefSeq | 2025-02-14 | 57,345 |
| *Dictyostelium discoideum* | Annotation submitted by NCBI RefSeq | NCBI RefSeq | 2023-04-03 | 13,315 |
| *Tetrahymena thermophila* | Annotation submitted by NCBI RefSeq | NCBI RefSeq | 2015-06-09 | 26,996 |
| *Plasmodium falciparum* | Annotation submitted by Plasmodium falciparum Genome Sequencing Consortium | Plasmodium falciparum Genome Sequencing Consortium | 2020-09-10 | 5,354 |
| *Schizosaccharomyces pombe* | Annotation submitted by PomBase | PomBase | 2024-08-16 | 5,166 |
| *Saccharomyces cerevisiae* | SGD R64-5-1 | SGD | 2026-07-10 | 6,039 |
| *Arabidopsis thaliana* | Annotation submitted by TAIR and Araport | TAIR and Araport | 2022-10-20 | 52,177 |
| *Caenorhabditis elegans* | WormBase WS298 | WormBase | 2025-12-01 | 30,562 |

### 2.1 What the panel spans

Read down the median-intron column of Table 2. It runs from **56 bp**
(*S. pombe*) to **1,781 bp** (human), a 32-fold range, and the 99th
percentile runs from **407 bp** to **107,644 bp**, a 264-fold range. That
tail is the thing that breaks GHMMs: Helixer's headline claim was that its
predictions are "much less sensitive to the length of the genome"
([10.1093/bioinformatics/btaa1044](https://doi.org/10.1093/bioinformatics/btaa1044)),
and BRAKER3's gains were "most pronounced for species with large and complex
genomes" ([10.1101/gr.278090.123](https://doi.org/10.1101/gr.278090.123)).
Both statements are about this axis.

GC runs from **19.5%** (*P. falciparum*) to **48.5%** (*N. crassa*).
GENSCAN's answer to composition in 1997 was to fit separate parameter sets
per C+G region ([10.1006/jmbi.1997.0951](https://doi.org/10.1006/jmbi.1997.0951)),
and minisplice still finds GC-rich introns specific to mammals and birds in
2026 ([10.1186/s13015-025-00293-7](https://doi.org/10.1186/s13015-025-00293-7)),
so composition remains entangled with clade at the splice signal itself. Any
model claiming clade independence has to be shown stratified by GC.

Gene density runs from **6.5 genes/Mb** (human) to **498.8** (*S. cerevisiae*),
and coding fraction from **3.6%** (maize) to **73.1%** (*S. cerevisiae*).
Maize and *S. cerevisiae* are the two ends of the "how much of the genome is
signal" axis; a model tuned on one will have a badly calibrated prior on the
other.

Structural complexity runs from **1.05 exons per transcript** in
*S. cerevisiae* (95.3% single-exon) to **15.08** in zebrafish. The
*S. cerevisiae* row is deliberately close to a degenerate case: a method that
cannot get a near-intronless genome right is not a gene finder.

Three picks need their own justification.

- ***Tetrahymena thermophila*** uses **nuclear genetic code 6**: TAA and TAG
  are read as glutamine and only TGA terminates. Every stop-codon prior in
  every tool in the review is wrong on this genome. It is in the panel as an
  explicit test of whether a model has learned "stop codon" as a fixed token
  or as a property it can infer. Its annotation is the weakest in the panel
  (§2.3), so it is scored but never used to rank.
- ***Plasmodium falciparum*** at 19.5% GC, and ***Dictyostelium discoideum***
  at 22.5%, are outside the compositional range of every clade-specific model
  the review found. EGAPx declares fungi, protists and nematodes out of scope
  in its own README, naming its supported taxa as Chordata, Arthropoda,
  Echinodermata, Mollusca, Cnidaria, monocots and eudicots; Tiberius's 2026
  multi-clade extension reaches "92% of currently available eukaryotic
  assemblies" ([10.64898/2026.04.24.720536](https://doi.org/10.64898/2026.04.24.720536)).
  A panel drawn from that 92% cannot measure species independence, so this
  one deliberately oversamples the excluded tail: three protists, three
  fungi, one nematode, one cnidarian, one tunicate.
- ***Ciona intestinalis*** and ***Nematostella vectensis*** are the
  early-branching animals. Ciona is the only invertebrate chordate; it is the
  outgroup to every vertebrate in the panel, and at 109.5 genes/Mb with a
  345 bp median intron it is a compact genome with vertebrate-like gene
  content.

### 2.2 How the numbers were produced

`benchmark/annotation_stats.py` streams the RefSeq `*_genomic.gff.gz` for the
accession in the row. Introns are derived from the exons of `mRNA` features
and de-duplicated by (seqid, start, end, strand), so an intron shared by
several isoforms counts once; quantiles are nearest-rank on that unique list.
`cds_fraction_pct` is the one column `annotation_stats.py` does not reproduce
on its own: its numerator is the script's `cds_total_bp`, but its denominator
is the Datasets API's `genome_bp`, not the total of the GFF's own
`##sequence-region` lines, which also counts the organelles. For
*S. cerevisiae* that is 8,825,064 bp over 12,071,326 bp = 73.1% in the panel
against 8,825,064 over 12,157,105 = 72.6% from the file, the difference being
the 85,779 bp mitochondrion (*S. pombe*: 57.3 against 57.2). The script now
prints `cds_fraction_of_sequence_region_pct` beside `sequence_region_bp` so
both denominators are visible. Neither is the denominator §4 scores against,
which is the selected sequence set reported in `sequence_selection`.

`protein_coding_genes` in Table 1 is NCBI's own count from the Datasets API;
`gff_gene_features` in `panel.tsv` is the count of `protein_coding` gene
features in the GFF, which is larger for assemblies carrying alt loci and
patch scaffolds (human: 20,076 versus 23,306). Both are kept because the
difference is exactly the kind of thing that silently inflates a denominator.

The whole panel's annotations are 318 MB of gzipped GFF and took 21.5 s of
CPU to summarize; the peak resident set was 428 MB, on human.
Reproduce with:

```
python3 benchmark/fetch.py --all --what gff --dest /tmp/panel
for g in /tmp/panel/*/*_genomic.gff.gz; do python3 benchmark/annotation_stats.py "$g"; done
```

### 2.3 Which annotation is truth, and where it is wrong

Truth for each species is the annotation named in Table 3, pinned by the MD5
in `panel.tsv`. Twelve are NCBI RefSeq; the rest are the community
annotations RefSeq itself redistributes (FlyBase, WormBase, SGD, PomBase,
TAIR/Araport, Broad). Community annotations are preferred where they exist
because they are curated against experimental evidence rather than produced
by the class of pipeline under test — scoring a predictor against Gnomon
output would measure agreement with a competitor, not correctness.

The reference is not ground truth, and the size of the error is now
measurable:

- Vertebrate selenoprotein genes, where UGA is recoded rather than a stop,
  are well annotated for only **11% of genes in Ensembl and 5% in NCBI
  GenBank** ([10.1371/journal.pcbi.1013885](https://doi.org/10.1371/journal.pcbi.1013885)).
  This affects human, mouse, chicken, *Xenopus*, zebrafish and fugu here.
- Manual curation of one *Pristionchus pacificus* strain corrected **more
  than 7,500 gene models, about 24% of the annotation**
  ([10.64898/2026.02.18.706511](https://doi.org/10.64898/2026.02.18.706511)),
  attributed to assembly errors, artificial transcript fusions and
  unexpressed genes. *C. elegans* is the panel's nematode and is far better
  curated, but the magnitude is the warning.
- GAP-MS finds hundreds of peptide-supported coding loci absent from
  reference annotations across nine crops
  ([10.64898/2026.03.17.712294](https://doi.org/10.64898/2026.03.17.712294)).
  Rice and maize are in the panel.

Per-species caveats that a scorer must carry:

| species | caveat |
|---|---|
| *Homo sapiens* | GRCh38.p14 carries alt loci and patch scaffolds; score on the primary assembly only, or the same gene is counted several times. |
| *Gallus gallus* | annotation release 106 (2021-12-23) is the oldest vertebrate annotation in the panel. |
| *Xenopus tropicalis* | release 104 (2019-12-06) is the oldest annotation of any panel species still used for ranking. |
| *Tetrahymena thermophila* | scaffold-level assembly (1,158 scaffolds), annotation submitted 2015-06-09, no non-coding or pseudogene calls at all — the 26,996 gene models are all `protein_coding`, which is itself implausible. Report, never rank. |
| *Saccharomyces cerevisiae* | 296 unique introns genome-wide; exon-level metrics have almost no support and gene-level metrics dominate. |
| *Danio rerio* | GRCz12ab annotation is from 2026-07-20 and will move; the MD5 in `panel.tsv` is the pin. The assembly is deliberately the one no public alignment is built on; see §2.4. |
| all | isoform choice matters: a locus with 40 annotated transcripts is easy to hit at gene level and hard at transcript level (§4.4). |

**A high-confidence subset is required and does not exist yet.** Above some
accuracy level the benchmark stops measuring the model and starts measuring
its ability to reproduce known annotation errors. §7 carries this as the
largest open item.

### 2.4 Assembly choice when the annotation and the alignment disagree (*Danio rerio*)

Zebrafish is the one panel species whose current reference annotation and
only public multiple alignment sit on different assemblies. marx raised it
as `docs/data-sources.md` §7 item 6; the call is mine and the panel
**keeps GRCz12ab** (`GCF_052040795.1`). The reasoning is recorded because
the same question recurs every time a telomere-to-telomere assembly lands
ahead of the alignments built on its predecessor.

The alternative is GRCz11 (`GCF_000002035.6`). It is a real option: it
carries a full RefSeq annotation (`GCF_000002035.6-RS_2024_08`, 2024-08-15,
27,158 protein-coding genes) and it is still the zebrafish assembly in
Ensembl 116 (`rest.ensembl.org/info/assembly/danio_rerio`, 2026-09-09).
Moving to it buys exactly one thing: the Ensembl fish EPO.

**That alignment is unusable here regardless of assembly.** Both fish
species sets — EPO with 32 species and EPO-extended with 65, from
`rest.ensembl.org/info/compara/species_sets/{EPO,EPO_EXTENDED}` on
2026-09-09 — contain `takifugu_rubripes`, and fugu is the `heldout_paired`
species paired with zebrafish. EPO is jointly inferred, so §3.2 channel 3
requires it to be *rebuilt* with fugu removed before a zebrafish training
window may be cut from it; dropping the row does not qualify. Zebrafish is
rebuild-only for the comparative arm on GRCz11 exactly as it is on GRCz12ab.
The contamination is not incidental either: the zebrafish/fugu pair exists to
test a 6.8-fold median-intron shift inside one subclass (Table 2), which is
the axis a fugu-informed alignment leaks on.

Nothing else is on either assembly. UCSC `danRer11` (= GRCz11) has **no
multiz, phastCons or phyloP track at all**: `api.genome.ucsc.edu
/list/tracks?genome=danRer11` returns 40 tracks, none of them conservation,
and `hgdownload.soe.ucsc.edu/goldenPath/danRer11/` has no `*way` directory,
only the `vsHg38`, `vsMm39`, `vsMm10` and `vsGalGal6` pairwise chains (marx's
`docs/data-sources.md` §2.1 records the same absence). GRCz12ab has no GenArk
hub as of 2026-09-09.

Since a rebuild is required on either assembly, the assembly to rebuild on is
the better one, and it is not close (NCBI Datasets `dataset_report` for both
accessions, 2026-09-09):

| | GRCz11 | GRCz12ab |
|---|---:|---:|
| scaffolds | 1,917 | 25 |
| contigs | 19,725 | 25 |
| contig N50 | 1.42 Mb | 59.43 Mb |
| gap bases | 4,689,282 | 0 |
| RefSeq annotation | RS_2024_08 (2024-08-15) | RS_2026_07 (2026-07-20) |
| protein-coding genes | 27,158 | 28,415 |

Zebrafish carries the panel's third-longest introns — p99 62,699 bp, max
1,090,140 bp (Table 2) — and GRCz11 spreads 19,725 contigs and 4.69 Mb of N
through them. A 1 Mb intron interrupted by an assembly gap is not the test
case the long-intron axis is there to measure. Keeping GRCz12ab also keeps
the panel's rule uniform: every one of the twenty species is scored against
the current reference annotation of its current assembly.

The third option — GRCz12ab for annotation, GRCz11 for alignment, with a
declared liftover — is technically available. UCSC ships
`danRer11ToGCA_052040795.1.over.chain.gz` (9.4 MB, 2026-05-08), and
`GCA_052040795.1` is the GenBank pair of `GCF_052040795.1`, differing only in
that RefSeq adds chromosome MT (which §4 drops as an organelle anyway). It is
rejected: it adds a chain-projection channel, a second assembly per species in
`panel.tsv` and one more thing for every submission to declare, and what it
buys is the fugu-contaminated alignment that §3.2 already forbids.

Consequence for T-human-008: from this benchmark's point of view zebrafish
belongs with the twelve panel species that have no usable public multiple
alignment, not with the two whose alignment is merely on an older assembly,
so the count is **thirteen of twenty**. Fugu is unaffected — an informant set
at *inference* on a held-out target is permitted (§3.2), and only the
training side of the pair is constrained. `panel.tsv` `notes` and the §2.3
caveat table record the decision.

## 3. Splits and held-out rules

`panel.tsv` assigns every species one of three splits.

- **`train` (10 species).** May be used for training, for hyperparameter
  selection, and as comparative informants. Any development split is carved
  out of these species *by chromosome*, never by borrowing a held-out genome.
- **`heldout` (6 species: chicken, *Ciona*, *Nematostella*, *S. pombe*,
  *P. falciparum*, *Tetrahymena*).** The cross-clade generalization number is
  computed on exactly these. Nothing derived from them may enter training.
- **`heldout_paired` (4 species: human, fugu, honey bee, rice).** A close
  relative *is* in training, on purpose, so the pair isolates a regime shift
  from a clade shift. Reported separately and never averaged into the
  cross-clade number.

### 3.1 The minimum phylogenetic distance rule

> For every `heldout` species, the most recent common ancestor it shares with
> the nearest species in `train` must sit at or above **CLASS** in the NCBI
> Taxonomy.

This is checked mechanically, not asserted:

```
$ python3 benchmark/leakage_check.py
# split separation: MRCA(held-out, nearest training species) must be at or above CLASS
ok	Gallus_gallus	nearest=Mus_musculus,Xenopus_tropicalis	mrca=Sarcopterygii (SUPERCLASS)
ok	Ciona_intestinalis	nearest=Danio_rerio,Mus_musculus,Xenopus_tropicalis	mrca=Chordata (PHYLUM)
ok	Nematostella_vectensis	nearest=Caenorhabditis_elegans,Danio_rerio,Drosophila_melanogaster,+2 more	mrca=Metazoa (KINGDOM)
ok	Schizosaccharomyces_pombe	nearest=Neurospora_crassa,Saccharomyces_cerevisiae	mrca=Ascomycota (PHYLUM)
ok	Plasmodium_falciparum	nearest=Arabidopsis_thaliana,Caenorhabditis_elegans,Danio_rerio,+7 more	mrca=Eukaryota (DOMAIN)
ok	Tetrahymena_thermophila	nearest=Arabidopsis_thaliana,Caenorhabditis_elegans,Danio_rerio,+7 more	mrca=Eukaryota (DOMAIN)
```

`nearest` is every training species tied at the deepest rank, not one of them.
Ties are the normal case and naming only the first row of `panel.tsv` that
reaches the rank is misleading: chicken ties with mouse *and* with
*X. tropicalis*, and the frog is the closer relative in time. At Eukaryota all
ten training species tie, so the list is truncated with a count.

The lineages come from `benchmark/taxonomy.tsv`, a cache of NCBI Taxonomy
lineages fetched on 2026-09-09, so the check is reproducible offline and does
not silently change when NCBI reorganizes a clade.

Rank was chosen over divergence time because a calibrated time tree covering
ciliates, apicomplexans, fungi, plants and animals at one scale does not
exist, and a rule that cannot be evaluated is not a rule. The cost is that
rank is uneven across kingdoms: `Ascomycota (PHYLUM)` separating the two
yeasts is a much older split than `Sarcopterygii (SUPERCLASS)` separating
chicken from mouse. The rule is therefore a floor, not a distance measure,
and §7 lists replacing it with a substitutions-per-site estimate as open.

The four `heldout_paired` species violate the rule deliberately, and the same
script prints them with the relative that puts them in violation:

```
paired	Homo_sapiens	nearest=Mus_musculus	mrca=Euarchontoglires (SUPERORDER)
paired	Takifugu_rubripes	nearest=Danio_rerio	mrca=Neopterygii (SUBCLASS)
paired	Apis_mellifera	nearest=Drosophila_melanogaster	mrca=Pterygota (SUBCLASS)
paired	Oryza_sativa	nearest=Zea_mays	mrca=Poaceae (FAMILY)
```

Each pair is a controlled contrast, which is why they earn their place:

| pair | shared ancestor | what changes |
|---|---|---|
| fugu vs zebrafish | Neopterygii | median intron 170 vs 1,150 bp (6.8×), p99 13,498 vs 62,699 bp, at nearly the same gene count |
| honey bee vs fruit fly | Pterygota | GC 32.5% vs 42%, exons/transcript 10.36 vs 5.96 |
| rice vs maize | Poaceae | genome 386 Mb vs 2.18 Gb, coding fraction 15.8% vs 3.6% |
| human vs mouse | Euarchontoglires | almost nothing — human is here because every published number is reported on it, so it is needed for comparability, and because no pretrained model can honestly claim not to have seen it |

That last row is the point of the split existing: a "held-out human" claim
from any model pretrained on public sequence is not credible, and the split
name is how this benchmark says so out loud.

### 3.2 Leakage channels

Phylogenetic distance is one channel of four, and it is the least dangerous.

1. **Label leakage.** The reference annotation of a held-out species, or any
   projection of it, entering training. This includes liftover and
   projection products — TOGA, CAT, Liftoff — which are annotation, not
   sequence ([10.1126/science.abn3107](https://doi.org/10.1126/science.abn3107),
   [10.1101/gr.233460.117](https://doi.org/10.1101/gr.233460.117)). It also
   includes protein databases: a UniProt or RefSeq protein set containing the
   held-out species' proteins is its annotation in another format, and any
   protein-evidence pipeline must state which release it used and that the
   held-out species was removed from it.
2. **Sequence leakage.** The held-out genome in a pretraining corpus. For any
   model built on a public DNA language model this is the default state, not
   an accident, and it cannot be undone after the fact. The only honest
   handling is declaration: a submission states, per held-out species,
   whether its pretraining corpus contained that assembly or another assembly
   of the same species. Vipsania's own documentation is the cautionary case:
   parsing `docs/test_species.tsv` and `docs/training_species.tsv` from the
   repository HEAD on 2026-09-09 gives 133 unique (model_id, assembly) test
   pairs, of which **75 also appear in the training manifest**
   ([10.64898/2026.08.26.747235](https://doi.org/10.64898/2026.08.26.747235)).
   That is legitimate under an unsupervised protocol — it is sequence
   exposure, not label leakage — and it still has to be said, because a
   reader who assumes otherwise will over-read the reported scores.

   This channel is not hypothetical on this panel. Helixer publishes the
   training and validation species of every shipped checkpoint in
   `docs/model_overview.md`; parsing that file at commit
   [`d17bb49`](https://github.com/usadellab/Helixer/blob/d17bb496b8842542590644b2437798711c0125d5/docs/model_overview.md)
   (109,914 bytes, sha256 `9492017552c6a11d7ede5d72d0188c1202ec94cc2f72463f96181e5d170e2516`,
   read 2026-09-09) and matching the scientific names literally against
   `benchmark/panel.tsv` gives:

   | panel species | Helixer model | listed as | listed accession |
   |---|---|---|---|
   | *Takifugu rubripes* | vertebrate | **training** | GCF_901000725.2_fTakRub1.2 |
   | *Ciona intestinalis* | invertebrate | **training** | GCF_000224145.3_KH |
   | *Drosophila melanogaster* | invertebrate | **training** | GCF_000001215.4_Release_6_plus_ISO1_MT |
   | *Apis mellifera* | invertebrate | **training** | GCF_003254395.2_Amel_HAv3.1 |
   | *Caenorhabditis elegans* | invertebrate | **training** | GCF_000002985.6_WBcel235 |
   | *Arabidopsis thaliana* | land_plant | **training** | TAIR10 |
   | *Oryza sativa* | land_plant | **training** | v7.0 |
   | *Zea mays* | land_plant | **training** | RefGen_V4 |
   | *Saccharomyces cerevisiae* | fungi | **training** | GCF_000146045.2_R64 |
   | *Schizosaccharomyces pombe* | fungi | **training** | GCF_000002945.1_ASM294v2 |
   | *Mus musculus* | mammal | **training** | GCF_000001635.26_GRCm38.p6 |
   | *Mus musculus* | vertebrate | validation | GCF_000001635.27_GRCm39 |
   | *Homo sapiens* | vertebrate | validation | GCF_000001405.39_GRCh38.p13 |
   | *Gallus gallus* | vertebrate | validation | GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b |
   | *Danio rerio* | vertebrate | validation | GCF_000002035.6_GRCz11 |
   | *Xenopus tropicalis* | vertebrate | validation | GCF_000004195.4_UCB_Xtro_10.0 |
   | *Nematostella vectensis* | invertebrate | validation | GCF_000209225.1_ASM20922v1 |

   Sixteen of the twenty panel species are named in a shipped checkpoint's
   lists, eleven of them in a *training* list. The four that appear nowhere
   are *Neurospora crassa*, *Plasmodium falciparum*, *Dictyostelium
   discoideum* and *Tetrahymena thermophila*, and Helixer ships no lineage
   model that covers the three protists, so **on this panel *N. crassa* is
   the only species Helixer can be run on without declaring exposure.**
   Six of the training matches — fugu, *S. pombe*, *A. thaliana*,
   *O. sativa*, *Z. mays*, and mouse under the mammal model — are listed at
   an *older assembly* than the panel's; that is a different assembly of the
   same species and §3.2's rule already calls it exposure. This is a
   documentary membership check against the vendor's own lists, not a
   measurement of how much the exposure is worth: the point is that a
   Helixer number on nineteen of these twenty species is an upper bound
   unless something else is said, and §3.3's `heldout_seen_in_pretraining`
   is where it gets said. Every run in §6 does exactly that.
3. **Alignment leakage.** This is the one specific to comparative methods and
   the one the charter's design directly invites. A whole-genome multiple
   alignment (multiz, Cactus —
   [10.1038/s41586-020-2871-y](https://doi.org/10.1038/s41586-020-2871-y))
   built to include a held-out species leaks that species' *sequence* into
   every training window that overlaps it, and if the alignment was filtered,
   scored, or projected using annotation, it leaks labels too. Rules:
   - **Jointly inferred alignments** — Cactus, Ensembl EPO, anything whose
     columns and ancestral sequences were estimated over the whole taxon set
     at once — must be *rebuilt* with the held-out species removed. A row
     dropped after the fact still shaped the columns it was dropped from.
   - **Reference-anchored alignments** — multiz, which stacks pairwise
     to-reference alignments in guide-tree order — may instead have the row
     dropped at cut time, because dropping a row there is a filter and not a
     realignment. The run declares which rows were dropped, in
     `alignment_rows_dropped` (§3.3). This matters for exactly two panel
     alignments: dm6 multiz124way contains honey bee (`apiMel4`,
     `heldout_paired`) and mm39 multiz35way contains human. Rebuilding either
     is a cluster job; dropping the row is free. The distinction is marx's,
     from `docs/data-sources.md` §7, and a submission that cannot say which
     kind its alignment is must assume the stricter rule.
   - An informant set used at *inference* on a held-out target is permitted
     and must be declared: `benchmark/leakage_check.py --informants` takes a
     `target<TAB>informant` TSV and reports each informant's MRCA with the
     target, marking any closer than CLASS as `declared-close`.
   - Informant *annotation* may never be projected onto a held-out target.
     This is the boundary between "comparative prediction" and "annotation
     transfer", and TOGA-class methods sit on the far side of it by design.
   - Conservation tracks (phyloP, phastCons) computed from an alignment that
     included the held-out species carry the same leak as the alignment.
4. **Evidence leakage.** RNA-seq from the held-out species is evidence, not
   annotation, and evidence-based pipelines cannot run without it. It is
   permitted at inference, must be declared with accessions, and its
   assembly cost (STAR, StringTie —
   [10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122)) counts against the
   cost budget in §4.7. It may never be used to fit model parameters.

### 3.3 What a submission declares

A run that is not accompanied by this block is not comparable and is not
scored:

```yaml
model: <name and version or commit>
training_species: [...]            # every genome whose sequence or annotation was used
pretraining_corpus: <name/version> # or "none"
heldout_seen_in_pretraining: {Gallus_gallus: no, Homo_sapiens: yes, ...}
protein_db: <release>              # or "none"; state held-out removal
alignment: <how built, which taxa> # or "none"; say whether it is jointly
                                   # inferred (Cactus, EPO) or
                                   # reference-anchored (multiz)
alignment_rows_dropped: [...]      # rows removed at cut time, per §3.2; or "none"
informants: <path to target/informant TSV>  # or "none"
rnaseq: {species: [accessions]}    # or "none"
hardware: <CPU model, cores, RAM, GPU model, VRAM>
```

## 4. Metrics

All metrics are computed against the reference annotation of the same
assembly, on the **nuclear primary assembly only** (no alt loci, no patch
scaffolds, no organelle genomes, no unplaced scaffolds under 10 kb), and are
reported per species. Sensitivity
is TP/(TP+FN) and precision is TP/(TP+FP) throughout. The word *specificity*
is deliberately not used anywhere in this benchmark: in the gene prediction
literature "SP" denotes precision in some papers and the true-negative rate
in others, and reusing the word guarantees that two correct implementations
of this spec will disagree. F1 is the harmonic mean of sensitivity and
precision.

The scored sequence set is derived from the reference and reported in
`sequence_selection` of every result, with the dropped sequences grouped by
reason, so no run has to be taken on trust. Where the reference carries RefSeq
`region` features the rule is exact: drop `genome=mitochondrion`,
`genome=chloroplast` and the other organelle values, and drop the
`genome=genomic` regions that carry a cytogenetic `map=` band, which is what
an alt locus or a patch scaffold is. **The discriminator is `map=`, not
`chromosome=`.** Unlocalized scaffolds carry `map=unlocalized` and unplaced
ones carry `chromosome=Unknown` with no `map=` at all: in GRCh38.p14 every one
of the 680 `genome=genomic` regions has a `chromosome=`, as do all 675 of
maize's and all 32 of *Nematostella*'s, so a rule keyed on `chromosome=`
discards every unplaced scaffold in the panel. Where the reference has no
`region` features — Ensembl output, and most predictor output — organelles and
alt loci are recognised from the sequence name instead (`MT`, `MtDNA`, `Pt`,
`chrM`, `*_PATCH`, `*_ALT`); that is a heuristic,
`reference_has_region_features` in the result says which path was taken, and
`--seqids` overrides both.

Organelles are dropped because they are out of the charter's scope and because
they do not use the nuclear genetic code: human chrM is translation table 2,
so every start and stop codon check on it would be wrong. Before this rule
existed the scorer silently included the mitochondrion of human, yeast,
*Arabidopsis* and maize, and the *Arabidopsis* and maize chloroplasts.

Predictions land on sequences that survive this filter or they are not scored.
The result reports `predicted_transcripts_not_scored` and
`predicted_sequences_absent_from_reference`, and the scorer warns on stderr
when more than half the prediction falls outside the scored set, because a
prediction in the wrong naming convention otherwise scores 0.0 everywhere and
is indistinguishable from a bad predictor: Ensembl calls *C. elegans*
chromosome I `I` and RefSeq calls it `NC_003279.8`.

**A CDS row in the reference is not automatically a protein-coding gene.**
Two classes are excluded from truth, counted by reason in
`reference_transcript_selection`, and excluded from the prediction by the same
test so that an identity run still scores exactly 1.0:

- **Pseudogenes.** RefSeq annotates them with real CDS blocks carrying
  `pseudo=true`, under a `gene_biotype=pseudogene` parent: 32 transcripts in
  *S. pombe* — one of them described as "malic enzyme with 2 frameshifts" —
  and 314 in GRCh38.p14, of which 198 sit on a sequence this benchmark scores. A predictor that correctly declines to call a
  frameshifted pseudogene would otherwise be charged a false negative at
  nucleotide, exon, locus and transcript level. The `MIN_INTRON` rule (§4.3)
  already kept their frameshift gaps out of the splice metrics; this is the
  other half.
- **Gene fragments.** Immunoglobulin and T-cell receptor segments
  (`gene_biotype=V_segment`, `D_segment`, `J_segment`, `C_region`: 891
  transcripts in GRCh38.p14, 390 of them on a scored sequence) are pieces of a
  gene assembled somatically, with no start or stop codon of their own. The panel's other 18 species have
  none.

Both are recognised only from what the file states, so a prediction that
carries neither attribute loses nothing, and `--score-all-transcripts`
restores the old behaviour for auditing. The biotype is read from
`gene_biotype` (RefSeq), `gene_type` (GENCODE) or `biotype` (Ensembl),
whichever the gene row carries; reading only the first applies the filter to a
RefSeq-shaped submission and not to a GENCODE-shaped one, which charges the
latter false positives for exactly the rows the reference is forbidden to
score (§6). On the human reference the filter
drops 588 of 132,030 transcripts and the identity run stays at 1.0 on every
metric.

Predictions are compared on the CDS, not the transcript, unless a metric
says otherwise; UTRs are out of scope for this charter and a model that does
not predict them must not be penalized.

### 4.1 Nucleotide level

A genomic position is a positive if it lies in the CDS of any annotated
protein-coding transcript on that strand. Strand matters: a position coding
on the plus strand and predicted coding on the minus strand is a false
positive and a false negative, not a hit. Overlapping genes are unioned.
Report sensitivity, precision, F1, and the Matthews correlation coefficient;
MCC is reported because the positive class is 3.6% of the genome in maize and
73.1% in *S. cerevisiae*, and F1 is not comparable across that range.

### 4.2 Exon level

A predicted CDS exon is a true positive only if **both** boundaries match an
annotated CDS exon of some transcript exactly. Report sensitivity, precision
and F1, then stratify by exon type — initial (contains the start codon),
internal, terminal (contains the stop), single (the whole CDS) — because
these have different error profiles and a pooled number hides it. Also report
the fraction of predicted exons that overlap an annotated exon but match
neither boundary; that is the "roughly right, structurally wrong" class that
nucleotide metrics forgive and downstream users cannot.

### 4.3 Splice-site level

Each annotated intron contributes a donor and an acceptor position. A
predicted site is a true positive if it is at the exact position and strand.
Report separately for donors and acceptors, and stratify by:

- **canonical vs non-canonical** dinucleotides (GT-AG, GC-AG, AT-AC, other),
  **per intron rather than per site** — see below;
- **intron length decile**, computed by the scorer from the reference
  annotation of the species being scored rather than read from `panel.tsv`,
  so that the long-intron tail is visible instead of averaged away. The cuts
  are emitted with the result (`splice.intron_length_decile_cuts`);
  `panel.tsv` carries only p10/median/p90/p99, which is enough to justify the
  panel and not enough to bin against.
- **local GC** in a 200 bp window, in five bins, centred on the donor and on
  the acceptor separately.

The intron-length stratification is the single most informative panel in the
whole benchmark for the charter's question, because it is where clade-specific
models are expected to differ from a species-independent one.

**The three stratifications do not share a denominator, and the result says
so** (`splice.strata_units`). The donor and acceptor totals, the
intron-length deciles and the local-GC bins are all counted per *site*; the
dinucleotide table is counted per *intron*, and has to be. A dinucleotide
class is the pair (donor, acceptor), so it is a property of the intron: a
donor shared by two introns that end `AG` and `AC` is in two classes at once.
The gap is not cosmetic on a real genome — the human reference's 218,446
introns carry 188,913 distinct donors and 193,359 distinct acceptors, so
`by_dinucleotide` sums to 15.6% more than the donor total printed beside it,
and a reader who adds the row up expecting the donor denominator gets a
different number. `reference_sites_multiple_dinuc_classes` and its predicted
counterpart report how many sites are in more than one class, which is the
count that would have to be resolved before a per-site table could exist. The
local-GC table is given for donors *and* acceptors; before the human run it
existed for donors only, which quietly made one of the three strata
single-sided while this section asked for both.

One position can belong to introns of two different lengths, and then the
decile is a choice rather than a fact: alternative splicing shares a donor
between a short and a long intron. That is 0 of 281 *S. cerevisiae* donors,
73 of 15,845 in *N. crassa*, 8,086 of 220,976 (3.7%) in *T. rubripes* and
23,387 of 188,913 (**12.4%**) in human, where the shortest and the longest
intron at one site are over 1 Mb apart. The site is scored once — the totals
never depended on this — and stratified by its **shortest** intron, which is
arbitrary but fixed; the number of sites the rule was applied to is reported
as `splice.donor.reference_sites_multiple_intron_lengths` and its acceptor and
prediction counterparts, so a reader can see how much of the column rests on
it. Before this rule the length kept was whichever intron the internal set
happened to yield last, which made the whole stratification depend on the
interpreter's hash seed: scoring the same two files twice moved per-decile
counts by tens of sites while every total stayed identical.

A gap between consecutive CDS blocks is only treated as an intron if it is at
least 20 bp. This is not a tuning knob: RefSeq encodes a programmed ribosomal
frameshift as two CDS blocks separated by 1 bp, and 47 of the 343 CDS gaps in
the *S. cerevisiae* reference are Ty retrotransposon frameshifts rather than
splice junctions. Scoring them as splice sites would have put a 14% floor of
non-splice-sites into the yeast donor and acceptor counts. Sub-threshold gaps
are counted and reported (`splice.reference_cds_gaps_below_min`), not dropped
silently.

### 4.4 Gene and transcript level

Two numbers, because they answer different questions.

- **Transcript exact match.** A predicted transcript is a true positive if
  its full ordered CDS exon chain, including start and stop, is identical to
  an annotated transcript's. This is the strict number and it will be low;
  low is fine as long as it is comparable. Exact-gene agreement has been the
  field's strictest endpoint since TWINSCAN
  ([10.1093/bioinformatics/17.suppl_1.s140](https://doi.org/10.1093/bioinformatics/17.suppl_1.s140)),
  and the reference sets behind published values differ enough that only a
  shared panel makes them comparable.
- **Locus level.** A predicted locus is a true positive if it overlaps an
  annotated protein-coding gene on the same strand by at least 1 CDS
  nucleotide and no other predicted locus overlaps it better. This detects
  the two structural errors exact matching cannot distinguish from a plain
  miss: **fusion** (one prediction spanning ≥2 annotated genes) and **split**
  (≥2 predictions inside one annotated gene). Report the fusion and split
  counts explicitly. Overlap in the *annotation* is charged to the annotation
  and not to the predictor: build the graph whose nodes are annotated genes
  and whose edges join genes whose CDS overlap on the same strand, and count a
  prediction as a fusion only when the annotated genes it touches lie in two
  or more different **connected components** of that graph. Splits are the
  mirror image, over the prediction's own overlap graph. A pairwise "do these
  two genes overlap each other" test is not enough, because overlap chains:
  GRCh38.p14 has 82 places where gene B overlaps both A and C while A and C
  are disjoint, and under the pairwise rule a perfect prediction of the human
  reference scores 82 fusions and 82 splits (and 137 and 132 on
  *S. cerevisiae*, which has 91 same-strand CDS-overlapping gene pairs).
  Under the component rule an identity run scores exactly zero of both on
  every panel species tested. The reference's own overlapping-pair count is
  reported alongside (`locus.reference_overlapping_locus_pairs`), because it
  is the ceiling on how well any one-to-one locus matching can do.
  Fusions are the failure mode manual curation of
  *P. pacificus* found most of ([10.64898/2026.02.18.706511](https://doi.org/10.64898/2026.02.18.706511)).

Matching is **greedy, not optimal**, at both levels, and the tie-break is part
of the specification so that a second implementation of this document agrees
with `score.py` on the same input: candidate pairs are taken in descending
order of shared CDS bases, ties broken by reference id then predicted id
ascending, and a pair is kept only if neither side is already used. Optimal
bipartite matching would score marginally differently in the rare case where
a locally best pair blocks two better ones; greedy is cheap, deterministic and
what is implemented.

Isoform rule: a predictor emitting one transcript per locus is scored against
the **single best-matching annotated isoform** per locus, and the remaining
annotated isoforms are not counted as false negatives. A predictor emitting
several isoforms is matched the same greedy way within the locus, and
unmatched predictions are false positives. Without this rule,
species with deep isoform annotation (human, zebrafish) are systematically
penalized against species with one transcript per gene (*P. falciparum*).

Within a locus, **an exact chain match outranks shared bases**, and this is
not the same as breaking ties on shared bases. An annotated isoform that
*contains* the prediction's whole CDS shares exactly as many bases with it as
the isoform the prediction *equals*, so on shared bases alone the two are
tied and the accession that sorts first wins — which is not, in general, the
one that matches. *T. rubripes* `rab44` is the case: the Helixer prediction is
`rna-XM_029826178.1` base for base, `rna-XM_011613896.2` is that chain plus 51
bases at one end, both share 10,749, and `011` sorts before `029`. Under the
old order 343 of 5,907 exact fugu matches (5.8%) and 34 of 6,923 *N. crassa*
ones were scored as misses; fugu's transcript F1 was 0.238 and is 0.252
(§6). Exactness-first is also the only order that makes the metric a property
of the two annotations rather than of their accession strings. Pairs are
therefore taken in ascending order of (0 if the chains are identical else 1,
then descending shared CDS bases, then reference id, then predicted id).

The benchmark deliberately does **not** designate one representative isoform
per locus as truth. Doing so would pick a winner the reference does not pick,
and the isoform rule already removes the penalty deep annotation would
otherwise carry. A model that needs a single target per locus — a
per-base labelling of a genome window, for instance — is making a *training*
decision, not a scoring one, and it is free to take the union over isoforms,
the longest CDS, or MANE Select where it exists; the benchmark scores whatever
it emits against the best-matching annotated isoform either way. What such a
model cannot do is claim the transcript-level column measures isoform choice:
against a one-prediction-per-locus model, that column is bounded by how close
the chosen target is to *some* annotated isoform, and fugu's 0.252 against
*S. cerevisiae*'s 0.860 is mostly that bound moving, not the model.

### 4.5 Start and stop codons

Reported separately from exon boundaries, because they are the part of the
problem where the genetic code is a prior rather than a signal: exact-position
sensitivity and precision for start codons and for stop codons. *Tetrahymena*
is the diagnostic — under genetic code 6 a model that has hard-coded TAA/TAG
as terminators will show near-zero stop-codon precision there and normal
numbers everywhere else, which is a signature no aggregate score would show.

**An incomplete CDS end contributes no codon.** RefSeq marks one with
`start_range=` or `end_range=` on the CDS row — which biological end that is
depends on the strand — and 1,415 human, 194 maize and 6 *S. pombe*
transcripts are partial this way. A 5'-partial gene has no annotated start
codon to hit, so counting its first base as a reference start codon charges
every predictor a false negative it cannot avoid. The two ends are excluded
independently, and dropping the reference position alone would only move the
charge, so a predicted position inside the span of a reference transcript
partial at that end is dropped from the denominator too — unless it coincides
with a surviving reference position, which is a real match and stays a true
positive. The counts are reported as `reference_partial_5prime` and
`reference_partial_3prime`. Note that this reads the **CDS** row, not the
mRNA: 788 of *S. pombe*'s 5,166 mRNA rows carry `partial=true`, describing
incomplete UTRs, while only 6 of its CDS rows do.

**A prediction states its own incomplete ends differently, and by omission.**
`start_range=`/`end_range=` is a *reference* convention; no predictor writes
it, so a scorer that reads only those attributes leaves
`predicted_partial_5prime`/`_3prime` structurally at zero — which is what they
read on the first submission that actually contained partial genes (AUGUSTUS
`--genemodel=partial` on *T. thermophila*: 250 transcripts with no start
codon, 69 with no stop). What a GTF-lineage predictor does instead is emit no
`start_codon` or `stop_codon` feature for the end that ran off the contig, so
the omission *is* the statement, and it is legible only in a file that uses
those features at all — Helixer and Tiberius emit neither and must not be read
as wholly partial. `score.py` infers a predicted partial end that way when the
file uses codon features, drops that end from the codon denominators exactly
as it does for the reference, and records which convention the two counts came
from in `predicted_partial_source` (`range_attribute`,
`missing_codon_feature` or `none`). Feature omission is the better of the two
available signals: the first CDS block's `phase` is the other, and it misses
the 69 of those 250 truncated *T. thermophila* genes that happen to resume in
frame with phase 0. Left unread, each truncated end is charged as a wrong
codon the predictor never claimed: 213 of those chains and 30 survive §4
selection onto a scored sequence, and reading them removes **208
false-positive starts and 22 false-positive stops** from that run, with no
true positive touched at either end.

**Worse, an unread 3'-partial end can be handed a stop codon it never
predicted.** The bare-GTF case and the truncated-gene case both arrive at the
3 bp extension as a chain with no `stop_codon` feature, and they need opposite
treatment. Extending a truncated gene does not merely misplace 3 bp: the
extended position can land on the reference's real stop and be scored a
**true positive** on a gene the predictor never finished. In the fixture for
this the pre-fix scorer returns stop `tp` 3 of 3 where 2 of 3 is right.
`score.py` skips the extension for any chain the file declares 3'-partial and
counts the skips in `transcripts_stop_not_extended_partial`; a file with no
`stop_codon` feature anywhere declares nothing, so every chain is still
extended there. The *T. thermophila* run cannot exercise this, because the
`tetrahymena` parameter set puts the stop inside the CDS and the extension
never runs; the *A. mellifera* run in §6 can, because `honeybee1` puts it
outside, and it is the first real run where the skip fires. 13 chains are
skipped there: 8 false-positive stop codons and 2 false-positive terminal
exons go away, and 27 of the 39 invented bases were false-positive nucleotides
(the other 12 include 3 that had been scored as true positives, which is the
same mechanism as the fixture's spurious stop, one order of magnitude
smaller). Reading the 5' omissions on the same run removes 6 false-positive
starts.

**The stop codon is inside the CDS here, and half the field disagrees.** GFF3
says a CDS chain includes its stop codon, and every reference in Table 3 does.
GTF says the opposite, and every tool with a GTF lineage inherits it:
AUGUSTUS emits the stop codon as a separate `stop_codon` feature, and BRAKER,
GeneMark and SNAP are in the same family. It is worth being exact about
AUGUSTUS, because "AUGUSTUS excludes the stop codon" is not a property of
AUGUSTUS: `stopCodonExcludedFromCDS` is a line in each *species* parameter
file, and of the 166 species shipped with AUGUSTUS 3.5.0, **44 set it `true`
and 122 set it `false`**. The two the benchmark happened to run first,
`saccharomyces_cerevisiae_S288C` and `schizosaccharomyces_pombe`, are both in
the 44, together with `ciona`, `honeybee1`, `neurospora_crassa` and
`zebrafish`; `human`, `fly`, `chicken`, `caenorhabditis`, `arabidopsis`,
`maize`, `rice`, `nematostella_vectensis`, `toxoplasma` and `tetrahymena` are
in the 122. The census is one line, and it is worth running before assuming a
convention:

```
grep -h stopCodonExcludedFromCDS "$AUGUSTUS_CONFIG_PATH"/species/*/*_parameters.cfg | awk '{print $2}' | sort | uniq -c
```

Two runs of the same
binary at the same version, one species apart, therefore need opposite flags,
which is the strongest argument for the rule below: there is no tool-level
answer to fall back on. Three base pairs is not a rounding
error, because those three sit at the 3' end of *every* CDS chain: they are
every terminal exon, every single-exon gene, every stop codon and every exact
transcript match. Scoring an AUGUSTUS prediction of *S. cerevisiae* against
the RefSeq annotation without correcting for this gives **exon F1 0.03 and
stop-codon F1 0.00 while nucleotide F1 stays at 0.96** — a result that reads
like a real one, not like a bug.

So the convention is **detected from the prediction, not taken on trust**.
`score.py` compares each transcript's `stop_codon` features against its CDS
blocks and reports `stop_codon_convention_detected` as `inside`, `outside`,
`mixed` or `unknown` in every result, and warns on stderr when the detection
disagrees with the `--stop-outside-cds` flag or when the file mixes the two.
With the flag, the prediction is brought into the reference's convention
before anything is scored: where a `stop_codon` feature exists it is unioned
into the chain, which is right even when the stop is split across an intron
where none exists — a bare GTF-derived CDS set — the last block in the
direction of translation is extended by 3 bp instead and the guess is counted
in `transcripts_stop_extended_by_3bp` so it is visible. A gene truncated by
the end of a contig has no `stop_codon` feature either and must *not* be
extended; the two absences are told apart below. A submission states
which convention its output uses.

**Feature-based detection answers nothing for a tool that emits no
`stop_codon` feature, and the modern ones do not.** Helixer's GFF3 is
`gene`/`mRNA`/`exon`/`CDS`/`five_prime_UTR`/`three_prime_UTR` and nothing
else; Tiberius is the same. For those files `stop_codon_convention_detected`
is `unknown` and, before this was fixed, the convention was in practice *the
default of a command-line flag* — the one thing this section says it must not
be. With `--genome` it is decidable without any feature: for each predicted
chain, read the last codon in the direction of translation and the codon
immediately after it, and see which one is a stop. `score.py` now does that
and reports `stop_codon_convention_from_genome`, the three counts it rests on
(`stop_convention_genome_inside`, `_outside`, `_neither`), and
`stop_codon_convention_source`, which is `stop_codon_feature`, `genome`,
`flag`, or `assumed`. **`assumed` is the value to look for before believing a
terminal-exon, single-exon, stop-codon or exact-transcript score.** When the
prediction has no `stop_codon` features and the genome says `outside`, the
3 bp are put back automatically; evidence outranks a default, but never an
explicit `--stop-outside-cds`, and the note goes to stderr either way. This
correction is deliberately *not* extended to a prediction that does carry
`stop_codon` features: there the file states its own convention, disagreeing
with the flag is a user error that should be seen rather than papered over,
and the loud stderr warning above already covers it. Only the no-feature case
has nothing to disagree with and no evidence but the default. Chains
whose terminal CDS block is under 3 bp are skipped rather than walked back
across the intron, so a split stop codon is never counted as evidence.

The stop codons used are those of the species' NCBI translation table, taken
from `panel.tsv`'s `genetic_code` and passed as `--genetic-code`. This is not
decoration: under table 6 *Tetrahymena* reads TAA and TAG as glutamine, so a
scorer that hard-coded table 1 would find "stops" in the middle of two thirds
of its genes and could talk itself into extending every chain by 3 bp. Under
table 6 the check returns `unknown` on evidence that table 1 would call
decisive, and declines to guess. Only the tables the panel uses are
implemented; an unlisted one is an error rather than a silent fall back to
table 1.

Measured on the two Helixer runs of §6: 24,166 of 24,176 *T. rubripes* chains
and 10,303 of 10,303 *N. crassa* chains end in a stop codon and none is
followed by one (the other 10 fugu chains end in neither, which is what a
chain running into an assembly gap looks like), so Helixer includes the stop
in its CDS — the same convention as the references. That was previously the
default's assumption and is now the run's measurement.

The window plan for this check is made before the convention is known, and
merging a `stop_codon` feature into a chain can add a junction the plan never
saw: AUGUSTUS splits a stop codon across an intron, and one predicted intron
of the *S. pombe* cross-parameter run is exactly that. Its dinucleotide came
back `unknown` — the same silent class the §6 window-fetcher defect produced —
until the plan was extended after the merge with a second pass over the FASTA
for the windows the first pass did not ask for. The second pass runs only when
a merge changed a chain, so the usual case still reads the genome once.

### 4.6 Proteome completeness

Translate the predicted CDS and run BUSCO in protein mode
([10.1093/molbev/msab199](https://doi.org/10.1093/molbev/msab199)) against
the deepest lineage dataset that applies to the species, reporting complete,
duplicated, fragmented and missing. Run it on the *reference* proteome too
and report both: BUSCO scores a proteome against an expectation, not against
this panel's truth, so the reference value is the ceiling and the number that
matters is the gap. OMArk
([10.1038/s41587-024-02147-w](https://doi.org/10.1038/s41587-024-02147-w)) is
reported alongside where an OMA clade exists, because it also flags
contamination and fragmentation that BUSCO's completeness number hides.

### 4.7 Cost

Measured, not quoted, on one declared machine for the whole panel:

- **wall clock seconds per Mb** of assembly, end to end, including any
  evidence preprocessing (alignment, RNA-seq assembly, protein search) but
  excluding one-time model download;
- **CPU-seconds per Mb** (user+sys across all processes), which is the number
  that differs from wall clock by the parallelism, and the one that made the
  EGAPx figures interpretable: 71 CPU-hours / 3 wall-hours for 144 Mb of
  *Drosophila*, 425 CPU-hours / 5.5 wall-hours for 1.1 Gb of chicken;
- **peak resident set** in GB, and peak GPU memory in GB where used;
- **first-run penalty** in seconds, reported separately, because a pinned
  TensorFlow that must JIT from PTX on a current consumer GPU can cost half
  an hour before the first prediction;
- **failure rate**: runs that did not produce output, over runs attempted,
  per species. A tool that fails on a third of the panel has a cost that no
  successful-run timing describes.

Cost is reported per species and never averaged across the panel, since
genome size varies by 260×. T-human-009 owns the baseline measurements; this
document owns their definition.

### 4.8 What gets reported, and the aggregate

The primary result is a table: 20 species × {nucleotide F1, exon F1, splice
donor F1, splice acceptor F1, transcript exact F1, locus F1, fusion count,
split count, BUSCO complete %, wall s/Mb, peak GB}. Everything else is a
view on it.

Two aggregates, both unweighted means so a large genome cannot dominate:

- **Cross-clade score**: mean over the 6 `heldout` species, excluding
  *Tetrahymena* (reported, never ranked, §2.3). This is the number the
  charter's generalization claim is judged on.
- **Regime-shift score**: mean over the 4 `heldout_paired` species, reported
  next to the corresponding training relative's score. The interesting
  quantity is the *drop*, not the level.

A submission that reports only an aggregate is incomplete. The stratified
panels — exon type, intron-length decile, GC bin, canonical vs non-canonical
splice site — are the deliverable; the aggregate is the headline.

## 5. Files

| file | what it is |
|---|---|
| `benchmark/panel.tsv` | the manifest: 20 species, accessions, annotation release, MD5 of the annotation, all statistics in Tables 1–2, split assignment |
| `benchmark/taxonomy.tsv` | NCBI Taxonomy lineage per species, cached 2026-09-09, so the distance rule is reproducible offline |
| `benchmark/fetch.py` | downloads genomes/annotations from the NCBI FTP mirror into a directory outside the repository and verifies both NCBI's MD5 and the manifest's; stdlib only |
| `benchmark/annotation_stats.py` | recomputes every statistic in Tables 1–2 from a GFF3; stdlib only |
| `benchmark/leakage_check.py` | enforces §3.1 and reports §3.2 informant separation; exits non-zero on violation |
| `benchmark/score.py` | scores one predicted GFF3 against the reference and emits the §4 metrics as JSON; refuses to run without the §3.3 declaration; `--self-test` checks it against built-in fixtures |
| `benchmark/report.py` | joins the JSON rows to `panel.tsv` and emits the §4.8 table and the two aggregates |

No genome or annotation data is committed. `panel.tsv` records the MD5 of
each annotation file so a fetch is verifiable years later, and `ftp_dir`
records the resolved NCBI directory because the FTP directory name is the
submitter's assembly name, which is not always the name the Datasets API
returns (`GCF_034140825.1` is `AGIS1.0` in the API and `ASM3414082v1` on the
FTP site).

## 6. Scoring implementation

`benchmark/score.py`, standard library only, takes a reference GFF3, a
predicted GFF3 and a species name and emits one JSON object with every
annotation-derived field in §4.8: §4.1 nucleotide with MCC, §4.2 exon with the
four-way type stratification and the overlap-but-no-boundary class, §4.3
donors and acceptors with the intron-length-decile stratification (and, when
`--genome` is given, the dinucleotide-class and local-GC ones), §4.4
transcript exact match under the isoform rule and locus matching with fusion
and split counts, and §4.5 start and stop codons. It refuses to run without a
§3.3 declaration carrying every required key, and records that file's SHA-256
in the result. §4.6 (BUSCO/OMArk) and §4.7 (cost) come from other tools;
`benchmark/report.py` joins them in through `--cost` and emits the §4.8 table
and the two aggregates. It refuses to print an aggregate over an incomplete
panel unless `--partial` says so, because an unweighted mean over a subset is
a different number wearing the same name.

Four decisions the implementation forced, all above: the intron-length
deciles are computed from the reference at score time (§4.3), a CDS gap under
20 bp is not a splice junction (§4.3), the stop-codon convention is detected
from the prediction rather than trusted to a flag (§4.5), and a reference CDS
row is not truth until it has been checked for `pseudo=true`, a
non-`protein_coding` biotype and an incomplete end (§4, §4.5).

Verification so far:

- `python3 benchmark/score.py --self-test` scores built-in fixture pairs and
  checks 145 expected values — it prints the count it ran, so this sentence
  cannot drift from the code again; it said 116 while the code ran 104: a
  prediction with one exact transcript, one
  shifted minus-strand boundary, one overlapping-but-unaligned locus and one
  spurious locus; a fusion-and-split pair; a self-comparison that must score
  exactly 1.0 on every metric; two sequence-selection fixtures, one in RefSeq
  shape and one in Ensembl shape, whose every line is a sequence shape taken
  from a real panel reference (chromosome, alt locus, unlocalized scaffold,
  unplaced scaffold, bare scaffold, mitochondrion, chloroplast, sub-10 kb
  scaffold); the streaming FASTA window reader against a plain read; and a
  two-record FASTA in which a window overruns the end of a record, on the
  first record and on the last, checking that it is clipped and that the
  windows queued behind it are still served; an AUGUSTUS-shaped prediction
  (CDS 3 bp short, separate `stop_codon` features, on both strands) which must
  score exactly 1.0 everywhere with `--stop-outside-cds`, be *detected* as
  `outside` with or without the flag, and lose exactly its terminal exons,
  single-exon genes, stop codons and transcript matches without it; the same
  file with the `stop_codon` rows stripped, where the convention is
  undetectable from the file and the 3 bp are guessed back; that same stripped
  file against a genome fixture whose only stop codons sit at the three
  positions the reference chains end on, where the convention has to be read
  off the sequence, the 3 bp are put back *without* the flag and the run
  recovers exactly 1.0, while a prediction that already includes its stop
  codons is read as `inside` and left alone, and the same file under
  translation table 6 is called `unknown` rather than extended; a two-isoform
  gene sharing one donor between a 100 bp and a 6,900 bp intron, which must
  land in the decile of the shorter one and be reported as ambiguous (§4.3),
  and that same fixture with a genome in which its two introns end `AG` and
  `AC`, so the one donor is in two dinucleotide classes at once and the
  intron-level table must count 2 where the site-level total counts 1, while
  the local-GC table must exist for both sides and count 1 donor against 2
  acceptors (§4.3);
  a prediction whose `stop_codon` feature sits across an intron from its last
  CDS block, whose merged-in junction must still get a dinucleotide class
  (§4.5); a two-isoform gene in which one isoform is the other plus 51 bases
  at one end and the prediction equals the shorter, where the containing
  isoform must not take the pairing and the same file against itself must
  still be 1.0 (§4.4); a prediction that declares an immunoglobulin segment
  and a pseudogene the way GENCODE and Ensembl declare them, which must cost
  it nothing and must not be reported as an off-panel sequence; a prediction
  that reuses transcript ids across two sequences, which must score exactly
  what the same prediction scores after the ids are made unique, block for
  block; and a `--genemodel=partial`-shaped prediction with one gene truncated
  at a contig start and one at a contig end, where the omitted codon features
  must be read as the statement of partiality, the truncated 3' end must not
  be extended, and the pre-fix code returns stop `tp` 3 of 3 where 2 is
  right (§4.5).
  The self-test is run under several `PYTHONHASHSEED` values, because two of
  its checks exist to catch results that depended on it.
- **Identity runs re-verified after the fixes above** on
  *S. cerevisiae*, *N. crassa*, *T. rubripes*, *C. elegans* and human: F1 1.0
  and MCC 1.0 on every metric with 0 fusions and 0 splits, and 0, 73, 8,086,
  2,120 and 23,387 donors respectively whose decile came from the
  shortest-intron rule. Re-verified again after the id-namespacing and
  predicted-partial changes, with *T. thermophila* and *A. mellifera* added:
  all seven still 1.0 and MCC 1.0 everywhere with 0 fusions and 0 splits.
- **Identity runs on all twenty panel references.** Every species in
  `panel.tsv` scored against its own reference gives F1 = 1.0 and MCC = 1.0 on
  every metric, with **0 fusions and 0 splits** everywhere. Cost on one laptop
  core ranges from 1.0 s / 43 MB (*P. falciparum*, 23 Mb) to 52 s / 335 MB
  (*Z. mays*, 2.18 Gb); human is 38 s / 0.92 GB over 3.10 Gb. The filter keeps
  3 sequences (*S. pombe*) to 685 (*Z. mays*); human is 102 scored sequences
  and 3,101,538,863 bp (24 chromosomes plus 78 scaffolds; 511 alt loci and
  patches, 91 sub-10 kb scaffolds and the mitochondrion dropped), 131,442
  transcripts of 132,030 (198 pseudogene and 390 gene-fragment transcripts
  dropped, §4), 19,932 loci. Every one of the twenty took the RefSeq `region`
  path (`reference_has_region_features: true`).
- **`--genome` across seven species and at 3 Gb.** The §4.3 dinucleotide and
  local-GC strata run on *T. rubripes* (384 Mb, 229,939 reference introns,
  32 s and 0.69 GB), *D. melanogaster*, *C. elegans*, *S. pombe*,
  *P. falciparum*, *S. cerevisiae* and now **human**, which is the case the
  streaming window reader was written for and had never been given: the
  identity run reads 3,101,538,863 scored bp across 102 sequences out of a
  973 MB gzip in 67 s and 1.35 GB, and the GENCODE 50 run — 347,110 predicted
  introns against the reference's 218,446 — in 88 s and 2.12 GB. Not one
  window came back unserved: no `unknown` dinucleotide class and no dropped
  GC bin in either run, over the 705 records the FASTA carries. Every metric
  in the GENCODE result is identical to the run without `--genome`; what the
  genome adds is the two strata and a second, independent reading of the
  stop-codon convention, which agrees with the file's own `stop_codon`
  features (349,754 chains `inside`, 579 `outside`, 18,253 neither).
  The human reference's own splice census: 215,956 GT-AG, 1,939 GC-AG, 219
  AT-AC and 332 other, so **1.14% non-GT-AG**. Non-GT-AG donor fractions recovered
  from the genome: *T. rubripes* 1.57%, *D. melanogaster* 0.98%,
  *C. elegans* 0.97%, *S. cerevisiae* 3.91%, *P. falciparum* 0.18%,
  *S. pombe* 0.16%.
  **This is what found the window-fetcher defect below**, and it is the
  argument for reporting the strata at all: one *T. rubripes* intron ending
  52 bp from the end of a 43 kb scaffold was being counted in an `unknown`
  dinucleotide class rather than as the GT-AG it is.
- **Defect found and fixed by that run: the window fetcher silently dropped
  windows near the end of a sequence, and everything queued behind them.**
  `WindowFetcher` serves a sorted queue of windows and pops one when its *end*
  has been read, but the queue is ordered by window *start*. A local-GC window
  centred on a splice site within `GC_WINDOW // 2` (100 bp) of the end of a
  sequence can never satisfy that test, so it stalled at the head of the queue
  and hid every later window on the same sequence. The fetcher returned `None`
  for all of them, which `dinuc_class` reports as `unknown` and `gc_bin` drops
  from the denominator — a silent loss, not an error. The fix serves the
  remaining queue against what was read when a record ends, clipping the
  window to the sequence rather than discarding it; the buffer is never
  trimmed past the head window's start, so nothing has been lost by then.
  Impact here was one intron in 229,939 (*T. rubripes*) and zero elsewhere,
  because the stall only reaches windows *behind* it and these assemblies put
  few splice sites near scaffold ends; on a more fragmented assembly, or with
  a larger `GC_WINDOW`, it would take out the tail of every affected scaffold.
  The regression fixture above fails on the old code with four wrong values
  and passes on the new.
- **Found by the human run: the three §4.3 strata had three denominators and
  one of them was single-sided.** `by_dinucleotide` is counted per intron,
  `by_intron_length_decile` and `by_local_gc` per site, and §4.3 read as
  though all three were per donor and per acceptor. On a small genome the
  difference is invisible — *S. cerevisiae* has 281 introns over 281 donors —
  but human has 218,446 reference introns over 188,913 donors and 193,359
  acceptors, so the dinucleotide row sums to **15.6% more** than the donor
  total beside it. Per-intron is the only defensible unit for that table: a
  class is the pair, and 106 human reference donors and 644 acceptors sit in
  introns of *two* classes at once — 129 and 2,515 in the GENCODE prediction,
  185 and 260 in the *T. rubripes* reference — so a per-site table would have
  to pick one, and the choice would be arbitrary in exactly the way the
  shortest-intron rule above exists to prevent. The six-to-one acceptor-to-donor
  ratio is itself the mechanism: a shared *donor* usually keeps its class,
  because both its introns end `AG`, while a shared *acceptor* changes class
  whenever one of its introns starts `GC` instead of `GT`. The result
  now states the units (`splice.strata_units`), reports the two-class site
  counts (`reference_sites_multiple_dinuc_classes` and its predicted
  counterpart), and gives the local-GC table for **acceptors as well as
  donors** — before this it existed for donors only, which quietly delivered
  one of the three stratifications on one side of the junction. Fixture: the
  shared-donor gene of the self-test, given a genome in which its two introns
  end `AG` and `AC`.
- **Real predictor output, and what it found.** AUGUSTUS 3.5.0 (bioconda
  `augustus-3.5.0-pl5321h5653ebf_10`) was run ab initio on *S. cerevisiae*
  (`--species=saccharomyces_cerevisiae_S288C`, 47 s wall on 17 cores, 237 MB
  peak) and on *S. pombe* twice, once with its own parameters and once with
  the *S. cerevisiae* ones. Its GFF3 is the shape no RefSeq-versus-RefSeq run
  can reach: `transcript` rather than `mRNA`, no `exon` features, no `region`
  features, no UTR, and the stop codon *outside* the CDS. Three defects
  followed, all fixed and all invisible on a reference:
  1. **`--stop-outside-cds` was a no-op.** The flag was recorded in the result
     and never applied, so the 3 bp were never restored. On *S. cerevisiae*
     that is exon F1 **0.027**, terminal-exon F1 0.000, single-exon F1 0.000
     and stop-codon F1 0.000, against nucleotide F1 0.959 and locus F1 0.919 —
     the healthy-looking numbers are exactly the ones a 3 bp 3' shift does not
     move. Corrected, the same prediction scores exon F1 0.757 and stop-codon
     F1 0.912. The flag now merges `stop_codon` features into the chain, and
     the convention is detected and cross-checked against the flag (§4.5).
  2. **A reused transcript id silently welded chains together.** AUGUSTUS
     restarts its gene numbering at `g1` in every invocation, so concatenating
     per-chromosome output without renaming gives one `g1.t1` per chromosome.
     The loader keyed on the id alone, so 5,154 predicted transcripts became
     663 chimaeras spanning chromosomes, scoring nucleotide F1 0.198 with no
     complaint. Ids that appear on more than one sequence or strand are now
     counted in `predicted_conflicting_transcript_ids` and warned about; the
     naive concatenation reports 466 of them.
  3. **`report.py` dropped a duplicate species without a word.** Two results
     for the same species — a run and its ablation, which is exactly what the
     two *S. pombe* runs are — silently kept the last and printed
     "1 of 20 species scored". It now names both files and exits 2.
  With those fixed, the numbers are the first real ones this benchmark has
  produced, and they say something the design intended:

  | run | nucleotide F1 | exon F1 | donor F1 | transcript F1 | locus F1 |
  |---|---|---|---|---|---|
  | *S. cerevisiae*, own parameters | 0.958 | 0.757 | 0.394 | 0.781 | 0.919 |
  | *S. pombe*, own parameters | 0.955 | 0.774 | 0.854 | 0.702 | 0.923 |
  | *S. pombe*, *S. cerevisiae* parameters | 0.868 | **0.296** | **0.175** | 0.393 | 0.827 |

  Swapping the parameter set between two ascomycete yeasts of nearly the same
  size and GC costs **62% of exon F1 and 80% of donor F1 while nucleotide F1
  falls by 9%**. That is the whole argument for §4.8 reporting a vector rather
  than a headline number: a nucleotide-level score would have called this a
  small degradation. It is also a floor for the charter's clade-independence
  claim — any model claiming it has to beat a parameter swap between two
  yeasts before a mammal-to-fungus claim means anything. The declarations and
  full results are in `benchmark/validation/`.
- **Review feedback from marx (PR #5), acted on.** Pseudogene and
  gene-fragment CDS rows are no longer scored as truth and incomplete CDS ends
  no longer enter the §4.5 denominators (§4, §4.5); the self-test gained a
  fixture built from those RefSeq shapes, including a minus-strand
  `start_range=` that is a *3'* end, and then ran 86 checks, up from 69. Every number in
  this section was recomputed afterwards: the three AUGUSTUS results moved by
  at most 0.002 F1 (the yeasts have 6 and 32 pseudogenes and no partial CDS on
  a scored sequence), and identity runs on human, *C. elegans*,
  *A. thaliana*, *Z. mays*, *Tetrahymena*, *P. falciparum*, *Nematostella*,
  *S. cerevisiae* and *S. pombe* are still 1.0 and MCC 1.0 on every metric
  with 0 fusions and 0 splits. The size of the correction is species-specific
  and not always small: *C. elegans* has **1,958** pseudogene transcripts with
  CDS rows, 6.4% of its 30,548, against none at all in the *Arabidopsis*,
  maize, *Tetrahymena*, *P. falciparum* and *Nematostella* references, which
  annotate pseudogenes without CDS. Partial CDS ends are commonest in
  *Tetrahymena* (253 at the 5' end, 191 at the 3'), human (712 and 1,048 in
  the file) and maize. The committed JSON under `benchmark/validation/` is
  the regenerated output, not the earlier one. `leakage_check.py` now prints
  every training species tied at the deepest MRCA rank (§3.1), the alignment
  leakage rule is split by alignment type (§3.2), and the greedy matching in
  §4.4 is now documented as greedy.
- **A second tool and a second output shape: Helixer 0.3.7 on *T. rubripes*,
  *N. crassa* and *S. cerevisiae*.** Run from the published container
  (`gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8`) on one
  consumer GPU, vertebrate and fungi checkpoints, `--batch-size 8`: fugu
  384 Mb in **91.6 min**, *N. crassa* 41 Mb in 5.8 min, *S. cerevisiae* 12 Mb
  in 1.9 min. This is the input
  shape §7 item 2 asked for and no AUGUSTUS run reaches — `gene`/`mRNA`/
  `exon`/`CDS`/`five_prime_UTR`/`three_prime_UTR` and **no `stop_codon`
  feature at all**, so the convention is decided by the §4.5 genome probe
  rather than by a feature or a flag (`stop_codon_convention_source:
  genome`), and the file carries UTRs, which the CDS-only scoring has to
  ignore. Scoring cost on one laptop core: 34 s / 0.71 GB for fugu, 3.7 s /
  88 MB for *N. crassa*, 1.1 s / 53 MB for *S. cerevisiae*.

  | run | nucleotide F1 | exon F1 | donor F1 | transcript F1 | locus F1 | start F1 | stop F1 |
  |---|---|---|---|---|---|---|---|
  | *N. crassa*, fungi model | 0.962 | 0.772 | 0.841 | 0.689 | 0.906 | 0.795 | 0.842 |
  | *S. cerevisiae*, fungi model | 0.986 | 0.825 | **0.378** | 0.860 | 0.948 | 0.896 | 0.938 |
  | *T. rubripes*, vertebrate model | 0.919 | 0.776 | 0.853 | **0.252** | 0.887 | 0.499 | 0.742 |

  *N. crassa* is the one species on this panel Helixer can be run on without
  declaring pretraining exposure (§3.2 channel 2), so it is the only one of
  the three whose numbers are a measurement rather than an upper bound. Fugu's transcript F1 is not a fugu failure so much as a
  reminder of what the metric measures: the RefSeq reference has 46,771
  transcripts over 22,090 loci and Helixer emits exactly one per locus, so
  three quarters of the reference chains have no candidate to match. It is also
the run that found the third defect below, which was worth 0.014 of it.
Locus F1
  0.887 on the same run is the number to compare against AUGUSTUS's, and the
  1,154 splits (against 32 in *N. crassa*) are where a one-isoform-per-locus
  predictor meets a vertebrate reference. Fugu is also *in the vertebrate
  checkpoint's own training list* at an older assembly, so 0.919 nucleotide F1
  there is an upper bound, as is every column of the *S. cerevisiae* row;
  §3.3's `heldout_seen_in_pretraining` says so in both declarations.

  Three defects, all found by these runs and all invisible on every earlier
  one:
  1. **The §4.3 decile stratification was not reproducible.** Scoring the
     same two fugu files twice gave different per-decile donor and acceptor
     counts while every total was identical. A splice site shared by introns
     of two lengths took whichever the internal set yielded last, which is
     hash-seed dependent; 8,086 fugu donors and 23,387 human ones are shared
     this way. The shortest intron now decides, and the count of sites the
     rule touched is reported (§4.3).
  2. **The genome window plan was made before the stop-codon merge.** Merging
     a `stop_codon` feature that sits across an intron creates a junction the
     plan never asked for, and its dinucleotide was reported as `unknown`
     rather than the GT-AG it is — one predicted intron of the *S. pombe*
     cross-parameter run. The plan is now extended after the merge (§4.5).
  3. **An exact transcript match could lose the within-locus pairing to an
     isoform that merely contains it.** The pairing ordered candidates by
     shared CDS bases with the accession as tie-break, and a containing
     isoform ties with the one the prediction equals, so the tie-break — an
     accession string — decided a scored outcome. 343 of fugu's 5,907 exact
     matches and 34 of *N. crassa*'s 6,923 went to the wrong isoform and were
     counted as a false positive and a false negative each. Exactness now
     outranks overlap (§4.4). This is the one defect so far that a deeply
     annotated reference is required to see: it cannot occur where every
     locus has one isoform, which is why the two yeasts and all three
     AUGUSTUS runs are unaffected, and it moved only the transcript column.
  All six results in `benchmark/validation/` were regenerated with the fixed
  scorer; the three AUGUSTUS numbers above are unchanged to five decimals,
  as are all of the *S. cerevisiae* Helixer ones and every non-transcript
  column of the other two.

  The *S. cerevisiae* row is there to answer item 10, and it does. AUGUSTUS's
  yeast donor F1 of 0.394 is **not** an AUGUSTUS artefact: Helixer, a
  completely different architecture with no shared code or training data,
  scores 0.378 on the same genome, and both do it the same way — by
  predicting far more introns than exist. The reference has 281 scored
  introns; AUGUSTUS predicts 566 and Helixer 734, giving donor precision 0.295
  and 0.262 against sensitivity 0.594 and 0.683. Every other metric on the
  same Helixer run is the best of the six (nucleotide F1 0.986, exon 0.825,
  transcript 0.860), so this is not a bad run. Two tools over-predicting
  introns by 2 to 2.6x in a 95.3% single-exon genome is a property of the
  regime, not of either tool, and it is the argument for keeping an
  intron-poor genome in the panel: an intron-rich reference hides intron
  over-prediction inside a large true-positive count, and this one does not.
  *S. cerevisiae* is in the fungi checkpoint's training list, so the other
  columns of that row are an upper bound; the over-prediction is if anything
  understated by that.
- **Ensembl input.** `Caenorhabditis_elegans.WBcel235.gff3.gz` from the
  Ensembl FTP site — no `region` features at all, `gene:`/`transcript:`
  ID prefixes — parses and scores: 6 nuclear chromosomes kept, `MtDNA`
  dropped by the name fallback, `reference_has_region_features: false` in the
  result, 31,853 transcripts, 19,973 loci, all metrics 1.0 against itself.
- **Naming-mismatch detection.** Scoring the Ensembl *C. elegans* annotation
  against the RefSeq one — same assembly, different sequence names — returns
  F1 0.0 with `predicted_transcripts_not_scored: 31865`,
  `predicted_sequences_absent_from_reference: 7` and a stderr warning, rather
  than a silent zero.
- **Against a synthetically degraded copy of the reference, on all 20 panel
  species** (`benchmark/degrade.py`, seed 20260909, 10% of transcripts deleted
  and the downstream boundary in transcription order of 10% of CDS segments
  moved 3 bp further downstream). This is the control an identity run cannot
  be. An identity run makes every metric 1.0 by construction, so a scorer that
  dropped the same thing from both sides would still pass it; a degraded copy
  has a *known* perturbation, so the direction and the rough size of every
  metric's response can be written down before the run and checked afterwards.
  Full results in `benchmark/validation/degraded/`, one scored JSON and one
  `degrade.py` summary per species.

  The perturbation is deliberately asymmetric, and that asymmetry is the test.
  Which coordinate a "downstream" shift moves depends on the strand -- `end`
  on `+`, `start` on `-` -- so it lands on stop codons and **donor** sites and
  never touches start codons or **acceptor** sites. A scorer whose splice-site
  or codon assignment is not strand-aware cannot reproduce that. The panel
  reproduces it on every species: donor F1 0.811 to 0.870 against acceptor F1
  0.934 to 0.989, with the acceptor and start-codon losses accounted for by
  the deleted transcripts alone.

  | species | nucl | exon | donor | acceptor | transcript | locus | start | stop |
  |---|---|---|---|---|---|---|---|---|
  | *H. sapiens* | 0.985 | 0.820 | 0.811 | 0.989 | 0.378 | 0.986 | 0.975 | 0.844 |
  | *M. musculus* | 0.981 | 0.835 | 0.830 | 0.987 | 0.400 | 0.981 | 0.973 | 0.851 |
  | *G. gallus* | 0.979 | 0.840 | 0.837 | 0.983 | 0.369 | 0.978 | 0.971 | 0.854 |
  | *D. rerio* | 0.976 | 0.844 | 0.841 | 0.980 | 0.377 | 0.976 | 0.970 | 0.857 |
  | *T. rubripes* | 0.972 | 0.860 | 0.859 | 0.975 | 0.365 | 0.971 | 0.966 | 0.861 |
  | *X. tropicalis* | 0.971 | 0.859 | 0.859 | 0.975 | 0.399 | 0.969 | 0.964 | 0.864 |
  | *C. intestinalis* | 0.955 | 0.859 | 0.859 | 0.957 | 0.445 | 0.955 | 0.954 | 0.857 |
  | *N. vectensis* | 0.966 | 0.864 | 0.864 | 0.971 | 0.472 | 0.964 | 0.960 | 0.863 |
  | *D. melanogaster* | 0.976 | 0.862 | 0.859 | 0.978 | 0.621 | 0.972 | 0.969 | 0.869 |
  | *A. mellifera* | 0.976 | 0.856 | 0.855 | 0.977 | 0.454 | 0.973 | 0.966 | 0.856 |
  | *C. elegans* | 0.961 | 0.859 | 0.859 | 0.960 | 0.524 | 0.958 | 0.952 | 0.858 |
  | *A. thaliana* | 0.969 | 0.868 | 0.870 | 0.974 | 0.590 | 0.967 | 0.963 | 0.864 |
  | *O. sativa* | 0.961 | 0.864 | 0.866 | 0.965 | 0.601 | 0.960 | 0.958 | 0.859 |
  | *Z. mays* | 0.964 | 0.863 | 0.863 | 0.969 | 0.585 | 0.962 | 0.960 | 0.860 |
  | *S. cerevisiae* | 0.943 | 0.845 | 0.854 | 0.934 | 0.841 | 0.944 | 0.944 | 0.845 |
  | *S. pombe* | 0.946 | 0.849 | 0.847 | 0.949 | 0.771 | 0.948 | 0.948 | 0.851 |
  | *N. crassa* | 0.951 | 0.854 | 0.856 | 0.955 | 0.725 | 0.951 | 0.951 | 0.852 |
  | *P. falciparum* | 0.945 | 0.843 | 0.842 | 0.941 | 0.730 | 0.945 | 0.945 | 0.845 |
  | *D. discoideum* | 0.946 | 0.850 | 0.849 | 0.946 | 0.749 | 0.946 | 0.946 | 0.851 |
  | *T. thermophila* | 0.947 | 0.852 | 0.851 | 0.947 | 0.649 | 0.948 | 0.948 | 0.856 |

  Four things the panel-wide control establishes that the two-species one did
  not:

  1. **Locus precision is exactly 1.0 and start-codon precision is exactly 1.0
     on all 20 species** -- zero false-positive loci and zero false-positive
     start codons across 775,253 scored reference transcripts. Deleting
     and 3'-shifting can only lose things and move stops, so anything else
     would be the scorer inventing them. The earlier two-species run reported
     human start F1 0.912 *with* false positives, which was correct for the
     ad-hoc script that produced it: it moved the `end` field regardless of
     strand, so on roughly half the annotation it was perturbing start codons
     while the text said stop codons. `degrade.py` exists so that the
     perturbation is a reviewable file rather than a shell one-liner.
  2. **The transcript level is the metric that responds to per-segment error,
     and it does so across the whole panel in the order the panel predicts.**
     Transcript F1 runs from 0.841 on *S. cerevisiae* (1.05 exons per
     transcript) to 0.365 on *T. rubripes* (14.4), while nucleotide F1 over
     the same species moves only from 0.943 to 0.972 -- and *upwards*, because
     an intron-rich genome has more CDS bases to absorb the same 3 bp. A chain
     of `n` segments survives a per-segment shift with probability
     `(1-P)**n`, so the crude floor for transcript sensitivity is
     `0.9 * 0.9**n`. Measured against the `exons_per_tx_mean` column of
     `panel.tsv` the floor is tight where CDS and exon counts nearly coincide
     -- 0.99 to 1.04 of predicted on *S. cerevisiae*, *S. pombe*,
     *D. discoideum*, *P. falciparum* and *N. crassa* -- and is exceeded by up
     to 2.0x in the vertebrates, where that column counts UTR exons the
     perturbation never touches and where the §4.4 isoform rule lets a broken
     chain still match a different reference isoform. Read the ratio as a
     check on ordering, not as a fit.
  3. **Fusion and split counts are not invariant under a perturbation that is
     neither a fusion nor a split.** Six species report 1 or 2 fusions and six
     report 1 to 8 splits, on input that only deletes transcripts and moves
     boundaries 3 bp. Running each perturbation alone at the same seed
     separates the two cleanly and exactly reproduces the combined counts:
     *every* split comes from the deletions (human 8, *A. mellifera* 2,
     *X. tropicalis* 1, *T. rubripes* 0) and *every* fusion from the 3 bp
     shifts (human 1, *A. mellifera* 1, *X. tropicalis* 0, *T. rubripes* 0).
     Both are real and both are §4.4 behaving as specified -- deleting an
     isoform can pull a gene out of the prediction's own CDS-overlap graph and
     split the component it was holding together, and a 3 bp extension can
     bridge a same-strand gene pair that the reference kept 1 to 3 bp apart --
     but the practical consequence is that a fusion or split count in the low
     single digits is inside the noise a structurally *correct* submission can
     produce on a 20,000-gene annotation. Report it; do not rank on it.
  4. **A defect in `degrade.py` itself, found by the ablation and not by the
     combined run.** The first version drew the per-segment shift only for CDS
     lines that survived the deletion pass, which makes the shift stream a
     function of the drop rate: at one seed, a drop-only run and a combined run
     shifted different segments, so the ablation was not a decomposition of the
     combined run at all. *T. rubripes* is where it showed -- the combined run
     reported 1 fusion that neither single perturbation could produce, which is
     impossible, since a combined prediction's CDS blocks are a subset of the
     shift-only one's. The shift is now drawn for every CDS line and applied
     only to the survivors, and the self-test asserts over four seeds that the
     surviving CDS coordinates of a combined run equal the shift-only run's on
     the same segments. The check fails on the pre-fix code.

  Scoring cost, one laptop core per species: 0.68 s and 44 MB on *S. pombe*
  to 52 s and 1.04 GB on *D. rerio* (1.16 M CDS segments); human 43 s and
  1.06 GB. `degrade.py`'s own cost was not measured separately.
- `benchmark/fetch.py --what fasta` is now exercised: it downloads and
  checksum-verifies the 3.8 MB *S. cerevisiae* FASTA. It is still untested at
  3 Gb (§7).
- **A prediction with many isoforms per locus: GENCODE 50 scored against
  RefSeq on human.** Every run above emits at most one transcript per locus,
  so the §4.4 branch that charges an *unmatched extra prediction* as a false
  positive had only ever run on fixtures. GENCODE 50 (Ensembl 116,
  primary-assembly GFF3, released 2026-04-08, sha256
  `272f9972d1bfffa344889c7ae47973bbd79a75e933c55505f347cd25121876fb`) is the
  input that exercises it: 370,910 CDS-bearing transcripts over the same
  assembly the RefSeq reference uses, 2.8x the reference's 131,442. It is not
  a predictor and its numbers are not an accuracy measurement — both files are
  human curation of the same genome — but it is a *real* multi-isoform
  submission, and it is the shape an evidence-based pipeline produces.
  Preparation was sequence renaming only: `chr*` and the GenBank scaffold
  accessions to RefSeq accessions from
  `GCF_000001405.40_GRCh38.p14_assembly_report.txt` (sha256
  `64318ddff470b69b261a667d813210044f60d4ce654253a547db80ff73638d38`), which
  the report maps for all but `KI270721.1` and `KI270734.1`; rows on those two
  (423 of 11,253,000) were dropped before scoring. 64 s and 1.6 GB on one
  laptop core, identical under `PYTHONHASHSEED` 0 and 1. Declaration and full
  result: `benchmark/validation/gencode50-Homo_sapiens.{yaml,json}`.

  | nucleotide F1 | exon F1 | donor F1 | acceptor F1 | transcript F1 | locus F1 | start F1 | stop F1 |
  |---|---|---|---|---|---|---|---|
  | 0.933 | 0.697 | 0.908 | 0.826 | **0.290** | 0.968 | 0.750 | 0.408 |

  Two defects came out of it, both specific to a submission that is not
  RefSeq-shaped:
  1. **The §4 biotype filter read only `gene_biotype`.** RefSeq spells it
     that way; GENCODE spells it `gene_type` and Ensembl spells it `biotype`.
     The filter is applied to the prediction as well as the reference, so
     reading one spelling meant the reference's 390 immunoglobulin and T-cell
     receptor transcripts were dropped as unanswerable while GENCODE's 421
     equivalents were kept and charged as false positives — the submission was
     penalised for answering a question §4 forbids the reference to ask.
     All three spellings are now read. Locus F1 0.959 to 0.968 (false-positive
     loci 1,133 to 723) and nucleotide F1 0.931 to 0.933; every prior run is
     unaffected, because none of the six declares a biotype at all.
  2. **`predicted_transcripts_not_scored` pooled two unrelated causes.** It
     was every predicted transcript minus the scored ones, so biotype-dropped
     rows landed in the field whose documented meaning is "on a sequence the
     reference does not have", and whose warning tells the submitter to check
     their sequence names. On this run it read 434 where 13 transcripts are
     actually off-panel. It now counts sequence exclusion only; the filter's
     own drops are already reported by reason in
     `predicted_transcript_selection`.

  What the run says about the metrics themselves, none of it a defect:
  transcript precision is 0.190 against sensitivity 0.610 — 300,036 predicted
  chains match nothing in RefSeq, which is what §4.4 does to a submission that
  offers 2.8 isoforms for every one the reference has. Locus F1 0.968 on the
  same run is the contrast the two levels exist for: the *loci* agree almost
  perfectly and the *chains* do not. Stop-codon precision 0.272 is mostly the
  91,818 nonsense-mediated-decay transcripts, each ending at a stop RefSeq
  does not annotate as one; terminal-exon F1 0.340 has the same source. The
  one asymmetry worth recording is that GENCODE states 257,153 distinct CDS
  acceptors against 212,078 donors (21% more), where RefSeq states 193,359
  against 188,913 (2.4% more), which is why acceptor F1 0.826 sits well below
  donor F1 0.908. Counted independently from the GFF3 outside the scorer, both
  numbers agree exactly, so this is annotation depth and not a scoring
  artefact; §4.3 reporting donors and acceptors separately is what makes it
  visible.
- **A prediction that is itself partial: AUGUSTUS `--genemodel=partial` on
  *T. thermophila* and *A. mellifera*.** Every submission above is
  `--genemodel=complete` or a curated annotation, so §4.5's predicted-partial
  handling had only ever run on fixtures — and the fixtures were built around
  the *reference* convention, `start_range=`/`end_range=`, which no predictor
  writes. Two runs of AUGUSTUS 3.5.0 with `--genemodel=partial` are the case
  that was missing: *T. thermophila*, one invocation per scaffold over its
  **1,158 scaffolds** under translation table 6, and *A. mellifera*
  (`--species=honeybee1`, 177 sequences, 36.4 min wall on 22 cores, 570 MB
  peak). Both are in AUGUSTUS's own training set, so neither is a
  measurement; they are here for the shapes they produce.

  | run | nucleotide F1 | exon F1 | donor F1 | transcript F1 | locus F1 | start F1 | stop F1 |
  |---|---|---|---|---|---|---|---|
  | *T. thermophila*, `tetrahymena`, code 6 | 0.420 | 0.328 | 0.392 | 0.208 | 0.536 | 0.310 | 0.434 |
  | *A. mellifera*, `honeybee1` | 0.891 | 0.693 | 0.811 | 0.225 | 0.793 | 0.402 | 0.618 |

  Three things came out of them.
  1. **Counting a reused transcript id was not enough; it had to be
     resolved.** The earlier fix counted ids appearing on more than one
     sequence and warned, but still welded their chains together. On a
     1,158-scaffold assembly that is not an edge case: 214 ids collide, and
     the welding takes the *whole file* down to 218 scored chains from 10,355,
     nucleotide F1 **0.020** against 0.420, exon 0.007 against 0.328. A CDS
     chain lies on one sequence and one strand by construction, so the loader
     now keys every id by `(sequence, strand)` and the collision resolves
     instead of merging. The control is exact: the same AUGUSTUS run repeated
     with `--uniqueGeneId=true`, which renames at the source, scores
     identically to the colliding file in every block. The warning now says
     the file is not valid GFF3 and what was done about it, rather than that
     the metrics are meaningless.
  2. **`predicted_partial_5prime`/`_3prime` could not fire.** They read the
     reference's range attributes only, so they were structurally zero on
     every prediction. *T. thermophila* states 250 chains with no
     `start_codon` feature and 69 with no `stop_codon`; reading that omission
     as the statement it is drops 213 5' and 30 3' ends from the codon
     denominators on the scored sequences and removes 208 false-positive
     starts and 22 false-positive stops. `predicted_partial_source` records
     which convention the counts came from, and distinguishes a file that
     uses codon features and declares no partial end (`missing_codon_feature`
     with both counts zero, as the three complete AUGUSTUS runs read) from
     one that cannot state it at all (`none`, as all three Helixer runs read).
  3. **A 3'-partial gene must not be handed a stop codon.** *A. mellifera* is
     the run that exercises it, because `honeybee1` excludes the stop from the
     CDS while `tetrahymena` includes it, so the 3 bp extension actually runs.
     13 chains are declared 3'-partial and skipped; extending them costs 8
     false-positive stop codons, 2 false-positive terminal exons and 27
     false-positive bases, and turns 3 true-positive bases into false
     negatives. The self-test fixture shows the worse form, where the invented
     3 bp land on the reference's real stop and score a **true positive** on a
     gene the predictor never finished.

  **The inference was then checked against a file that states partiality
  both ways.** GENCODE 50 writes `start_codon`/`stop_codon` features *and*
  tags incomplete CDS ends `cds_start_NF`/`cds_end_NF`, so the two are
  independent. Read from the omissions: 13,203 5'-partial and 19,535
  3'-partial. Read from the tags: 13,226 and 19,858. They agree on 13,157 of
  13,203 (99.7%) and 19,371 of 19,535 (99.2%), with 46 and 164 inferred but
  untagged and 69 and 487 tagged but not inferred. That is the accuracy of
  the signal on 370,910 transcripts of a curated human annotation, and it is
  why the GENCODE row's codon columns moved when the inference went in: start
  F1 0.637 to 0.750 (precision 0.591 to 0.821) and stop F1 0.360 to 0.408,
  with every other column identical to five decimals. The tags are not read
  by the scorer — no predictor writes them — but they are the check that the
  omission means what §4.5 says it means.

  All eight results in `benchmark/validation/` were regenerated with the
  scorer that has these three changes. Apart from the two new fields, the
  three AUGUSTUS and three Helixer runs are unchanged in every value; only
  GENCODE's start and stop blocks move.

## 7. Open items

1. **The scorer does not compute §4.6 or §4.7.** BUSCO, OMArk, and the cost
   columns are external and are merged by `report.py --cost`; nothing yet
   produces that TSV. T-human-009 owns the cost half.
2. **Scored submissions cover seven species, three sources and nine runs.**
   AUGUSTUS on the two yeasts, on *T. thermophila* and on *A. mellifera*,
   Helixer on *T. rubripes*, *N. crassa* and *S. cerevisiae*, and GENCODE 50
   on human (§6) between them cover both stop-codon conventions, all three
   detection routes (feature, genome and flag), a prediction with UTRs, a
   reference with several isoforms per locus, a *prediction* with 2.8 isoforms
   per locus, a GENCODE-shaped attribute set, a non-standard genetic code, a
   1,158-scaffold assembly with colliding predicted ids, a prediction that
   declares its own partial ends, and a 384 Mb genome with 105 unplaced
   scaffolds. Ten defects came out of those runs. What is still not covered:
   the §4.5 **blind 3 bp extension** — the fallback for a prediction that has
   no `stop_codon` feature *and* whose genome says `outside` — is still
   fixture-only, because Helixer's convention is `inside` and both AUGUSTUS
   conventions come with the feature; a prediction whose sequence set
   genuinely diverges from the reference's (every submission so far was made
   against the reference assembly, and the GENCODE run was renamed onto it, so
   the divergence is zero); and an *evidence-based pipeline* run end to end,
   which is what would make the multi-isoform case a measurement rather than
   the annotation comparison the GENCODE run is.
   Degraded-copy control runs now cover **20 of 20** species (§6), which
   closes the coverage half of this item: every metric has been checked
   against a known perturbation on every panel genome, and the scorer's
   strand-awareness is pinned by the donor/acceptor and start/stop asymmetry
   that control produces. What remains open is real predictor output, not
   coverage. Tiberius is the obvious next tool: it is the one that would
   exercise the blind extension if its GTF turns out to exclude the stop.
3. **No high-confidence subset** (§2.3). Needed before any accuracy above
   roughly the annotation error rate means anything. Candidate construction:
   loci with MANE Select support in human, community-curated loci elsewhere,
   intersected with peptide or RNA-seq support. Possibly its own task.
4. **Distance rule is taxonomic rank, not distance** (§3.1). Replacing it
   with substitutions per site from a fixed marker set would make the floor
   comparable across kingdoms. Needs a tree; overlaps T-human-008.
5. **Ensembl and predictor input rely on a name heuristic** (§4). Without
   RefSeq `region` features the scorer recognises organelles and alt loci from
   the sequence name, which is verified on Ensembl *C. elegans*, on the
   AUGUSTUS runs and on fixtures, but is not a specification. Note that the
   heuristic never had to act on the AUGUSTUS runs: the *sequence selection is
   derived from the reference*, so a prediction naming the mitochondrion
   `NC_001224.1` is dropped because the reference calls it an organelle, not
   because the name looked like one. Any submission whose *reference* lacks
   `region` features should pass `--seqids` and say so in the declaration.
6. **No RNA-seq accessions chosen** (§3.2 channel 4). Evidence-based tools
   cannot be run on the panel until each species has a declared, fixed
   RNA-seq set. Overlaps T-human-008.
7. **BUSCO and OMArk lineage datasets not pinned** (§4.6), and neither tool
   is stdlib; they are the only external dependencies the benchmark needs.
8. ***Tetrahymena* annotation quality** (§2.3) may be poor enough that even
   "report, never rank" is generous. An alternative code-6 ciliate with a
   better annotation would be preferable if one exists.
9. **Fetch script tested on 20 of 20 annotations, 9 of 20 genome FASTAs,
   closed at 3 Gb.** `--what fasta` works end to end on *S. cerevisiae*,
   *S. pombe*, *P. falciparum*, *C. elegans*, *D. melanogaster*,
   *T. rubripes*, *N. crassa*, *T. thermophila*, *A. mellifera* and — this is
   what closes the item — **human**: `GCF_000001405.40_GRCh38.p14_genomic.fna.gz`,
   972,898,531 bytes for 3.1 Gb of sequence, MD5-verified against the NCBI
   manifest in 19.6 s. What that leaves untested is not scale but the
   remaining eleven species, which are the same code path.
10. **Both tools over-predict introns in *S. cerevisiae*, and the benchmark
   cannot yet say what the right number is.** This item asked for a second
   predictor and now has one (§6): against 281 scored reference introns in a
   95.3% single-exon genome (Table 2), AUGUSTUS predicts 566 (donor F1 0.394)
   and Helixer 734 (donor F1 0.378). Agreement between two unrelated
   architectures rules out the tool-specific explanation and leaves two: both
   are genuinely bad at deciding *not* to splice when introns are rare, or
   some of those 566 and 734 are real introns the reference does not
   annotate — *S. cerevisiae* introns are short and RefSeq's yeast annotation
   is CDS-first. Deciding between them needs the high-confidence subset of
   item 3 or RNA-seq junction support (item 6), not another predictor.
11. **The local-GC stratification (§4.3) degenerates on AT-rich genomes.**
   The five bins are fixed absolute GC bands, which is what makes the column
   comparable across species, but the panel deliberately spans 19.5% to 48.5%
   GC: 98.9% of *P. falciparum* donors land in the `<30%` bin and 83% of
   *S. pombe* donors in `30-40%`. Human is the counter-case that makes the
   fixed bands worth keeping: its 188,913 donors spread 7,913 / 50,492 /
   43,960 / 45,405 / 41,143 across the five, so the column carries real
   signal there. Read it across species, not within one. Per-species GC
   quantile bins would be the alternative and would not be comparable; this
   is a documented limitation, not a defect.
12. **Which non-`protein_coding` biotypes to exclude is a judgement call.**
   §4 drops immunoglobulin and T-cell receptor segments along with
   pseudogenes, because they are gene fragments assembled somatically and
   have no start or stop codon of their own. That is 390 of the 132,030 human
   transcripts on a scored sequence and nothing at all in the other 19 species, so the decision
   costs almost nothing today; it is recorded because the alternative —
   scoring them and accepting that no predictor can produce them — is one
   `--score-all-transcripts` run away, and because a proteome-completeness
   number (§4.6) may turn out to depend on it.
13. **Partial-CDS handling, closed.** §4.5 excludes an incomplete reference
   end from the codon denominators and excludes any predicted codon inside
   that span; what had never been run was the case that matters, a fragmented
   assembly where the *prediction* is also truncated at contig ends. The
   AUGUSTUS `--genemodel=partial` run on *T. thermophila*'s 1,158 scaffolds
   (§6) is that case, and it showed that `predicted_partial_5prime` and
   `predicted_partial_3prime` could not fire at all, because they read a
   reference attribute no predictor writes. They now read the predictor's own
   statement — an omitted `start_codon`/`stop_codon` feature — checked
   against GENCODE's independent `cds_start_NF`/`cds_end_NF` tags at 99.7%
   and 99.2% agreement over 370,910 transcripts (§6). The extension of §4.5
   that keeps a truncated gene from being handed a stop codon it never
   predicted came out of the *A. mellifera* run, which is where the
   stop-outside convention makes the extension run at all. What is still not
   covered: a predictor that emits partial genes *and* no codon features at
   all, which states nothing and for which nothing can be inferred.
