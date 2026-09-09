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
| *Danio rerio* | GRCz12ab annotation is from 2026-07-20 and will move; the MD5 in `panel.tsv` is the pin. |
| all | isoform choice matters: a locus with 40 annotated transcripts is easy to hit at gene level and hard at transcript level (§4.4). |

**A high-confidence subset is required and does not exist yet.** Above some
accuracy level the benchmark stops measuring the model and starts measuring
its ability to reproduce known annotation errors. §7 carries this as the
largest open item.

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
ok	Gallus_gallus	nearest=Mus_musculus	mrca=Sarcopterygii (SUPERCLASS)
ok	Ciona_intestinalis	nearest=Mus_musculus	mrca=Chordata (PHYLUM)
ok	Nematostella_vectensis	nearest=Mus_musculus	mrca=Metazoa (KINGDOM)
ok	Schizosaccharomyces_pombe	nearest=Saccharomyces_cerevisiae	mrca=Ascomycota (PHYLUM)
ok	Plasmodium_falciparum	nearest=Mus_musculus	mrca=Eukaryota (DOMAIN)
ok	Tetrahymena_thermophila	nearest=Mus_musculus	mrca=Eukaryota (DOMAIN)
```

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
3. **Alignment leakage.** This is the one specific to comparative methods and
   the one the charter's design directly invites. A whole-genome multiple
   alignment (multiz, Cactus —
   [10.1038/s41586-020-2871-y](https://doi.org/10.1038/s41586-020-2871-y))
   built to include a held-out species leaks that species' *sequence* into
   every training window that overlaps it, and if the alignment was filtered,
   scored, or projected using annotation, it leaks labels too. Rules:
   - Training alignments must be rebuilt with held-out species removed from
     the taxon set, not merely masked at read time. A row dropped after the
     alignment was computed still shaped the columns.
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
alignment: <how built, which taxa> # or "none"
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

- **canonical vs non-canonical** dinucleotides (GT-AG, GC-AG, AT-AC, other);
- **intron length decile**, computed by the scorer from the reference
  annotation of the species being scored rather than read from `panel.tsv`,
  so that the long-intron tail is visible instead of averaged away. The cuts
  are emitted with the result (`splice.intron_length_decile_cuts`);
  `panel.tsv` carries only p10/median/p90/p99, which is enough to justify the
  panel and not enough to bin against.
- **local GC** in a 200 bp window, in five bins.

The intron-length stratification is the single most informative panel in the
whole benchmark for the charter's question, because it is where clade-specific
models are expected to differ from a species-independent one.

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

Isoform rule: a predictor emitting one transcript per locus is scored against
the **single best-matching annotated isoform** per locus, and the remaining
annotated isoforms are not counted as false negatives. A predictor emitting
several isoforms is scored with an optimal one-to-one matching within the
locus, and unmatched predictions are false positives. Without this rule,
species with deep isoform annotation (human, zebrafish) are systematically
penalized against species with one transcript per gene (*P. falciparum*).

### 4.5 Start and stop codons

Reported separately from exon boundaries, because they are the part of the
problem where the genetic code is a prior rather than a signal: exact-position
sensitivity and precision for start codons and for stop codons. *Tetrahymena*
is the diagnostic — under genetic code 6 a model that has hard-coded TAA/TAG
as terminators will show near-zero stop-codon precision there and normal
numbers everywhere else, which is a signature no aggregate score would show.

### 4.6 Proteome completeness

Translate the predicted CDS and run BUSCO in protein mode
([10.1093/molbev/msab199](https://doi.org/10.1093/molbev/msab199)) against
the deepest lineage dataset that applies to the species, reporting complete,
duplicated, fragmented and missing. Run it on the *reference* proteome too
and report both: BUSCO scores a proteome against a expectation, not against
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

Two decisions the implementation forced, both above: the intron-length
deciles are computed from the reference at score time (§4.3), and a CDS gap
under 20 bp is not a splice junction (§4.3).

Verification so far:

- `python3 benchmark/score.py --self-test` scores built-in fixture pairs and
  checks 41 expected counts: a prediction with one exact transcript, one
  shifted minus-strand boundary, one overlapping-but-unaligned locus and one
  spurious locus; a fusion-and-split pair; a self-comparison that must score
  exactly 1.0 on every metric; two sequence-selection fixtures, one in RefSeq
  shape and one in Ensembl shape, whose every line is a sequence shape taken
  from a real panel reference (chromosome, alt locus, unlocalized scaffold,
  unplaced scaffold, bare scaffold, mitochondrion, chloroplast, sub-10 kb
  scaffold); the streaming FASTA window reader against a plain read; and a
  two-record FASTA in which a window overruns the end of a record, on the
  first record and on the last, checking that it is clipped and that the
  windows queued behind it are still served.
- **Identity runs on all twenty panel references.** Every species in
  `panel.tsv` scored against its own reference gives F1 = 1.0 and MCC = 1.0 on
  every metric, with **0 fusions and 0 splits** everywhere. Cost on one laptop
  core ranges from 1.0 s / 43 MB (*P. falciparum*, 23 Mb) to 52 s / 335 MB
  (*Z. mays*, 2.18 Gb); human is 38 s / 0.92 GB over 3.10 Gb. The filter keeps
  3 sequences (*S. pombe*) to 685 (*Z. mays*); human is 102 scored sequences
  and 3,101,538,863 bp (24 chromosomes plus 78 scaffolds; 511 alt loci and
  patches, 91 sub-10 kb scaffolds and the mitochondrion dropped), 132,030
  transcripts, 20,520 loci. Every one of the twenty took the RefSeq `region`
  path (`reference_has_region_features: true`).
- **`--genome` across six species and at vertebrate scale.** The §4.3
  dinucleotide and local-GC strata now run on *T. rubripes* (384 Mb, 230,048
  introns, 32 s and 0.69 GB), *D. melanogaster*, *C. elegans*, *S. pombe*,
  *P. falciparum* and *S. cerevisiae*. Non-GT-AG donor fractions recovered
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
  Impact here was one intron in 230,048 (*T. rubripes*) and zero elsewhere,
  because the stall only reaches windows *behind* it and these assemblies put
  few splice sites near scaffold ends; on a more fragmented assembly, or with
  a larger `GC_WINDOW`, it would take out the tail of every affected scaffold.
  The regression fixture above fails on the old code with four wrong values
  and passes on the new.
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
- Against a synthetically degraded copy of a reference (10% of transcripts
  deleted, 10% of CDS 3' boundaries shifted by 3 bp): on *S. cerevisiae*,
  nucleotide F1 0.947, exon F1 0.856, donor F1 0.892, transcript F1 0.852,
  locus F1 0.949, recovering 270 GT-AG, 8 GC-AG and 18 other donor
  dinucleotides from the genome FASTA; on human, nucleotide F1 0.984, exon F1
  0.821, donor F1 0.891, acceptor F1 0.895, transcript F1 0.366, locus F1
  0.986 with 0 fusions and 7 splits, start F1 0.912, stop F1 0.907. Human
  transcript F1 falls that far because it averages 13.9 CDS exons per
  transcript, so a 10% per-exon boundary shift breaks about three quarters of
  the chains; that sensitivity is the point of reporting the transcript level
  separately from the nucleotide level.
- `benchmark/fetch.py --what fasta` is now exercised: it downloads and
  checksum-verifies the 3.8 MB *S. cerevisiae* FASTA. It is still untested at
  3 Gb (§7).

## 7. Open items

1. **The scorer does not compute §4.6 or §4.7.** BUSCO, OMArk, and the cost
   columns are external and are merged by `report.py --cost`; nothing yet
   produces that TSV. T-human-009 owns the cost half.
2. **Degraded-copy runs cover two of twenty species.** Identity runs now
   cover all twenty and `--genome` covers six including a vertebrate (§6), so
   what is left untested is behaviour under *wrong* input: only
   *S. cerevisiae* and *H. sapiens* have been scored against a degraded copy,
   and no real predictor output has ever been scored. An identity run
   exercises every code path but pins only the fixed points; the degraded runs
   are what show the metrics move in the right direction and by how much.
   The 14 remaining species need one each, and the first real GFF3 out of
   AUGUSTUS or Helixer will exercise the §4 conventions (stop codon in or out
   of the CDS, `Parent` shapes, missing `region` features) that no
   RefSeq-versus-RefSeq run can reach.
3. **No high-confidence subset** (§2.3). Needed before any accuracy above
   roughly the annotation error rate means anything. Candidate construction:
   loci with MANE Select support in human, community-curated loci elsewhere,
   intersected with peptide or RNA-seq support. Possibly its own task.
4. **Distance rule is taxonomic rank, not distance** (§3.1). Replacing it
   with substitutions per site from a fixed marker set would make the floor
   comparable across kingdoms. Needs a tree; overlaps T-human-008.
5. **Ensembl and predictor input rely on a name heuristic** (§4). Without
   RefSeq `region` features the scorer recognises organelles and alt loci from
   the sequence name, which is verified on Ensembl *C. elegans* and on
   fixtures but is not a specification. Any submission whose reference lacks
   `region` features should pass `--seqids` and say so in the declaration.
6. **No RNA-seq accessions chosen** (§3.2 channel 4). Evidence-based tools
   cannot be run on the panel until each species has a declared, fixed
   RNA-seq set. Overlaps T-human-008.
7. **BUSCO and OMArk lineage datasets not pinned** (§4.6), and neither tool
   is stdlib; they are the only external dependencies the benchmark needs.
8. ***Tetrahymena* annotation quality** (§2.3) may be poor enough that even
   "report, never rank" is generous. An alternative code-6 ciliate with a
   better annotation would be preferable if one exists.
9. **Fetch script tested on 20 of 20 annotations, 6 of 20 genome FASTAs.**
   `--what fasta` works end to end on *S. cerevisiae*, *S. pombe*,
   *P. falciparum*, *C. elegans*, *D. melanogaster* and *T. rubripes*
   (391 Mb, the largest so far, checksum verified). It has still not been
   exercised at 3 Gb scale.
10. **The local-GC stratification (§4.3) degenerates on AT-rich genomes.**
   The five bins are fixed absolute GC bands, which is what makes the column
   comparable across species, but the panel deliberately spans 19.5% to 48.5%
   GC: 98.9% of *P. falciparum* donors land in the `<30%` bin and 83% of
   *S. pombe* donors in `30-40%`. Read that column across species, not within
   one. Per-species GC quantile bins would be the alternative and would not be
   comparable; this is a documented limitation, not a defect.
