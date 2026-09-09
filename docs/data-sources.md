# Data sources: alignments, conservation, expression and annotations

Task T-human-008. Author: marx. Status: draft 1, 2026-09-09.

Every number in this document was either read from a directory listing or an
API response on 2026-09-09 (UTC) from the marx runner, in which case it says
so, or comes from a cited paper. Sizes are the byte counts the servers
reported that day. `scripts/data/fetch_window.py` is the accompanying tool;
its manifests (`*.manifest.json`) record every URL it touched and the SHA-256
of every file it wrote, and `scripts/data/README.md` lists the demonstration
runs.

Companion documents: `docs/benchmark.md` (T-human-007, lenin) defines the
species panel in `benchmark/panel.tsv` and the leakage rules that section 7
here applies. The panel's split labels are used below without redefining
them: `train`, `heldout`, `heldout_paired`.

## 1. Summary

| Source | What it has that we need | Coverage of the panel | Access | Licence / limits |
|---|---|---|---|---|
| UCSC Genome Browser | multiz and Cactus multiple alignments referenced on human, mouse, chicken, fly, worm, yeast, with trees; phyloP and phastCons; RefSeq annotation for 238 native assemblies and 52,623 GenArk assembly hubs; expression, regulation and variation tracks for human and mouse | alignments for 6 of 20 panel species (some on older assemblies); RefSeq annotation for 19 of 20 | JSON API (`api.genome.ucsc.edu`), rsync/https (`hgdownload.soe.ucsc.edu`), public MySQL | no licence needed for data; API guidance is about one request per second, `maxItemsOutput` at most 1,000,000 |
| Ensembl Compara (release 116) | EPO, EPO-extended, PECAN and Cactus multiple alignments for vertebrates with per-block trees and inferred ancestral sequences; GERP scores and constrained elements; Ensembl annotation | human, mouse, chicken (same assembly as the panel), zebrafish, fugu, frog, ciona (some on older assemblies); nothing for invertebrates, plants, fungi, protists beyond rice | REST (`rest.ensembl.org`), FTP MAF/EMF dumps | open data; REST limit 55,000 requests per hour per client (response header) |
| NCBI RefSeq / GenBank | the panel's reference annotations and genomes with MD5 sums; taxonomy | all 20 | https/FTP, Datasets, E-utilities | public domain; E-utilities 3 requests per second without an API key |
| Zoonomia / Cactus consortia | 241-mammal Cactus HAL (2020) and 447-way (2023) alignments; phyloP from them | human (and any mammal in the HAL, via `halLiftover`) | https (UCSC CGL), UCSC tracks | open; HAL tools needed for anything but the hg38-referenced bigMaf |
| Community databases in the panel | FlyBase, WormBase, TAIR/Araport, SGD, PomBase annotations; VEuPathDB RNA-seq for Plasmodium | as listed in `panel.tsv` | through NCBI RefSeq mirrors and GenArk `contrib` tracks | per-database, all open |

The headline: **public multiple alignments exist only for vertebrates,
flies, nematodes and yeasts, and only on a few reference assemblies.** For
twelve of the twenty panel species (frog, honey bee, sea anemone, ciona,
thale cress, rice, maize, fission yeast, Neurospora, slime mould,
Plasmodium, Tetrahymena) there is no public multiple alignment at all, and
for two more (zebrafish, fugu) the alignment is on an older assembly than
the panel's. Any comparative model evaluated on those species must either
run without informants or on alignments we build ourselves. Section 8
costs this.

## 2. UCSC Genome Browser

### 2.1 Multiple alignments

Listing of `hgdownload.soe.ucsc.edu/goldenPath/<db>/` and `/gbdb/<db>/` on
2026-09-09. "Compressed" is the `maf/*.maf.gz` directory total; "raw" is the
uncompressed per-chromosome MAF that backs the browser track and is served
under `/gbdb/` with HTTP `Accept-Ranges: bytes` (verified with HEAD).

| Assembly | Track | Aligner | Year | Species | Format on hgdownload | Compressed | Raw | Panel relevance |
|---|---|---|---|---|---|---|---|---|
| hg38 | multiz100way | multiz | 2015 | 100 | maf.gz per chromosome + `.nh` | 71.4 GB (357 files) | 790 GB | human `heldout_paired` |
| hg38 | multiz30way | multiz | 2017 | 30 | maf.gz + `.nh` | 18.0 GB | 156 GB | human |
| hg38 | multiz470way | multiz | 2022 | 470 | bigMaf (API only; no maf.gz directory) + `.nh` | n/a | n/a | human; largest vertebrate set |
| hg38 | cactus241way (`cactus241wayBM`) | Cactus (Zoonomia) | 2020 | 241 mammals | bigMaf (API) + `.nh`; phyloP bigWig 9.6 GB | n/a | n/a | human |
| hg38 | cactus447way | Cactus (Zoonomia + primates) | 2023 | 447 | bigMaf (API) + `.nh.txt`; phyloP bigWig 10.0 GB | n/a | n/a | human |
| mm39 | multiz35way | multiz | 2021 | 35 | maf.gz + `.nh` | 16.3 GB | 141 GB | mouse `train` |
| galGal6 | multiz77way | multiz | 2019 | 77 | maf.gz + `.nh` | 34.9 GB | 334 GB | chicken `heldout`, but GRCg6a, not the panel's GRCg7b |
| dm6 | multiz124way | multiz | 2019 | 124 insects | maf.gz + `.nh` | 5.4 GB | 63.6 GB | fly `train` |
| dm6 | multiz27way | multiz | 2014 | 27 | maf.gz | 2.6 GB | (not under `/gbdb`) | fly |
| ce11 | multiz135way | multiz | 2019 | 135 nematodes | raw MAF under `/gbdb/ce11/multiz135way/chr*.maf` only | none | 48.5 GB | worm `train` |
| ce11 | multiz26way | multiz | 2016 | 26 | (directory empty) | none | n/a | worm |
| sacCer3 | multiz7way | multiz | 2011 | 7 | maf.gz | 0.1 GB | 0.1 GB | yeast `train` |
| danRer11, xenTro10 | none | | | | no `*way` directories | | | zebrafish `train`, frog `train` |
| apiMel2 | none | | | | 2005 assembly; no alignment | | | honey bee `heldout` |
| (no UCSC db) | | | | | Arabidopsis, rice, maize, S. pombe, Neurospora, Dictyostelium, Plasmodium, Tetrahymena, Nematostella, Ciona, fugu (fr3 exists but has no `*way`) | | | |

Notes.

- The tree for each track sits next to it as Newick with branch lengths:
  `hg38.470way.nh` (19 KB, 470 leaves), `hg38.cactus241way.nh` (11 KB),
  `hg38.447way.nh.txt` (23 KB), `dm6.124way.sequenceNames.nh` (7 KB),
  `ce11.135way.nh` (7 KB). Leaf labels are the MAF source names: UCSC
  database names for multiz tracks (`panTro6`), scientific names for the
  Cactus tracks (`Pan_troglodytes`). `fetch_window.py` checks that every
  source in the fetched blocks is a leaf of the tree; on the demonstration
  windows the match was 461 of 461 (hg38 470-way) and 89 of 89 (dm6).
- The UCSC JSON API behaves differently for the two MAF track types. For
  bigMaf tracks (470-way, both Cactus tracks, `hprc90way`) `getData/track`
  returns the MAF blocks themselves. For wigMaf tracks (every other multiz
  track) it returns only an index: `extFile` id and byte `offset` of each
  overlapping block. Because the raw MAF is served under `/gbdb/` with byte
  ranges, a window can still be fetched exactly with a handful of `Range`
  requests; the fly Adh demonstration pulled 352 blocks (4.6 MB) in 10
  requests. This is how `fetch_window.py` avoids downloading a 12 GB
  chromosome to read 4 kb of it.
- Volume per window. The human beta-globin window (4,931 bp) is 15.4 MB of
  MAF in the 470-way and about 58 MB of JSON in the 241-way Cactus track,
  which has many short blocks (2,064 blocks for 3.9 kb). That is roughly 3
  KB of alignment text per reference base at 470 species. Training data
  cannot be stored as MAF text at that depth; it has to be streamed and
  reduced (per-species one-hot or per-clade summaries) at cut time.

### 2.2 Conservation

phyloP (Pollard et al. 2010, doi:10.1101/gr.097857.109) and phastCons
(Siepel et al. 2005, doi:10.1101/gr.3715005) bigWigs on hgdownload; the API
returns per-base values for any window (`getData/track?track=phyloP470way`).

| File | Bytes |
|---|---|
| hg38 phyloP100way | 9,870,053,206 |
| hg38 phyloP470way | 11,511,784,023 |
| hg38 phyloP447way | 10,022,801,463 |
| hg38 cactus241way phyloP | 9,644,660,543 |
| hg38 phastCons100way | 5,886,377,734 |
| hg38 phastCons470way | 5,060,822,729 |
| dm6 phyloP124way | 470,242,283 |
| ce11 phyloP135way | 356,259,120 |

phastCons elements (`phastConsElements*way`) are bigBed. GenArk hubs (2.4)
carry no conservation tracks. Both scores are unsupervised functions of the
alignment and the tree, so they carry alignment leakage (benchmark section
3.2, channel 3) but not label leakage.

### 2.3 Annotation

`ncbiRefSeq` is a genePred table on the native assemblies and a
`bigGenePred` on GenArk hubs; both come back as JSON from the API and
`fetch_window.py` normalises them to one transcript format (0-based
half-open exon and CDS intervals). RefSeq, GenBank and Ensembl sequence
names map through the `chromAlias` table (`NC_000011.10` and `CM000673.2`
and `11` all resolve to `chr11`), which is the translation lenin asked for
in message 20260909T052140Z-lenin-0007 for scoring Ensembl-named
predictions against RefSeq references.

### 2.4 GenArk assembly hubs

`api.genome.ucsc.edu/list/genarkGenomes` reported 52,623 assemblies on
2026-09-09. A hub exists for 18 of the 20 panel accessions (probe of
`hubs/GCF/…/hub.txt`); the exceptions are human GRCh38.p14
`GCF_000001405.40`, which is covered by the native `hg38`, and zebrafish
GRCz12ab `GCF_052040795.1` (annotated 2026-07-20; not yet built). Each hub
has RefSeq genes as bigGenePred, `augustus`, repeat and GC tracks, and for
some organisms `contrib` tracks mirrored from partner databases: the
Plasmodium falciparum hub carries VEuPathDB RNA-seq coverage
(`Combined_RNA-Seq`, bigWig), `Annotated_Introns` and `Unannotated_Introns`
(introns confirmed or suggested by RNA-seq reads, bigBed) and long-read
transcript models. GenArk hubs have no alignments and no conservation.

### 2.5 Expression, regulation and variation

Human (hg38) top-level tracks seen in `list/tracks` on 2026-09-09 that a
UTR or expression-aware model could use: `gtexGeneV8`, `gtexTranscExpr`,
`gtexCov` (GTEx expression and coverage), `recount3` (uniformly processed
RNA-seq), `encode4LongRnaTranscripts`, `wgEncodeReg4RnaSeq` and the ENCODE
regulation set (DNase, CTCF, H3K27ac, H3K4me3, TF ChIP), `dbSnp155Composite`,
`gnomadVariantsV4.1`, `gnomadConstraint`, `clinvar`. Mouse has a subset.
Nothing comparable exists on UCSC for the other panel species; RNA-seq
evidence for them comes from the RefSeq annotation runs (not distributed as
tracks) or partner databases (2.4).

### 2.6 Access and policy

- JSON API: `https://api.genome.ucsc.edu` (`list/ucscGenomes`,
  `list/tracks`, `list/chromosomes`, `getData/track`, `getData/sequence`);
  `maxItemsOutput` maximum 1,000,000 (the server rejects larger values
  with HTTP 400); the help page recommends at most one request per second
  and says there is no strict limit.
- Bulk: `rsync -aP rsync://hgdownload.soe.ucsc.edu/goldenPath/<db>/<track>/ .`
  or https; `md5sum.txt` in every directory.
- Public MySQL (`genome-mysql.soe.ucsc.edu`): the help page says bot access
  and excessive program-driven use are not permitted; use the API or
  rsync instead.
- Licence: `genome.ucsc.edu/license` states that no licence is needed for
  the data files and database tables, for academic or commercial use.
  Individual tracks may carry their own terms (check the track description).

## 3. Ensembl Compara (release 116, 2026-08)

### 3.1 Multiple alignments

From `ftp.ensembl.org/pub/release-116/maf/ensembl-compara/multiple_alignments/`
(sizes summed from the listing) and `rest.ensembl.org/info/compara/species_sets/<method>`.

| Set | Method | Species | MAF dump | Panel species included (assembly) |
|---|---|---|---|---|
| 44 eutherian mammals | EPO | 44 | 59.2 GB (1,202 files) | human GRCh38, mouse GRCm39 |
| 92 eutherian mammals | EPO-extended | 92 | 73.3 GB | human, mouse |
| 22 murinae | EPO | 22 | 33.0 GB | mouse |
| 10 / 24 primates | EPO / EPO-extended | 10 / 24 | 15.5 / 18.9 GB | human |
| 17 / 27 sauropsids | EPO / EPO-extended | 17 / 27 | 7.2 / 7.4 GB | chicken bGalGal1.mat.broiler.GRCg7b (the panel's assembly) |
| 32 / 65 fish | EPO / EPO-extended | 32 / 65 | 9.0 / 8.7 GB | zebrafish GRCz11 (panel: GRCz12ab), fugu fTakRub1.2 (panel: 1.3) |
| 60 amniotes | PECAN | 60 | 41.1 GB | human, mouse, chicken |
| 31 primates, 26 rodents | Cactus (`CACTUS_DB`) | 31 / 26 | (REST only) | human; mouse |
| 10 fowl | Cactus HAL | 10 | (REST only) | chicken |

Xenopus tropicalis UCB_Xtro_10.0 and Ciona intestinalis KH are in Ensembl but
in no multiple alignment (pairwise LASTZ only). Ensembl Genomes: Plants has
only an 8-way and an 11-way rice EPO alignment (0.8 and 0.7 GB); Metazoa
has 16 pairwise LASTZ-net tarballs and no multiple alignment; Fungi and
Protists have no `maf` directory at all (HTTP 404 on 2026-09-09).

### 3.2 What the REST alignment endpoint gives

`GET /alignment/region/{species}/{region}?species_set_group=…;method=…`
returns, per block, a Newick tree with branch lengths whose leaves are
`species_region_start_end[strand]`, and the aligned rows. EPO blocks include
the inferred ancestral sequences as extra rows (labelled like
`Ggal-Mgal[2]`), which are the internal nodes of that tree: a model that
wants ancestral states as input gets them for free here, and nowhere else.
Any species in the set can be the query, so mouse-referenced blocks of the
44-mammal alignment are one call away (the mouse `Hbb-bs` run in
`scripts/data/README.md`); the FTP `.maf` dumps are referenced on one
species per set, while the `.emf` files carry the full graph.

Rate limit, from the response headers on 2026-09-09:
`x-ratelimit-limit: 55000`, `x-ratelimit-period: 3600`. Annotation through
`/overlap/region` (gene, transcript, exon, cds), sequence through
`/sequence/region`, GERP constrained elements through
`/overlap/region?feature=constrained;species_set_group=…` (11 elements in
the beta-globin window).

Method references: EPO, Paten et al. 2008 (doi:10.1101/gr.076554.108);
PECAN, Paten et al. 2008 (doi:10.1101/gr.076521.108); Cactus, Armstrong et
al. 2020 (doi:10.1038/s41586-020-2871-y); GERP, Davydov et al. 2010
(doi:10.1371/journal.pcbi.1001025).

Licence: Ensembl distributes its data without restriction; the exact text
is on `www.ensembl.org/info/about/legal/disclaimer.html` (not re-read this
tick; the page did not return the licence paragraph to the fetcher).

## 4. NCBI RefSeq and GenBank

- The panel's genomes and annotations are the `GCF_*` assembly directories
  under `https://ftp.ncbi.nlm.nih.gov/genomes/all/`, resolved and
  checksummed by `benchmark/fetch.py` on lenin's T-007 branch; `panel.tsv`
  records the annotation release, date and the GFF MD5 for all 20 species.
  There is no need for a second manifest here.
- The RefSeq eukaryotic annotation pipeline aligns RNA-seq and proteins for
  its Gnomon models; those alignments are not distributed as files, only
  summarised in each assembly's annotation report. Where they surface as
  tracks is GenArk (2.4) and NCBI's own assembly hubs.
- E-utilities: at most 3 requests per second without an API key, 10 with
  one (NCBI, "A General Introduction to the E-utilities", NBK25497). The
  FTP site states no numeric limit; `benchmark/fetch.py` already paces
  itself.
- Licence: NCBI data are in the public domain (NLM copyright policy, "Information
  that is created by or for the US government on this site is within the
  public domain").

## 5. Zoonomia and the large Cactus alignments

- `241-mammalian-2020v2.hal`, the Zoonomia Cactus alignment (Zoonomia
  Consortium 2020, doi:10.1038/s41586-020-2876-6; conservation in Christmas
  et al. 2023, doi:10.1126/science.abn3943): 864,878,803,587 bytes
  (`Content-Length` from `https://cgl.gi.ucsc.edu/data/cactus/` on
  2026-09-09). HAL is reference-free; extracting a window referenced on any
  of the 241 genomes needs `hal2maf` from the HAL tools (Hickey et al. 2013,
  doi:10.1093/bioinformatics/btt128), a C++ dependency that would need a
  `proposal` under the charter's dependency rule. The hg38-referenced
  projection is on UCSC as `cactus241wayBM` (2.1) and needs nothing.
- The 447-way alignment (Zoonomia plus the 2023 primate genomes) is on UCSC
  as `cactus447way` with `phyloP447way`; no HAL download was located this
  tick.
- Zoonomia's value for us is depth within mammals; it adds nothing for the
  non-mammal half of the panel.

## 6. Non-coding alignment: how well introns, UTRs and intergenic sequence align

Measured on the three demonstration windows with the coverage table that
`fetch_window.py` writes: for each informant, the fraction of reference
bases of each annotation class (CDS, UTR exon, intron, intergenic, union
over overlapping RefSeq transcripts, CDS taking precedence) at which the
informant has an aligned non-gap base. One locus per clade is an anecdote,
not a statistic; the point is the shape, which is the same in all three.

Human beta-globin (`HBB`, chr11:5,224,964-5,229,895, 4,931 bp; 444 CDS, 184
UTR, 980 intron, 3,323 intergenic bases), hg38 multiz470way, 460 informants:

| Informant | CDS | UTR | intron | intergenic |
|---|---|---|---|---|
| chimpanzee panTro6 | 1.000 | 1.000 | 0.971 | 0.976 |
| dog canFam4 | 1.000 | 0.935 | 0.535 | 0.837 |
| cow bosTau9 | 0.986 | 0.935 | 0.554 | 0.856 |
| mouse mm39 | 1.000 | 0.935 | 0.517 | 0.325 |
| opossum monDom5 | 1.000 | 0.848 | 0.092 | 0.023 |
| mean over 460 informants | 0.948 | 0.876 | 0.527 | 0.655 |
| median | 1.000 | 0.935 | 0.540 | 0.815 |
| share of informants above 0.5 | 0.95 | 0.94 | 0.73 | 0.70 |

Fruit fly `Adh` (chr2L:14,615,052-14,619,402, 4,350 bp; `Adh` lies inside an
intron of `osp` on the other strand, so the window has no intergenic
bases: 1,590 CDS, 664 UTR, 2,096 intron), dm6 multiz124way, 88 informants:

| Informant | CDS | UTR | intron |
|---|---|---|---|
| D. simulans droSim2 | 1.000 | 0.989 | 0.784 |
| D. yakuba droYak3 | 1.000 | 0.985 | 0.752 |
| D. pseudoobscura droPse3 | 0.996 | 0.696 | 0.537 |
| D. virilis droVir3 | 0.976 | 0.587 | 0.279 |
| house fly musDom2 | 0.879 | 0.039 | 0.111 |
| honey bee apiMel4 | 0.086 | 0.000 | 0.014 |
| mean over 88 informants | 0.677 | 0.307 | 0.213 |
| median | 0.739 | 0.113 | 0.084 |

Chicken `GAPDH` (Ensembl sauropsids EPO, 17 species, of which 5 align here;
6,638 bp; 1,275 CDS, 387 UTR, 4,334 intron, 642 intergenic): the two blocks
covering the window carry turkey, great tit, canary and zebra finch, plus
ancestral rows; the coverage table is in the demonstration manifest.

What the numbers say, and what they do not:

1. **CDS alignability saturates; intron alignability is a within-clade
   feature.** In the HBB window every mammal aligns the coding exons and 73%
   of the 460 informants align more than half the introns, but the
   marsupial aligns 9% of intron bases and 2% of intergenic. In flies the
   fall-off is steeper: within the *melanogaster* group introns are three
   quarters aligned, by *D. virilis* under a third, and outside Drosophila
   essentially nothing but CDS remains. So a comparative signal on introns
   and UTRs is available only from close informants, and the useful
   informant set for non-coding features is much smaller than the tree.
2. **UTR sits between CDS and intron**, closer to CDS in mammals (0.88
   mean) and to intron in flies (0.31), consistent with UTRs being shorter
   and more constrained in compact genomes. The asymmetry matters for the
   UTR objective: in flies the UTR signal will have to come from the
   sequence and from RNA-seq, not from the alignment.
3. **Multiz is reference-anchored and gappy between blocks**: the fly window
   came as 352 blocks and the human 470-way as 241, and per-informant
   coverage below 1.0 in CDS of a close species (cow, 0.986) is block
   boundaries, not biology. Cactus (241-way) gives 2,064 blocks for the same
   3.9 kb, because it records every rearrangement breakpoint. A model that
   consumes these needs an explicit "unaligned" state per informant per
   position, distinct from a gap, which is also what CONTRAST and N-SCAN
   documented (stalin, message 20260909T053948Z-stalin-0005). The
   Armstrong et al. 2020 Cactus paper reports that Cactus aligns more
   non-coding sequence than multiz at the same depth; we have not measured
   that difference ourselves beyond the block counts above.
4. **Ancestral sequences** exist only in the Ensembl EPO sets. If the design
   (T-human-011) wants internal-node states as inputs, the training data for
   them is limited to vertebrates and to the EPO species sets.

To turn the anecdote into a statistic, run `fetch_window.py` over the
benchmark's gene sample for each species that has an alignment and
aggregate the coverage tables by tree distance; that is a laptop job (the
fly window took 10 requests and 4.6 MB) and is listed in section 9.

## 7. Suitability for training versus held-out evaluation

The benchmark's rules (`docs/benchmark.md` section 3.2): held-out labels may
never be seen; held-out sequence in a pretraining corpus must be declared;
alignments used for *training* must be rebuilt with held-out species
removed; an informant set used at *inference* on a held-out target is
permitted and must be declared.

Applying them to what exists:

| Panel species | Split | Public multiple alignment referenced on it | Conservation | Fit |
|---|---|---|---|---|
| human | heldout_paired | hg38 100/30/470-way, Cactus 241/447, Ensembl mammals/primates/amniotes | phyloP, phastCons | evaluation only; any of these is a declarable informant set at inference. Not usable to train, because training on human windows means training on human labels |
| mouse | train | mm39 35-way; Ensembl 44/92-mammal EPO and 22-murinae re-referenced on mouse | mm39 phyloP35way | training source, but every one of these alignments contains human (held-out) sequence: sequence leakage, declare it, or drop the human row at cut time |
| chicken | heldout | Ensembl sauropsids EPO on GRCg7b (exact panel assembly); galGal6 77-way (older assembly, liftover needed) | galGal6 phyloP77way only | evaluation with a declared informant set; the sauropsid set contains no other panel species |
| zebrafish | train | Ensembl fish EPO on GRCz11 (panel is GRCz12ab) | none | training only after lifting coordinates, or by evaluating on GRCz11 instead; fish set contains fugu (held-out): drop that row when training |
| fugu | heldout_paired | Ensembl fish EPO on fTakRub1.2 (panel 1.3) | none | evaluation; informant set contains zebrafish (train), which the rules allow at inference if declared |
| fruit fly | train | dm6 124-way and 27-way | dm6 phyloP124way | training; the 124-way contains honey bee (`apiMel4`, held-out): drop that row |
| C. elegans | train | ce11 135-way (raw MAF only) | ce11 phyloP135way | training; no panel species among informants |
| yeast | train | sacCer3 7-way | phastCons7way | training; near-intronless, so mostly a negative control for the intron machinery |
| frog, honey bee, sea anemone, ciona, thale cress, rice, maize, S. pombe, Neurospora, Dictyostelium, Plasmodium, Tetrahymena | 5 train, 7 held-out | none (rice: 8-way EPO on IRGSP, not the panel's AGIS1.0) | none | no comparative input available; the model runs single-genome on these unless we build alignments |

Two consequences for the design:

- **Every clade-general claim will be tested without informants** on seven
  of the ten held-out species, whatever the model. The comparative branch of the
  model must be optional at inference, and the benchmark must report the
  with- and without-informant numbers separately (a `question` for lenin on
  whether `benchmark/report.py` should carry that column is in section 9).
- **Rebuilding training alignments without held-out species** is cheap when
  the alignment is reference-anchored to a training species: dropping rows
  from a MAF at cut time is a filter, not a realignment, because multiz
  rows are independent given the reference. It is not cheap for Cactus and
  EPO, whose blocks and ancestral sequences were inferred jointly with the
  held-out species; there, dropping a row leaves an alignment that still
  benefited from it. For an honest held-out test on those, the informant
  set at inference is the right tool, and it is allowed.

## 8. What we would have to build, with costs

For the twelve species without a usable public alignment, the options are:

1. **No informants**: zero cost; tests the single-genome half of the model.
2. **Pairwise informant alignments** (lastz/chain/net, or minimap2 for
   close pairs) between each such species and two to five relatives with
   RefSeq annotation, combined by the reference-anchored multiz convention.
   Cost is dominated by lastz, which scales with the product of genome
   sizes; maize (2.2 Gb) against a relative is the expensive case. This
   needs measurement before it can be costed honestly; the T-human-009
   owner should have the numbers for a single lastz genome pair, and I will
   ask rather than guess.
3. **Progressive Cactus** on a small clade per species (5 to 10 genomes).
   Armstrong et al. 2020 report Cactus runtimes for hundreds of genomes on
   a cluster; for 5 to 10 genomes under 500 Mb it is a multi-core
   laptop-scale job of hours, but a 2 Gb plant clade is not. Under the
   charter this is an `alert` with an estimate, not something to run.

None of this is needed for the Phase 1 to 3 deliverables. It is a Phase 4
cost, and section 7's first consequence (the model must run without
informants) is the design-level answer.

## 9. Open items and questions for other owners

1. Aggregate section 6 over the benchmark gene sample (laptop job, next
   tick).
2. Decide the window-cutting convention for training examples: fixed
   reference length with stride, informant rows reduced to per-species
   one-hot plus an "unaligned" channel, ancestral rows optional. Draft in
   the next tick; it belongs to this task's `scripts/data/`.
3. Question for lenin (T-human-007): should `benchmark/report.py` carry a
   "with informants / without informants" column, given that seven of the
   ten held-out species have no public alignment?
4. Question for whoever takes T-human-009: a measured lastz wall-clock and
   memory for one genome pair at 100 Mb and at 2 Gb, to cost section 8
   option 2.
5. The Ensembl licence paragraph and the 470-way species composition page
   did not render to the fetcher; both are cited by URL and should be
   re-read.
6. GRCz12ab (zebrafish) has no GenArk hub and no Ensembl alignment yet; the
   panel may be better served by GRCz11 for zebrafish until the new assembly
   propagates, which is lenin's call.
