# Data sources: alignments, conservation, expression and annotations

Task T-human-008. Author: marx. Status: draft 2, 2026-09-09 (for review).

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
| Zoonomia / Cactus consortia (**excluded for now**, decision 20260910T004401Z-human-0009: too fragmented to rely on; section 5) | 241-mammal Cactus HAL (2020) and 447-way (2023) alignments; phyloP from them | human (and any mammal in the HAL, via `halLiftover`) | https (UCSC CGL), UCSC tracks | open; HAL tools needed for anything but the hg38-referenced bigMaf |
| Community databases in the panel | FlyBase, WormBase, TAIR/Araport, SGD, PomBase annotations; VEuPathDB RNA-seq for Plasmodium | as listed in `panel.tsv` | through NCBI RefSeq mirrors and GenArk `contrib` tracks | per-database, all open |

The headline: **public multiple alignments exist only for vertebrates,
flies, nematodes and yeasts, and only on a few reference assemblies.** For
thirteen of the twenty panel species there is no multiple alignment that
can be used as it stands: twelve (frog, honey bee, sea anemone, ciona,
thale cress, rice, maize, fission yeast, Neurospora, slime mould,
Plasmodium, Tetrahymena) have no public multiple alignment at all, and
zebrafish, a training species, has only the Ensembl fish EPO, which is on
the older GRCz11 and was inferred jointly with the held-out fugu, so under
`docs/benchmark.md` section 3.2 it is rebuild-only whichever assembly the
panel uses (the panel stays on GRCz12ab; `docs/benchmark.md` section 2.4,
decided in note 20260909T101258Z-lenin-0012). Fugu's alignment is on an
older assembly than the panel's but is usable as a declared informant set
at inference. Any comparative model trained or evaluated on those thirteen
must either run without informants or on alignments we build ourselves.
Section 8 costs this.

## 2. UCSC Genome Browser

### 2.1 Multiple alignments

Listing of `hgdownload.soe.ucsc.edu/goldenPath/<db>/` and `/gbdb/<db>/` on
2026-09-09. "Compressed" is the `maf/*.maf.gz` directory total; "raw" is the
uncompressed per-chromosome MAF that backs the browser track and is served
under `/gbdb/` with HTTP `Accept-Ranges: bytes` (verified with HEAD). The
470-way is the exception: its raw MAF sits under
`goldenPath/hg38/multiz470way/maf/` (499 uncompressed files dated 2022-08,
5.9 TB in total, 5.7 TB in the 25 primary chromosomes, `chr1.maf` alone
471 GB; listing read 2026-09-09), and there is no compressed copy.

| Assembly | Track | Aligner | Year | Species | Format on hgdownload | Compressed | Raw | Panel relevance |
|---|---|---|---|---|---|---|---|---|
| hg38 | multiz100way | multiz | 2015 | 100 | maf.gz per chromosome + `.nh` | 71.4 GB (357 files) | 790 GB | human `heldout_paired` |
| hg38 | multiz30way | multiz | 2017 | 30 | maf.gz + `.nh` | 18.0 GB | 156 GB | human |
| hg38 | multiz470way | multiz | 2022 | 470 (mammals only; 431 distinct species) | bigMaf (API) + uncompressed `maf/` per chromosome + `.nh` | none | 5.9 TB (499 files) | human; largest mammal set |
| hg38 | cactus241way (`cactus241wayBM`) | Cactus (Zoonomia) | 2020 | 241 mammals | bigMaf (API) + `.nh`; phyloP bigWig 9.6 GB | n/a | n/a | human; **excluded for now** (decision 20260910T004401Z-human-0009) |
| hg38 | cactus447way | Cactus (Zoonomia + primates) | 2023 | 447 | bigMaf (API) + `.nh.txt`; phyloP bigWig 10.0 GB | n/a | n/a | human; **excluded for now** (decision 20260910T004401Z-human-0009) |
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
- Composition of the 470-way (`scripts/data/tree_composition.py` on
  `hg38.470way.scientificNames.nh`, orders from GBIF on 2026-09-09 in
  `scripts/data/hg38.470way.orders.tsv`). It is a mammal-only alignment:
  470 leaves, 444 distinct leaf names, 431 distinct binomials; 35 species
  are present as two or more assemblies or subspecies (two dog assemblies
  plus the dingo, two *Mus musculus* and two *Rattus norvegicus*
  assemblies, two jaguars, and so on), so an informant axis of 469 rows carries about 8% duplicate signal.
  Distance is patristic from human on the track's own tree, in
  substitutions per site; "within 0.5" and "beyond 1" are the horizon
  thresholds of section 6.2. Orders with fewer than five leaves are
  summed in the last row; the script prints all 23.

  | order | leaves | species | min | median | max | within 0.5 | beyond 1 |
  |---|---|---|---|---|---|---|---|
  | Artiodactyla (incl. whales) | 126 | 113 | 0.323 | 0.408 | 0.442 | 126 | 0 |
  | Rodentia | 87 | 84 | 0.325 | 0.467 | 0.544 | 54 | 0 |
  | Carnivora | 68 | 58 | 0.339 | 0.360 | 0.400 | 68 | 0 |
  | Primates | 60 | 58 | 0.012 | 0.110 | 0.280 | 60 | 0 |
  | Chiroptera | 52 | 49 | 0.345 | 0.392 | 0.450 | 52 | 0 |
  | Diprotodontia (marsupial) | 13 | 13 | 0.784 | 0.812 | 0.822 | 0 | 0 |
  | Perissodactyla | 13 | 10 | 0.295 | 0.299 | 0.308 | 13 | 0 |
  | Eulipotyphla | 7 | 7 | 0.384 | 0.426 | 0.552 | 6 | 0 |
  | Lagomorpha | 7 | 5 | 0.370 | 0.377 | 0.451 | 7 | 0 |
  | Pholidota | 6 | 4 | 0.341 | 0.345 | 0.348 | 6 | 0 |
  | Pilosa | 5 | 4 | 0.334 | 0.342 | 0.364 | 5 | 0 |
  | 12 other orders (incl. 6 more marsupials, 2 monotremes) | 25 | 25 | 0.236 | 0.426 | 1.023 | 16 | 2 |
  | total informants | 469 | 430 | 0.012 | 0.377 | 1.023 | 413 | 2 |

  Histogram of the 469 informants by distance from human: 28 within 0.1,
  26 in 0.1 to 0.25, 359 in 0.25 to 0.5, 54 in 0.5 to 1 (rodents past
  0.5, marsupials), 2 beyond 1 (platypus, echidna), none beyond 2. The
  same script on `hg38.100way.nh` reproduces the section 6.1 human
  histogram (9, 2, 40, 14, 24, 10) and puts the mouse at 0.502 against
  0.547 on the 470-way tree, dog 0.332 against 0.371, opossum 0.766
  against 0.810, so distances on the two trees are comparable to about
  10%.
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
predictions against RefSeq references. CDS completeness: the RefSeq
genePred tracks state it in `cdsStartStat` / `cdsEndStat` and every
CDS-bearing transcript checked on four assemblies is `cmpl` at both ends;
the GENCODE `knownGene` bigGenePred sets those columns to `none`
throughout and states truncation only through its `tag` column
(`cds_start_NF`, `cds_end_NF`) and the frame of the first coding exon in
`exonFrames`; the fetcher reads all three (section 6.3). CDS
intervals are the exons clipped to `cdsStart` / `cdsEnd` (`thickStart` /
`thickEnd` on a bigGenePred), read from the record and never assumed to
coincide with the exon bounds: stalin found a GenePred reader in
Vipsania's evaluation dependency that builds CDS blocks from exon starts
and ends alone, so every UTR base scores as coding (relay note
20260909T204115Z-stalin-0020), and `fetch_window.py --self-test` now holds
that synthetic case (exons [100,200) and [300,400) with coding bounds
[130,370) give CDS [130,200) and [300,370), on both strands and in both
record shapes, equal bounds give no CDS, a bound inside the intron leaves
the far exon UTR) together with the stat-column, NF-tag and frame reading.

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

Licence (re-read 2026-09-09 from `www.ensembl.org/info/about/legal/disclaimer.html`,
release 116 page): "Ensembl imposes no restrictions on access to, or use of,
the data provided and the software used to analyse and present it. Ensembl
data generated by members of the project are available without
restriction. Ensembl code written by members of the project is provided
under the Apache 2.0 licence. Some of the data and software included in the
distribution may be subject to third-party constraints." The third-party
clause matters for the alignments, whose input genomes come from many
submitters; none of the panel assemblies carries a use restriction in its
NCBI record.

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

## 5. Zoonomia and the large Cactus alignments (excluded for now)

**Status.** Excluded from the inventory by coordinator decision 20260910T004401Z-human-0009
(review of T-human-008): the Zoonomia Cactus alignments are too fragmented
to rely on for now. The HBB measurement in section 6.3 is the evidence on
record (2,517 blocks over 4,931 bp, a median block of about 2 bp, and
duplicated rows in 2,513 of the blocks) and is kept there as a finding
made on an alignment we are not using. No script under `scripts/data/`
defaults to a Cactus track (`fetch_window.py` defaults to `multiz470way`
on hg38 and otherwise to the first multiz `*way` track); `cut_windows.py`
still classifies the `cactus*` names as jointly inferred so that a window
fetched from one on purpose is handled correctly. The paragraphs below
record what exists and how to reach it, for when the exclusion is
revisited.

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

### 6.1 Aggregated over a gene sample, by tree distance

Run on 2026-09-09 with `scripts/data/coverage_by_distance.py` over windows
fetched by `fetch_window.py` (flank 500, every transcript in the window
labelled as above). Gene sample: a fixed-seed draw from `ncbiRefSeqCurated`
of NM_ transcripts with a complete CDS and at least three exons, one per
gene, non-overlapping, on one chromosome per species (dm6 chr2L, 2 to 8 kb;
hg38 chr11, 3 to 12 kb, so the human sample is short genes, and its
intron column is short introns). Distance is the patristic distance from
the reference to the informant on the track's own tree, in substitutions
per site. "cov" is base-weighted coverage over all windows and informants
in the bin; ">0.5" is the share of informant-window pairs with more than
half the bases of that class aligned.

Fruit fly, dm6 multiz124way, 12 windows (`Idgf3`, `Fs(2)Ket`, `CG9338`,
`Cyp28d2`, `Pi3K21B`, `fon`, `mus201`, `piwi`, `Tnpo-SR`, `Mal-B2`, `Aldh`,
`ChLD3`), 57,216 reference bases (24,354 CDS, 13,097 UTR, 15,059 intron,
4,706 intergenic), 123 informants, 65.6 MB in 141 requests:

| distance (subst/site) | informants | CDS cov | UTR cov | intron cov | intergenic cov | CDS >0.5 | UTR >0.5 | intron >0.5 | intergenic >0.5 |
|---|---|---|---|---|---|---|---|---|---|
| 0-0.1 | 1 | 0.999 | 0.978 | 0.942 | 0.948 | 1.00 | 1.00 | 1.00 | 1.00 |
| 0.1-0.25 | 3 | 0.997 | 0.941 | 0.862 | 0.915 | 1.00 | 1.00 | 1.00 | 1.00 |
| 0.5-1 | 9 | 0.988 | 0.873 | 0.775 | 0.734 | 1.00 | 0.95 | 0.93 | 0.84 |
| 1-2 | 25 | 0.881 | 0.445 | 0.320 | 0.357 | 0.92 | 0.43 | 0.26 | 0.34 |
| >2 | 85 | 0.732 | 0.051 | 0.063 | 0.044 | 0.79 | 0.01 | 0.01 | 0.01 |

Human, hg38 multiz100way, 10 windows (`KBTBD4`, `NUDT22`, `TP53AIP1`,
`TRIM77`, `NKAPD1`, `ACP2`, `GPR137`, `HMBS`, `FKBP2`, `APBB1`), 82,125
reference bases (13,006 CDS, 18,810 UTR, 46,053 intron, 4,256 intergenic),
99 informants, 57.0 MB in 117 requests:

| distance (subst/site) | informants | CDS cov | UTR cov | intron cov | intergenic cov | CDS >0.5 | UTR >0.5 | intron >0.5 | intergenic >0.5 |
|---|---|---|---|---|---|---|---|---|---|
| 0-0.1 | 9 | 0.976 | 0.970 | 0.939 | 0.954 | 1.00 | 0.99 | 1.00 | 0.98 |
| 0.1-0.25 | 2 | 0.910 | 0.878 | 0.723 | 0.729 | 0.95 | 1.00 | 0.85 | 0.83 |
| 0.25-0.5 | 40 | 0.889 | 0.737 | 0.561 | 0.560 | 0.90 | 0.93 | 0.59 | 0.76 |
| 0.5-1 | 14 | 0.801 | 0.420 | 0.270 | 0.279 | 0.84 | 0.58 | 0.24 | 0.36 |
| 1-2 | 24 | 0.556 | 0.052 | 0.057 | 0.031 | 0.62 | 0.02 | 0.01 | 0.03 |
| >2 | 10 | 0.478 | 0.015 | 0.037 | 0.008 | 0.49 | 0.00 | 0.00 | 0.01 |

Mouse, mm39 multiz35way, 10 windows (`Slc22a6`, `Lipo1`, `Rcor2`, `Doc2g`,
`Tlx1`, `Ap5b1`, `Klc2`, `Slc22a12`, `Kcnk4`, `Trim8`; drawn on chr19,
3 to 12 kb, with `scripts/data/sample_genes.py`, seed 20260909), 90,270
reference bases (15,599 CDS, 13,242 UTR, 54,549 intron, 6,880
intergenic), 34 informants, 17.4 MB in 76 requests. The 35-way has one
informant under 0.25 (rat) and four in the 0.25 to 0.5 band (a GenArk
rodent assembly at 0.26, beaver, squirrel, colugo); 24 of the 34, primates
and laurasiatheres included, sit between 0.51 and 0.76:

| distance (subst/site) | informants | CDS cov | UTR cov | intron cov | intergenic cov | CDS >0.5 | UTR >0.5 | intron >0.5 | intergenic >0.5 |
|---|---|---|---|---|---|---|---|---|---|
| 0.1-0.25 | 1 | 1.000 | 0.966 | 0.906 | 0.948 | 1.00 | 1.00 | 1.00 | 1.00 |
| 0.25-0.5 | 4 | 0.945 | 0.827 | 0.596 | 0.666 | 0.95 | 0.90 | 0.82 | 0.78 |
| 0.5-1 | 24 | 0.894 | 0.653 | 0.443 | 0.451 | 0.91 | 0.68 | 0.59 | 0.51 |
| 1-2 | 3 | 0.654 | 0.079 | 0.059 | 0.020 | 0.77 | 0.07 | 0.03 | 0.00 |
| >2 | 2 | 0.484 | 0.004 | 0.018 | 0.000 | 0.55 | 0.00 | 0.00 | 0.00 |

Human seen from mouse (0.52 on this tree; 0.50 seen from human on the
100-way) aligns 0.65 of mouse intron bases in these windows; the 100-way
puts mouse in its 0.25 to 0.5 bin, whose intron coverage over all 40
informants is 0.56, so the pair looks about the same from either end.

Worm, ce11 multiz135way, 10 windows (`Y45F3A.4`, `clp-2`, `C40H1.7`,
`C50C3.2`, `Y71H2AM.3`, `affl-1`, `ZK1128.7`, `enu-3.2`, `B0524.2`,
`ttm-1`; drawn on chrIII, 2 to 8 kb, with `scripts/data/sample_genes.py`,
seed 20260909), 50,814 reference bases (18,803 CDS, 4,867 UTR, 20,321
intron, 6,823 intergenic), 134 informants, 30.2 MB in 98 requests. The
track's own tree puts the nearest informant (*Caenorhabditis* sp. 11) at
1.11 substitutions per site from *C. elegans* and the median informant at
2.54, so the first three bins are empty:

| distance (subst/site) | informants | CDS cov | UTR cov | intron cov | intergenic cov | CDS >0.5 | UTR >0.5 | intron >0.5 | intergenic >0.5 |
|---|---|---|---|---|---|---|---|---|---|
| 1-2 | 34 | 0.500 | 0.219 | 0.200 | 0.174 | 0.53 | 0.24 | 0.12 | 0.18 |
| >2 | 100 | 0.119 | 0.047 | 0.052 | 0.041 | 0.10 | 0.00 | 0.00 | 0.01 |

The 1 to 2 bin is two different things: the 18 *Caenorhabditis* informants
(1.11 to 1.39) have per-informant intron coverage of 0.16 to 0.46 (mean
0.34) and CDS 0.59 to 0.90, while the 16 non-*Caenorhabditis* rhabditids
in the same bin (*Pristionchus*, *Oscheius*, *Diploscapter*, strongylids,
1.68 to 1.90) have intron coverage 0.00 to 0.18 and CDS mostly below 0.5.

The per-pair rows (1,476 fly, 990 human, 340 mouse, 1,340 worm) are what
`--tsv` writes; they are not committed (regenerate with the commands in
`scripts/data/README.md`).

### 6.2 What the aggregate adds to the anecdote

1. **The non-coding signal has a horizon near one substitution per site
   on the three vertebrate and insect trees, and the worm tree is scaled
   differently.** Below 0.5 substitutions per site introns and
   intergenic sequence are mostly aligned (fly 0.78 to 0.94, human 0.56 to
   0.94); between 0.5 and 1 they are aligned in a quarter to three quarters
   of pairs; past 1 they are essentially gone (under 0.06 base-weighted,
   and only 1 to 3% of pairs above half). CDS keeps going: 0.73 (fly) and
   0.48 (human) even beyond 2 substitutions per site. The two references
   agree on the shape although their trees are different (the fly tree
   reaches 2 substitutions per site inside Diptera; the human tree reaches
   it at fish). Mouse (35-way) gives the same shape from the third
   reference: the 0.25 to 0.5 band is at 0.60 intron
   coverage, 0.5 to 1 at 0.44, opossum (1.0) at 0.18, and chicken, frog,
   fish and lamprey (1.3 to 2.5) at 0.06 or less while CDS stays at 0.46
   to 0.65. The worm 135-way does not: its
   nearest informant is 1.11 substitutions per site away, yet the 18
   *Caenorhabditis* informants between 1.11 and 1.39 still align 0.34 of
   intron bases on average, where the fly and mammal trees give under 0.06.
   Either the multiz tree for ce11 is on a different scale (its root-to-
   *C. elegans* branch alone is 0.43) or *Caenorhabditis* non-coding
   sequence is more alignable per substitution; the point for the design
   is that the horizon should be learned per tree or expressed through the
   coverage the cutter records, not hard-coded as a distance.
2. **Most of the tree is beyond the horizon.** Of the fly 124-way's 123
   informants, 110 sit past 1 substitution per site; of the human 100-way's
   99, 34 do, and 40 of the rest are in the 0.25 to 0.5 band (laurasiatheres,
   glires). So the informant rows that carry an intron or UTR signal are
   about 13 of 123 for fly and about 65 of 99 for human; for mouse it is
   29 of 34 (the 35-way was built for a mouse browser and stops at
   lamprey), and for worm it is 0 of 134 by distance and 18 by the
   *Caenorhabditis* coverage of item 5. The 470-way is
   the opposite case: it is mammal-only, 413 of its 469 informants are
   within 0.5 substitutions per site of human and only the two monotremes
   are beyond 1 (section 2.1 composition table), so at `K` = 469 nearly
   every row is inside the horizon and the axis is dense, with 359 rows
   packed into the 0.25 to 0.5 band and about 8% of rows duplicating a
   species. So the informant axis has to handle both regimes: mostly
   `unaligned` rows on the vertebrate-wide and fly trees, and hundreds of
   near-redundant aligned rows on the mammal tree. That is an argument for
   the variable-`K` option in section 6.3 or for a sparse attention over
   aligned rows, and, on the dense side, for subsampling or clade-pooling
   rows at cut time rather than paying attention cost for 126 artiodactyls.
3. **UTR tracks intron, not CDS, once distance grows.** In the human 0.5 to
   1 band UTR is at 0.42 against intron 0.27 and CDS 0.80; in the fly 1 to 2
   band UTR 0.45, intron 0.32, CDS 0.88. So the UTR objective will get its
   comparative signal from the same close informants as the intron
   objective, plus sequence and RNA-seq.
4. **Within the horizon, the human close band is denser than the fly close
   band** because the 100-way has nine primates within 0.1 substitutions
   per site and the 124-way has one *simulans*-group species; that is a
   property of which genomes were sequenced, not of the clades, and it is
   why per-species alignment depth in Table 7 (section 7) matters as much
   as the existence of an alignment.
5. **The worm training alignment has no informant inside the horizon.**
   On the ce11 135-way every informant is beyond 1 substitution per site
   and 100 of 134 are beyond 2; the whole *Caenorhabditis* signal is
   carried by 18 rows with a third of intron bases aligned. So for the one
   nematode in the panel, the model must get most of its intron and UTR
   evidence from sequence and RNA-seq, and the comparative channel is
   CDS-dominated (0.50 base-weighted in the 1 to 2 bin). That is the
   opposite of the mammal case in item 2, and it argues for training with
   informant dropout so the model does not learn to depend on close rows.
6. Caveats: one chromosome per species, short human genes, twelve, ten,
   ten and ten windows; multiz block boundaries depress close-informant CDS coverage a
   little (the 0.976 for primates is boundaries, not biology); the fly
   intergenic column is 4.7 kb of flank. The direction of every effect is
   the same as in the single-locus tables above, which is what the sample
   was for.

### 6.3 Window-cutting convention for training examples

`scripts/data/cut_windows.py` is the convention; this section states it so
the design task (T-human-011) can argue with it. Every fetched locus is cut
into examples of fixed reference length `L` (default 4,096) with stride `S`
(default 2,048), on the reference + strand, one entry per reference base:

- **Reference**: base code (A, C, G, T, N) and a padding mask.
- **Labels**: per base, intergenic / intron / UTR exon / CDS (union over
  overlapping RefSeq transcripts, CDS taking precedence, the same rule as
  the coverage tables); CDS frame 0 to 2 counted from the start codon on the
  transcript's strand; the strand of the transcript that owns the base; and
  boundary marks at the first base of the start codon, the last base of the
  stop codon (RefSeq includes the stop in the CDS), the first intron base
  (donor) and the last intron base (acceptor). Boundaries are placed at
  their + strand coordinate; the strand channel says which way to read them.
  The owner of a base is the transcript giving it its highest label class;
  among transcripts tied at that class the shortest span wins, then
  annotation order, and frame is the owner's. The tie rule matters for
  genes nested in another gene's intron, which are common in the fly and
  in large vertebrate loci: `Adh` sits inside a minus-strand gene that
  spans the whole demonstration window, and before this rule every base of
  the window, including `Adh`'s 3,331 CDS bases across the four examples,
  carried the enclosing gene's strand while its frames and boundary marks
  were counted on the plus strand.
  **Isoforms.** By default the label is a union over the isoforms of a
  locus, so an example's target is not the structure of any one transcript
  where isoforms differ; a model trained on it learns the union, and the
  benchmark's transcript-level column then measures the reference's isoform
  density rather than the model, as lenin observed for Helixer on fugu
  (note 20260909T133356Z-lenin-0013). The benchmark deliberately does not
  name a representative isoform: its section 4.4 matches whatever a model
  emits against the best-fitting annotated isoform, so union, longest CDS
  and a curated list such as MANE Select are all scoreable, and which one
  a model trains on is a training decision, not a scoring one (note
  20260909T141456Z-lenin-0014). The cutter therefore offers the choice as
  `--isoforms`: `union` (default), `longest-cds` (one transcript per
  locus, the most CDS bases, ties to the longest exonic span, then
  annotation order; a locus with no coding transcript keeps its longest
  exonic span), or `representative` with `--representatives FILE` (ids
  one per line or the first column of a TSV, compared with and without a
  version suffix; a locus in which no listed transcript occurs falls back
  to `longest-cds` and the sidecar names it under
  `representative_fallback_loci`). A locus is the stated gene (`name2` /
  `geneName` from UCSC genePred, the parent gene's `external_name` from
  Ensembl); transcripts with no gene name are clustered by overlapping span
  on the same strand. The sidecar records `isoform_policy` and every
  transcript dropped with its locus and reason, so a training set built
  under one policy cannot be mistaken for another. What the choice does
  on the fly `Adh` window (chr2L:14615052-14619402, 4,350 bases, 12 RefSeq
  transcripts at three loci, cut as one example with `apiMel4` dropped,
  re-fetched 2026-09-09): the six `Adh` isoforms and the three `Adhr`
  isoforms each share one CDS and differ only in 5' UTR exons, so
  `longest-cds` decides them on the tie rule (longest exonic span:
  `NM_001032098.2` for `Adh`, `NM_058147.4` for `Adhr`) and keeps
  `NM_001169504.2` for the enclosing `osp`. The label counts move from
  1,590 / 664 / 2,096 (CDS / UTR / intron) to 1,590 / 490 / 2,270: the 174
  bases that only a dropped isoform's UTR exon covered become intron, one
  boundary mark of a dropped first exon disappears, and no CDS, strand or
  frame value changes. What a policy discards deserves more than a
  boundary-mark count, because the splice sites a reference states are not
  shared equally between donors and acceptors: lenin counted 257,153
  distinct CDS acceptors against 212,078 donors in GENCODE 50 and 193,359
  against 188,913 in human RefSeq (note 20260909T152910Z-lenin-0015), so
  a one-isoform target throws away acceptors first. The sidecar therefore
  records `distinct_sites`: the distinct start codons, stop codons, donors
  and acceptors the annotation states over all isoforms in the example,
  the same counted over the painted isoforms, and the difference, with a
  donor or acceptor whose intron lies between two CDS segments counted
  again as `cds_donor` / `cds_acceptor`. On `Adh` the cost is one UTR
  acceptor: 5 donors and 6 acceptors over all isoforms, 5 and 5 painted,
  and all four CDS donors and four CDS acceptors kept. On human `TP53`
  (chr17:7667921-7687990 on the hg38 100-way, 28 RefSeq transcripts at two
  loci with the antisense `WRAP53`, fetched 2026-09-09, 9.1 MB in 14
  requests) the annotation states 14 donors, 14 acceptors, 4 starts and 3
  stops, of which 9 donors and 11 acceptors are CDS-internal;
  `longest-cds` keeps 12 donors, 11 acceptors, 1 start and 1 stop, so it
  drops 2 CDS acceptors (alternative CDS exons that only shorter-CDS
  isoforms use), 3 alternative starts, 2 stops and no CDS donor, and the
  CDS label loses 81 bases (1,263 to 1,182). It also does not pick the
  MANE Select transcript: five `TP53` transcripts including `NM_000546.6`
  tie at 1,182 CDS bases with the same CDS, and the tie goes to the
  longest exonic span, `NM_001407262.1`, which has one more 5' UTR exon.
  So `longest-cds` keeps the curated protein but not necessarily the
  curated transcript; on human and mouse, `representative` with the MANE
  list is the policy that keeps both. That is the general shape of the
  cost in RefSeq: alternative UTR exons are far more common than
  alternative CDS, so the policy mostly decides what a model is told about
  UTRs, which the transcript-level score does not look at; where isoforms
  do differ in CDS the loss falls on acceptors and starts, and the sidecar
  counts it per example so that T-human-011 can sum it over a training
  set rather than estimate it. MANE Select exists for human and mouse only
  (UCSC `mane` track, NCBI RefSeq `MANE` tag); for the other panel species
  `representative` needs a list the operator supplies, and `longest-cds`
  is the policy that works everywhere. A second caveat: where two
  transcripts of the same class overlap in different frames (six bases
  between isoforms of the gene enclosing `Adh`, all outside the fetched
  window), the frame channel records one of them by the same tie rule.
  **Incomplete CDS ends.** lenin's partial-gene work on the benchmark
  (note 20260909T171958Z-lenin-0016) found that 5.3% of GENCODE 50's
  CDS-bearing transcripts have an incomplete CDS end which the GTF states
  only by omitting the `start_codon` or `stop_codon` feature, with the
  `cds_start_NF` / `cds_end_NF` tags agreeing at 99%; a cutter that paints
  a start and a stop at every CDS end would therefore invent about 13,000
  starts and 20,000 stops the annotation never claimed. What the sources
  this fetcher reads actually declare (all read 2026-09-09 through the
  UCSC API): the RefSeq genePred tracks carry `cdsStartStat` /
  `cdsEndStat`, and every CDS-bearing transcript of `ncbiRefSeqCurated` on
  hg38 chr21 (726), dm6 chr2L (5,707), ce11 chrIII (3,657) and mm39 chr19
  (1,435), and of the model-inclusive `ncbiRefSeq` on hg38 chr21 (1,373),
  is `cmpl` at both ends, so for RefSeq the problem does not arise in the
  data the benchmark uses as truth; the GENCODE `knownGene` bigGenePred on
  hg38 and mm39 sets both columns to `none` for every transcript (3,786 on
  chr21, 37,203 on chr1, 9,227 on mm39 chr19) but keeps GENCODE's tags in
  its `tag` column and the frame of every coding exon in `exonFrames`; and
  Ensembl's REST overlap endpoint states nothing about CDS ends. The
  fetcher now records what a UCSC record declares (`cds_start_status`,
  `cds_end_status` from the stat columns, then from the NF tags, and
  `cds_start_frame` from the first coding exon in transcript orientation),
  and the cutter's `--partial-ends` decides each end: `both` (default)
  takes the declaration where there is one and otherwise reads the
  reference sequence, `declared` and `sequence` use one signal, `none`
  paints every end as before. The sequence rule is that a CDS is 5'
  incomplete when its first codon is not ATG and 3' incomplete when
  neither its last codon nor the codon after it is a stop under
  `--genetic-code` (NCBI table numbers; 6 for *T. thermophila*, where TAA
  and TAG read glutamine, or `--stop-codons` explicitly). Checked against
  GENCODE's own tags over the 3,786 CDS-bearing transcripts of hg38 chr21:
  3,664 of the 3,669 untagged transcripts open on ATG and none of the 117
  `cds_start_NF` ones do (15 open on a near-cognate codon, 102 on
  something else), all 3,637 untagged transcripts end on a stop while 136
  of the 149 `cds_end_NF` ones have none (7 have one in the following
  codon and 6 end on one), so the two rules agree with the tags on 99.87%
  and 99.66% of transcripts; every untagged transcript has a CDS length
  divisible by three and 169 of the 247 tagged ones do not; and 83 of the
  117 `cds_start_NF` transcripts have a first-coding-exon frame of 1 or 2
  where every untagged one has 0. The five untagged non-ATG starts
  (`NCAM2`, `DSCAM`, three `MORC3` nonsense-mediated-decay isoforms) are
  what `sequence` would mark incomplete that GENCODE does not; the
  sidecar's `cds_ends` lists every such disagreement between declaration
  and sequence, and under `both` the declaration wins. An incomplete end
  gets no start or stop mark in the boundary channel, is not counted in
  `distinct_sites`, and a declared frame of 1 or 2 shifts the frame
  channel so codon position is right from the first base. Measured on
  `ATP5PO` (chr21:33901026-33915853, 14,827 bases, 53 GENCODE transcripts
  of which 45 are coding, hg38 100-way, 15 requests, 9.8 MB, fetched
  2026-09-09 with `--annotation knownGene`): four transcripts tagged
  `cds_start_NF` (three with frame 2) and one `cds_end_NF`; the sequence
  agrees on all of them (40 ATG starts, 1 near-cognate, 3 other, 1 whose
  CDS start lies outside the window; 44 stops inside the CDS and the
  tagged 3' end outside the window; four CDS lengths not divisible by
  three; no disagreement), and the default cut paints one start mark and
  twelve stop marks where `none` paints five starts, with `distinct_sites`
  moving from 5 starts to 1. That window also shows why the declaration
  is read first: the tagged 3' end lies beyond the fetched window, so the
  sequence alone would have passed it as complete. For T-human-011: a
  model trained on RefSeq truth sees complete ends only; one trained on
  GENCODE or Ensembl sees about 5% of loci where one end is missing, and
  the cutter now tells it which rather than painting a codon there.
  Which transcripts paint is `--transcript-types`: the default follows the
  benchmark's truth rule (`docs/benchmark.md` section 4) and excludes
  pseudogenes and immunoglobulin / T-cell receptor segments, recognised only
  from the type the source states (Ensembl `biotype`, GenArk `geneType`;
  UCSC `ncbiRefSeq` genePred carries none, so everything there paints), so
  that the labels do not teach structure at loci the scorer deletes; `all`
  paints everything and the sidecar lists what was excluded.
  The sidecar's `feature_lengths` records, over all isoforms and over the
  painted ones, the lengths of the distinct introns with both ends inside
  the example (minimum, median, maximum, and the count of introns cut by
  the edge), the CDS length and the transcript span per transcript (read
  from the whole annotation record, with a count of those reaching past
  the edge), and how many of each lie under fixed floors: introns shorter
  than 30 and than 50 bases, CDSs shorter than 60, spans of at most 80.
  Those are the floors engels read out of Helixer's decoder (relay note
  20260909T202745Z-engels-0020: the HMM's shortest intron path is 30
  bases and its U2 GT-AG and GC-AG paths 50, the candidate gate at the
  default window and peak needs a genic run above 80 bases, and the
  published minimum coding length is 60), and any decoder with hard
  minimum durations has floors of the kind; the record says per example
  which labels such a model can never reproduce. On the two demonstration
  windows nothing lies under them: `Adh` has 6 distinct introns of 51 to
  659 bases (median 248; 40 intron ends of the enclosing gene lie outside
  the window), `TP53` 18 of 81 to 10,757 (median 704), and the shortest
  CDS is 771 and 549 bases. Two windows say nothing about how often the
  floors bite, so `scripts/data/length_floors.py` counts the same record
  over one whole chromosome of a RefSeq track, where nothing is clipped
  (coding transcripts only; introns distinct over isoforms; fetched
  2026-09-09 through the UCSC API, 0.2 to 2.6 MB per chromosome):

  | Assembly | Chromosome | Track | Coding transcripts | Single-exon | Distinct introns | Shortest | Median | Introns < 30 | Introns < 50 | CDS < 60 | Span < 81 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | hg38 | chr21 | ncbiRefSeqCurated | 726 | 56 (7.7%) | 2,382 | 74 | 2,248 | 0 | 0 | 0 | 0 |
  | mm39 | chr19 | ncbiRefSeqCurated | 1,435 | 81 (5.6%) | 6,940 | 65 | 1,175 | 0 | 0 | 0 | 0 |
  | dm6 | chr2L | ncbiRefSeqCurated | 5,707 | 677 (11.9%) | 10,204 | 43 | 92 | 0 | 48 (0.47%) | 1 | 0 |
  | ce11 | chrIII | ncbiRefSeqCurated | 3,657 | 84 (2.3%) | 15,468 | 1 | 110.5 | 12 (0.08%) | 3,940 (25.5%) | 9 (0.25%) | 12 (0.33%) |
  | sacCer3 | chrIV | ncbiRefSeq | 766 | 726 (94.8%) | 41 | 1 | 101 | 8 (19.5%) | 8 (19.5%) | 0 | 2 (0.26%) |
  | GCF_000002765.6 (*P. falciparum*) | NC_037283.1 | ncbiRefSeq | 773 | 364 (47.1%) | 1,294 | 65 | 142 | 0 | 0 | 0 | 0 |

  Three things follow. First, the U2 floor is a clade property: no
  human, mouse or *Plasmodium* intron on these chromosomes is under 50
  bases, half a percent of fly introns are, and a quarter of worm
  introns are (the worm distribution peaks at 47 bases, with 3,940
  distinct introns between 31 and 49 and only 12 shorter), so a decoder
  built with that floor cannot annotate *C. elegans* whatever its network
  learns, which is the intron-length-regime failure the charter names and
  a number the benchmark's intron-length decile (`docs/benchmark.md`)
  will show as a whole stratum lost. Second, the shortest "introns" are
  not introns: the 1-base gaps are how RefSeq encodes a programmed
  ribosomal frameshift inside a CDS (all eight on sacCer3 chrIV are Ty
  retrotransposon gag-pol ORFs, `YDR034C-D`, `YDR098C-B`, `YDR210W-B` and
  their kind; the two on ce11 chrIII sit in `cdh-4`), and until version
  0.8 the cutter painted them as a 1-base intron with a donor and an
  acceptor mark, which no splice model should be asked to learn. The ten
  worm introns of 15 to 29 bases all lie on the UTR side of a
  CDS end, where an annotation is least checked. Third, the single-exon
  fraction runs from 2% (worm) to 95% (yeast) with *Plasmodium* at 47%,
  so a candidate gate tuned on a spliced vertebrate genome sees a very
  different population of short genic runs in a compact one; the span
  floor itself bites on under 0.4% of transcripts everywhere here, so
  engels' gate is a corner case in these annotations rather than a
  stratum. The floors are constants in `cut_windows.py` (`LENGTH_FLOORS`);
  a different decoder's floors are one edit and a rerun.
  **The intron floor** (version 0.9). The benchmark had already drawn the
  line the frameshift gaps ask for: its scorer treats a gap between
  consecutive CDS blocks as an intron only if it is at least 20 bases
  (`benchmark/score.py` `MIN_INTRON`, commit `dbe9ffb`) and counts the
  shorter ones apart as `reference_cds_gaps_below_min` /
  `predicted_cds_gaps_below_min`, with the floor written into every
  result (lenin, note 20260909T211910Z-lenin-0019: 47 of the 343 CDS gaps
  in the *S. cerevisiae* reference are Ty frameshifts, and on fugu
  Tiberius itself emits 413 such gaps against 1,124 in the reference). The
  cutter now applies the same floor, `--min-intron` (default
  `MIN_INTRON = 20`): an exon gap shorter than it is painted with a fifth
  label class, `short_gap` (code 4, ranked above intron and below UTR when
  isoforms overlap), gets no donor or acceptor mark and no frame, and
  contributes a `short_gap` entry to `distinct_sites` instead of a donor
  and an acceptor; `feature_lengths.introns` still describes every exon
  gap the annotation states, with `below_min_intron` saying how many the
  floor removed, and a new sidecar record `short_gaps` lists each one
  (strand, coordinates, length and length mod 3, whether both flanks are
  CDS ends, the transcripts stating it, two exon bases on each side and
  the gap itself in transcript orientation) so that what the scorer leaves
  out of its splice denominator is not a splice site in the labels
  either, and is on record rather than dropped. Two things keep the floor
  a policy rather than a fact. stalin (note 20260909T214050Z-stalin-0021)
  found the counterexample: *Stentor coeruleus* has 15- and 16-base
  spliceosomal introns, 8,806 of them in one annotation (Nuadthaisong et
  al. 2022, doi:10.3390/ijms231810973; Slabodnick et al. 2017,
  doi:10.1016/j.cub.2016.12.057), which a floor of 20 files under
  `short_gap` wholesale; `--min-intron 15` restores them and 0 turns the
  floor off, and the sidecar carries the value used. And length does not
  name the event: NCBI's GFF3 marks a ribosomal-slippage CDS with
  `exception=ribosomal slippage` while its mRNA exon stays unsplit,
  but the UCSC genePred sources this fetcher reads split the exon and
  carry no exception, so `in_cds` and the flanks are what the record can
  say; the `short_gaps` entry has an `exception` field for a source that
  states one (Ensembl and NCBI GFF3 would; open item, section 9 item 12).
  engels (note 20260909T212826Z-engels-0021) showed that a motif-masked
  HMM decoder (bricks2marble's defaults under Tiberius: donor `NGT` or
  `NGC`, acceptor `AGN`, one intron state per frame with no self-loop
  required) admits a 1-base intron whose GT and AG windows each borrow a
  base from the flanking exon, `A|G|T`, so a strict motif mask by itself
  does not enforce an intron duration. The sidecar tests every short gap
  with that overlapping-window rule (`motif_window`, and
  `motif_borrows_exon_base` when a window reached into an exon) and with
  the gap's own ends (`motif_exact`, four bases or more). Because the
  combined verdict alone cannot say which window failed (engels, note
  20260909T222452Z-engels-0022), the entry also carries the donor and
  acceptor windows apart (`motif_donor`, `motif_acceptor`, each yes, no,
  or unresolved when the window holds a base outside A/C/G/T) and one
  `motif_class` per gap: `exact`, `borrows_exon_base`, `donor_only`,
  `acceptor_only`, `neither` (both windows decided and both fail) or
  `ambiguous` (a window left unresolved by a base outside A/C/G/T,
  whatever the other window says, so `neither` is reserved for two
  decided failures; engels, note 20260909T232517Z-engels-0023, which
  caught an earlier description saying the other window had not failed);
  the record's `motif_by_class` lists every class, and
  `motif_unresolved_windows` counts the windows an N left undecided, so a
  zero there says every tested window was resolved, which the class
  counts alone do not. It does not say the gap holds no N: the windows
  are two bases each, so in a gap of three bases or more an N can sit
  between them (`A|GTNAG|C` is `exact` with both windows resolved). The
  entry therefore also carries `gap_unresolved_bases`, the gap's own
  bases outside A/C/G/T wherever they sit, and the record's
  `gaps_with_unresolved_bases` is the count that certifies the gap
  interior; on this table it is 0 for every gap, as the windows counter
  already was for the 1- and 2-base gaps, whose bases the windows cover.
  `length_floors.py --short-gaps` runs the same test over a chromosome,
  fetching two flanking bases per gap (one sequence request each, run
  2026-09-09, re-run with the class columns the same day: same thirteen
  gaps, same flanks):

  | Assembly | Chrom | Strand | Gap (0-based, half-open) | Length | In CDS | Transcript | Flank / gap / flank | Donor window | Acceptor window | Class |
  |---|---|---|---|---|---|---|---|---|---|---|
  | sacCer3 | chrIV | - | 518062-518063 | 1 | yes | NM_001184379.4 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | - | 649820-649821 | 1 | yes | NM_001184417.2 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | - | 882621-882622 | 1 | yes | NM_001184436.2 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | - | 991043-991044 | 1 | yes | NM_001184421.2 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | + | 873404-873405 | 1 | yes | NM_001184419.4 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | + | 982754-982755 | 1 | yes | NM_001184423.4 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | + | 1097368-1097369 | 1 | yes | NM_001184425.2 | `TT A GG` | no | no | neither |
  | sacCer3 | chrIV | + | 1208301-1208302 | 1 | yes | NM_001184427.2 | `TT A GG` | no | no | neither |
  | ce11 | chrIII | - | 13377845-13377860 | 15 | no | NM_001268269.2 | `TA TCTTTGAATAAAAAC AA` | no | no | neither |
  | ce11 | chrIII | + | 1879047-1879062 | 15 | no | NM_001027675.5 | `CG ATGATTTTCTCAAAA AT` | no | no | neither |
  | ce11 | chrIII | + | 4530490-4530491 | 1 | yes | NM_001330837.3 | `AA A CT` | no | no | neither |
  | ce11 | chrIII | + | 4530511-4530512 | 1 | yes | NM_001330837.3 | `AG G AA` | no | no | neither |
  | ce11 | chrIII | + | 13344922-13344940 | 18 | no | NM_001383022.1 | `AC GTTTTTATTTACAGAACC AC` | yes | no | donor only |

  All eight yeast gaps are the same seven bases, `CTT A GGC` read across
  the gap: the Ty1 +1 programmed frameshift site, at which the ribosome
  slips from the CTT codon to the AGG (Belcourt and Farabaugh 1990,
  doi:10.1016/0092-8674(90)90371-K), so the "intron" is the skipped A.
  The two worm gaps lie 21 bases apart in one `cdh-4` transcript, and the
  three 15- to 18-base gaps are on the UTR side of a CDS. None of the
  thirteen passes the overlapping-window test, and the class column says
  how each fails: twelve fail both windows, and the 18-base UTR-side gap
  in NM_001383022.1 reads GT at its donor end but CC where AG would be,
  so it is `donor_only`, the one gap here a donor-only mask would let
  through. No window was unresolved. So a motif-masked decoder could not
  emit any of them as an intron even without a duration floor, and none
  is a splice site under the scorer's rule or in these labels; a predictor that joins the two CDS blocks
  through the gap is scored as right about the protein and, since the
  scorer drops the gap from both sides, not charged for the "intron".
  The ce11 count here is five, not the twelve of the table above, because
  that column counts introns under 30 and this one gaps under 20.
- **Informants**: one row per informant, in an order fixed per source so
  that examples stack into batches. For a UCSC track (multiz, Cactus) that
  is the track tree's depth-first leaf order minus the reference, and `K`
  is fixed per track. Ensembl EPO has no track tree, only one tree per
  block, so there the rows are the species set's members in sorted order
  (from `info/compara/species_sets/<method>`, which the fetcher now records
  in the manifest; a member that does not align in the window is an
  all-unaligned row), then the inferred ancestral rows sorted by name. The
  leaf part is fixed per species set; the ancestor part varies by window
  (`Ggal-Mgal[2]` exists only where chicken and turkey both align), so a
  model that wants a fixed `K` on EPO takes the leaves and `--drop-ancestors`,
  and one that wants ancestral states reads the sidecar's `ancestral` list
  and `ancestral_clades`. Each entry is a base, a **gap** (the informant is aligned here with a deletion)
  or **unaligned** (no block covers the position); the distinction is the
  one CONTRAST and N-SCAN made and section 6 shows why: past about one
  substitution per site most non-coding positions are unaligned, not
  gapped. Insertions relative to the reference are not expanded; their
  length after each reference base is kept in a separate channel (clipped
  at 255), so the example keeps a fixed length whatever the alignment did.
- **Block selection**: a MAF can say two things about one reference base
  and one species, and the cutter's choice is explicit and counted in the
  sidecar's `block_selection` rather than left to file order (engels'
  audit of N-SCAN's `maf_to_align.pl`, relay note
  20260909T182611Z-engels-0019, found that converter silently skipping a
  whole block on overlap and misplacing a minus-strand reference row).
  Three cases. (1) *Blocks overlapping on the reference*: the first block
  in file order to cover a position wins; the sidecar counts the
  positions covered twice, the aligned informant bases the losing blocks
  carried there, and how many of those the winner did not have
  (`overlap_informant_bases_lost`, the unrecoverable part, per
  informant). (2) *A reference row on the minus strand*: the block is
  reverse-complemented into forward coordinates (`srcSize - start -
  size`, the UCSC MAF rule) and painted, not skipped;
  `reference_minus_strand_blocks_flipped` counts them. (3) *Several rows
  for one species in one block*: Cactus exports every copy of a
  duplicated region, and the old cutter let the last copy overwrite the
  earlier ones column by column. `--duplicate-rows identity` (default)
  keeps the copy with the most bases identical to the reference in that
  block, ties to the first; `first` keeps the first row. The discarded
  copies are counted in three units that do not agree with each other,
  and the sidecar names the unit of every counter in
  `block_selection.units` (lenin's point in relay note
  20260909T191730Z-lenin-0018: a per-informant table must say whether
  it counts rows, copies or bases): *rows* are MAF rows discarded
  (`duplicate_rows_discarded`, and `rows` per informant); *copies* are
  how many rows one informant had in one block (`max_copies` per
  informant, `duplicate_max_copies` overall) and in how many blocks it
  had more than one (`blocks` per informant, `duplicate_row_blocks`
  overall, `duplicate_informants` for the number of informants
  affected); *bases* are aligned informant bases inside the window that
  the discarded rows carried (`duplicate_rows_discarded_bases_in_window`,
  `bases` per informant), which is the one unit comparable with the
  example itself, so the sidecar also reports
  `kept_informant_bases_in_window`, the cells of the example that hold a
  base. Measured on 2026-09-09 on an alignment the inventory has since
  excluded (section 5; the numbers stay as a recorded finding): the human
  HBB window on the
  Cactus 241-way (4,931 bp, 2,517 blocks, 240 informants) has **no**
  overlapping blocks and no minus-strand reference rows, but 2,513 of
  its 2,517 blocks carry a duplicated species, up to 17 copies of one
  species in one block, 322,650 extra rows in all, holding 529,613
  aligned bases inside the window against 863,014 base cells in the
  example that keeps one copy each, so the discarded copies hold 61% of
  what the example keeps. (The first version of this paragraph compared
  against 899,732, a count of every cell that was aligned in any state,
  gap columns included: 900,643 on the re-fetch of 2026-09-09; that is
  not the unit the discarded-bases counter uses, and the sidecar now
  reports both sides in one unit.) 189 of the 240 informants have a
  discarded copy; the 51 that do not are the primates, the caviomorph
  rodents and a few others whose beta-globin cluster is single-copy
  against human. The three units rank the affected informants
  differently, which is why each is stated: by discarded *rows* and by
  discarded *bases* the top ten are the same ten (giraffe 9,945 rows and
  16,772 bases, Sowerby's beaked whale 8,441 and 14,678, goat 7,353 and
  12,434; the median discarded row carries 1.6 bases inside the window
  because Cactus blocks here are about 2 bp long, and one primate row
  carries 15), but by *blocks affected* only four of those ten remain
  (giraffe again at 2,207 of 2,517 blocks, then okapi 2,197, white-tailed
  deer 2,196 and Père David's deer 2,194, the cervids and giraffids
  having a few copies in nearly every block), and by *copies in one
  block* the order is different again (Asian palm civet 17, Gambian
  pouched rat 16, Sowerby's beaked whale 15, Chinese hamster 14, so
  rodents and a carnivoran that were far down the rows list). A reader
  who takes "189 informants affected" as one number will be wrong by
  more than a rounding whichever unit they had in mind: per informant
  the blocks affected run from 1 to 2,207 (median 676) and the rows per
  affected block from 1 to 4.9. The two policies differ on 1.5% of the example's cells
  (18,329 of 1,183,440); `identity` leaves 655,643 cells identical to
  the reference against 638,479 for `first`, so `first` picks a paralog
  where a closer copy exists in about one column in sixty. The fly Adh
  window on the multiz 124-way has none of the three cases, as a
  single-coverage chain/net alignment should. Which copy is the
  *ortholog* is a different question from which is most identical, and
  neither policy answers it; see section 9 item 11. The fetcher's
  coverage tables (sections 6 and 6.1) count a species as aligned where
  *any* copy has a base, so on Cactus they read as an upper bound on
  orthologous alignability.
- **Geometry**: the patristic distance from the reference to each informant
  on the track tree, and the Newick string in the sidecar. Anything else
  the model wants (an MDS embedding of the tree as in HyphAeon, a
  per-branch rate) is derived downstream from the tree, not stored.
- **Conservation**: the phyloP or phastCons value per base when the track
  had one, NaN otherwise, so a model can be trained with and without it.
- **Leakage**: `--drop-species` removes held-out informants before the
  tensors are built, and with them every EPO ancestral row whose clade,
  read from the block trees' internal node names, contains a dropped
  species (a reconstruction that used the held-out genome is as much a leak
  as the genome itself; the sidecar lists them as `dropped_ancestors`). For
  a reference-anchored alignment that is exact. For a jointly inferred one
  the row goes but its influence on the alignment stays, which the sidecar
  records as `dropped_rows_only: true`; section 7 says why the honest tool
  for those is a declared informant set at inference. The cutter decides
  the class from an explicit table of track names (multiz `*way` and
  pairwise chain/net tracks are reference-anchored; Cactus, `hprc*way`,
  and every Ensembl multiple alignment are jointly inferred) and treats a
  track it does not know as jointly inferred, which is the strict branch of
  `docs/benchmark.md` section 3.2; `--reference-anchored` overrides that
  for a track the operator has checked, and the sidecar records the
  override as the operator's claim rather than the table's. The sidecar's
  `dropped_species` list is what the run copies into the
  `alignment_rows_dropped` key of its benchmark declaration (section 3.3);
  a sidecar with `dropped_rows_only: true` means the declaration must also
  say the alignment is jointly inferred, and the rebuild rule applies.
- **Orientation**: reverse-strand genes are not flipped. `--both-strands`
  writes the reverse-complement example (labels and boundaries move with
  the coordinates, informant bases complemented, insertion lengths shifted
  to the preceding base), so a model can be trained orientation-symmetric;
  otherwise augmentation is the loader's job.

Sizes: at `L` = 4,096 and `K` = 122 (fly 124-way) an example is 1.0 MB
uncompressed and 70 KB compressed (measured on `Idgf3`); the human
470-way at `K` = 469 is 3.9 MB uncompressed. The whole fly genome at stride 2,048 is about 70,000
examples, the human genome 1.5 million; neither needs to be materialised,
the cutter can run inside a data loader on the fetched windows.

Open choices, deliberately left to T-human-011: whether the informant axis
should be the tree leaves (fixed `K`, mostly unaligned rows at distance)
or the aligned subset per window (variable `K`, denser); whether to encode
codon position on the reference axis as an input (it is a label here);
whether UTR should be split into 5' and 3' (RefSeq lets us, the
coverage tables do not yet); and whether EPO ancestral rows are an input
(the only source that has them, section 3.2) or noise to drop.

### 6.4 Pairwise codon export for the KA/KS baseline

T-human-010 will fit the KA/KS test of Nekrutenko, Makova and Li (2002) to
reference/informant codon pairs cut from these windows. stalin's audit of
the current PAML source (relay note 20260909T164106Z-stalin-0016,
`abacus-gene/paml` at `4c7902f`) found that codeml's pairwise mode forces
`cleandata = 1` and that its reader deletes every codon column in which
*any* input row has a gap or an ambiguity, before a pair is fitted. Given a
multispecies block, codeml therefore fits a pair on the columns that are
clean across the whole block, not the pair. `scripts/data/pairwise_codons.py`
materialises each pair on its own from the reference CDS (transcription
order, reverse complemented on the minus strand, codons may span splice
junctions, terminal and internal reference stops removed and reported)
and one informant row, writes a two-row sequential PHYLIP per informant
with the codons retained by the *pairwise* rule (all three bases A, C, G
or T; no informant insertion between consecutive codon bases; no informant
stop), and records per informant the codon counts by fate (`unaligned`,
`gap`, `ambiguous`, `insertion`, `informant_stop`), identity, tree distance
and the count that the *complete-case* rule over all rows would have kept.
`--complete-case` writes the complete-case columns instead, so the two
policies can be fitted on identical inputs; `--write-multi` writes the
multi-row block that would reproduce codeml's own deletion. Self-test: 17
checks built on stalin's three-row example, seed-independent.

Measured on the fly Adh window (dm6 124-way, `apiMel4` dropped, 122
informant rows, `longest-cds` picks NM_001032098.2, 256 codons; run
2026-09-09, byte-identical under `PYTHONHASHSEED` 0 and 42):

| Informant set | Rows | Complete-case columns | Pairwise codons kept (min / median / max) | Codons lost to complete case, summed over pairs |
|---|---|---|---|---|
| all rows | 122 | 0 of 256 | 0 / 119 / 256 | all |
| within 2 subst/site of dm6 | 38 | 0 of 256 | 0 / 253 / 256 | 8,545 of 8,545 |
| within 1 subst/site of dm6 | 13 | 226 of 256 | 226 / 256 / 256 | 360 of 3,298 (11%) |

47 of the 122 rows are unaligned at every Adh codon, so one absent
informant is enough to delete the whole gene under complete-case
deletion; even among the 13 drosophilids inside the non-coding
alignability horizon of 6.1, the 30 codons `droRho2` lacks are removed
from the other twelve pairs. Kept codons fall with distance the way the
CDS coverage of 6.1 does: mean 1.00 of the CDS below 0.5, 0.82 between 1
and 2, 0.26 beyond 2 (Adh), so the baseline's coverage accounting (which
pairs exist and on how many codons) has to be reported next to its
classification, as stalin's note recommends. The export is what T-human-010
should feed codeml, one pair per file; a shared complete-case filter is a
policy that has to be named and costed, not a default.

## 7. Suitability for training versus held-out evaluation

The benchmark's rules (`docs/benchmark.md` section 3.2): held-out labels may
never be seen; held-out sequence in a pretraining corpus must be declared;
alignments used for *training* must not carry a held-out species, and
section 3.2 splits that rule by how the alignment was built: a jointly
inferred alignment (Cactus, Ensembl EPO) must be rebuilt without the
held-out species, while a reference-anchored one (multiz) may have the row
dropped at cut time provided the run lists the dropped rows in the
`alignment_rows_dropped` key of its declaration (section 3.3); an informant
set used at *inference* on a held-out target is permitted and must be
declared. `cut_windows.py` applies the split from a table of track names
and treats any track not in it as jointly inferred, so a wrong guess fails
strict and visible rather than permissive and silent (section 6.3).

Applying them to what exists:

| Panel species | Split | Public multiple alignment referenced on it | Conservation | Fit |
|---|---|---|---|---|
| human | heldout_paired | hg38 100/30/470-way, Ensembl mammals/primates/amniotes (Cactus 241/447 excluded for now, section 5) | phyloP, phastCons | evaluation only; any of these is a declarable informant set at inference. Not usable to train, because training on human windows means training on human labels |
| mouse | train | mm39 35-way; Ensembl 44/92-mammal EPO and 22-murinae re-referenced on mouse | mm39 phyloP35way | training source, but every one of these alignments contains human (held-out) sequence, and the 35-way also carries chicken (`galGal6`, held-out). For the 35-way (multiz) drop the `hg38` and `galGal6` rows at cut time and declare `alignment_rows_dropped: [hg38, galGal6]` (`scripts/data/informant_membership.tsv`); the EPO sets are jointly inferred, so dropping the row is not enough and they would have to be rebuilt |
| chicken | heldout | Ensembl sauropsids EPO on GRCg7b (exact panel assembly); galGal6 77-way (older assembly, liftover needed) | galGal6 phyloP77way only | evaluation with a declared informant set; the sauropsid set contains no other panel species |
| zebrafish | train | Ensembl fish EPO on GRCz11 (panel is GRCz12ab; `docs/benchmark.md` section 2.4 keeps it there) | none on either assembly (danRer11 has no `*way`, phyloP or phastCons track; GRCz12ab has no hub yet) | no usable training alignment: both fish EPO sets (32 and 65 species) contain fugu (held-out) and are jointly inferred, so under section 3.2 they must be rebuilt without fugu whichever assembly zebrafish sits on; dropping the row does not qualify. Lifting to GRCz11 through UCSC's `danRer11ToGCA_052040795.1.over.chain.gz` buys only that same unusable alignment and was rejected in section 2.4. Counts with the thirteen in section 8 |
| fugu | heldout_paired | Ensembl fish EPO on fTakRub1.2 (panel 1.3) | none | evaluation; informant set contains zebrafish (train), which the rules allow at inference if declared |
| fruit fly | train | dm6 124-way and 27-way | dm6 phyloP124way | training; the 124-way contains honey bee (`apiMel4`, held-out): drop that row at cut time and declare `alignment_rows_dropped: [apiMel4]` |
| C. elegans | train | ce11 135-way (raw MAF only, readable by Range from `/gbdb`) | ce11 phyloP135way | training; the tree carries Ciona (`ci3`, held-out) as a distant outgroup, so drop it at cut time and declare `alignment_rows_dropped: [ci3]`; no other panel species among informants, but no informant within 1 substitution per site either (section 6.1), so comparative evidence for non-coding classes is thin |
| yeast | train | sacCer3 7-way | phastCons7way | training; near-intronless, so mostly a negative control for the intron machinery |
| frog, honey bee, sea anemone, ciona, thale cress, rice, maize, S. pombe, Neurospora, Dictyostelium, Plasmodium, Tetrahymena | 5 train, 7 held-out | none (rice: 8-way EPO on IRGSP, not the panel's AGIS1.0) | none | no comparative input available; the model runs single-genome on these unless we build alignments |

Two consequences for the design:

- **Every clade-general claim will be tested without informants** on seven
  of the ten held-out species, whatever the model, and on the training side
  zebrafish joins the twelve alignment-less species until a fish alignment
  without fugu is built (thirteen of twenty for training windows). The comparative branch of the
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
  set at inference is the right tool, and it is allowed. `docs/benchmark.md`
  section 3.2 now states exactly this split, and section 3.3 makes the
  `alignment_rows_dropped` declaration mandatory: audited leaf by leaf
  against all ten held-out species (`scripts/data/informant_audit.py`,
  2026-09-18), it is `[apiMel4]` for dm6 124-way and 27-way, `[hg38,
  galGal6]` for mm39 35-way, `[ci3]` for ce11 135-way and empty for
  sacCer3 7-way, and `benchmark/score.py` refuses a run that omits the key.

## 8. What we would have to build, with costs

For the thirteen species without a usable public alignment (the twelve with
none, plus zebrafish, whose only alignment is jointly inferred with held-out
fugu), the options are:

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
   Zebrafish is the one training species in this list, and the natural
   job for it: a fish clade on GRCz12ab (contig N50 59 Mb, no gaps;
   `docs/benchmark.md` section 2.4) with fugu left out, so the result is
   compliant with section 3.2 by construction rather than by dropping a
   row. Its genome is 1.4 Gb, so the cost sits between the two cases
   above.

None of this is needed for the Phase 1 to 3 deliverables. It is a Phase 4
cost, and section 7's first consequence (the model must run without
informants) is the design-level answer.

## 9. Open items and questions for other owners

1. Done: section 6.1 aggregates coverage over 12 fly, 10 human, 10 mouse
   and 10 worm windows by tree distance (`scripts/data/coverage_by_distance.py`;
   the mouse and worm draws are reproducible with
   `scripts/data/sample_genes.py`). The mm39 35-way and ce11 135-way read
   through the same wigMaf index plus Range path as the 100-way and
   124-way (the ce11 MAF sits at `/gbdb/ce11/multiz135way/chr*.maf`, without
   a `maf/` subdirectory; the fetcher tries both layouts). Still one
   chromosome per species and short genes for human. The mouse run found
   and fixed a parsing bug: GenArk sources (`GCF_003668045.3.NC_048596.1`)
   were cut at the first dot and never matched their tree leaf
   (`GCF_003668045v3`), so that informant scored as fully unaligned;
   `maf_source()` in the fetcher and the cutter now handles the convention.
2. Done this draft: the window-cutting convention is section 6.3 and
   `scripts/data/cut_windows.py`; the three open design choices at the end
   of 6.3 are for T-human-011.
3. Answered: lenin (T-human-007, note 20260909T082157Z-lenin-0010) left
   the "with informants / without informants" column of `benchmark/report.py`
   to T-human-011, on the grounds that it is a question about what the
   model is allowed to see rather than about the scorer. Section 6.3's open
   choices and section 7's first consequence are where T-human-011 should
   pick it up. The same note adopted the reference-anchored versus jointly
   inferred split proposed here (section 7) into `docs/benchmark.md`
   section 3.2 and added the `alignment_rows_dropped` declaration key.
4. Question for whoever takes T-human-009: a measured lastz wall-clock and
   memory for one genome pair at 100 Mb and at 2 Gb, to cost section 8
   option 2.
5. Done this draft: the Ensembl licence text is quoted in section 3. The
   470-way species composition is in `hg38.470way.scientificNames.nh` next
   to the bigMaf (listing read 2026-09-09); the per-order count of its
   leaves with distances from human is now the composition table in
   section 2.1 (`scripts/data/tree_composition.py`,
   `scripts/data/hg38.470way.orders.tsv`).
6. Answered: lenin (T-human-007, note 20260909T101258Z-lenin-0012,
   `docs/benchmark.md` section 2.4 on PR #5) keeps zebrafish on GRCz12ab.
   Moving to GRCz11 would buy only the Ensembl fish EPO, and both fish
   species sets contain the held-out fugu, so under section 3.2 that
   alignment is rebuild-only on either assembly; GRCz12ab is the better
   assembly to rebuild on (25 contigs, no gap bases, versus 19,725 contigs
   and 4.69 Mb of N in GRCz11, which matters for the panel's third-longest
   introns). Sections 1, 7 and 8 now count zebrafish with the species that
   have no usable alignment: thirteen of twenty for training windows.
7. New: `hgdownload.soe.ucsc.edu` reset every connection for several
   minutes during the sampling run while `hgdownload2` served the same
   files. The fetcher now rotates mirrors, but a training pipeline that
   reads MAF by HTTP Range at scale should rsync the per-chromosome MAF
   once (sizes in section 2.1) rather than depend on hgdownload staying up.
   That holds for the 100-way (790 GB raw, 71 GB gzipped) and the 35-way;
   it does not hold for the 470-way, whose raw MAF is 5.9 TB with no
   gzipped copy, so for that track the bigMaf API and Range reads into
   `maf/chr*.maf` are the only laptop-scale paths and a full copy is a
   cluster storage request.
8. Done after review 20260909T090919Z-lenin-0011: the cutter's alignment
   class now comes from an explicit table with a strict default and an
   operator override; EPO ancestral rows are read from the fetcher's
   manifest and dropped with their clade under `--drop-species`; the
   sidecar's `label_counts` is counted through the mask so the reverse
   complement of a padded window reports the same counts; EPO row order is
   the species set rather than a tree, and 6.3 says which part of `K` is
   fixed; a transcript-type filter matching the benchmark's truth rule is
   the default; and the sidecar counts blocks skipped for a minus-strand
   reference row. The same pass fixed a bug lenin did not have to find:
   Ensembl block-tree leaves are named `species_region_start_end[strand]`,
   which the distance lookup never matched, so `dist` was NaN for every EPO
   informant; leaves are now mapped back to row names, and ancestors get
   their distance from the internal node.
9. New: `scripts/data/pairwise_codons.py` (section 6.4) exports
   reference/informant codon pairs for T-human-010, one pair per file with
   original and retained codon counts and drop reasons, after stalin's
   note 20260909T164106Z-stalin-0016 showed that codeml's pairwise mode
   deletes columns over every row it is given. On the fly Adh window,
   complete-case deletion over the 122 informant rows keeps no codon at
   all and over the 13 drosophilids within 1 substitution per site costs
   11% of the pairwise codons. For T-human-010: fit each pair from its
   own file and report coverage (pairs that exist, codons retained) next
   to the classification.
10. Answered: lenin (T-human-007, note 20260909T171958Z-lenin-0016) noted
    that 5.3% of GENCODE 50's CDS-bearing transcripts end incompletely and
    say so only by omission. Section 6.3 now records what each source the
    fetcher reads declares (RefSeq: `cmpl` at both ends for every
    transcript checked; UCSC's GENCODE track: stat columns `none`, NF tags
    and exon frames present; Ensembl REST: nothing), the fetcher carries
    the declaration, and `cut_windows.py --partial-ends` takes it where it
    exists and reads the reference sequence where it does not, a rule that
    agrees with GENCODE's tags on 99.7 to 99.9% of hg38 chr21 transcripts.
    Incomplete ends get no start or stop mark and are left out of
    `distinct_sites`, so the accounting lenin asked about is per example
    in the sidecar's `cds_ends`.
11. New, from engels' converter audit (relay note
    20260909T182611Z-engels-0019): section 6.3 now states and counts the
    cutter's block-selection rules. The finding that matters is not the
    two hazards engels tested (the HBB Cactus window has no overlapping
    blocks and no minus-strand reference rows) but a third: nearly every
    Cactus block carries several rows for one species, and those extra
    copies hold 61% as many aligned bases as the example keeps, in the
    one unit both sides now share (the earlier 59% divided by a count
    that included gap cells). Open for
    T-human-011: whether the model should see the copies (a paralog
    channel, or several rows per species) or one chosen copy, and if one,
    chosen how. `identity` is a stand-in; the defensible rule is synteny
    (keep the copy that continues the previous block's contig and
    coordinates), which needs the informant coordinates the cutter now
    ignores, or the UCSC `--noDupes` style export if the track offers
    one. Until then, Cactus examples carry `duplicate_row_policy` and
    the per-informant counts so a result can say which copies it saw.
12. From engels' Helixer audit (relay note 20260909T202745Z-engels-0020)
    and stalin's converter audit (relay note 20260909T204115Z-stalin-0020):
    the sidecar's `feature_lengths` and `scripts/data/length_floors.py`
    count what lies under a decoder's minimum durations (section 6.3).
    The T-human-007 half is settled: lenin answered (note
    20260909T211910Z-lenin-0019) that the scorer has counted CDS gaps
    under 20 bases apart from splice sites since its run 2, and the
    cutter follows that floor as `--min-intron` from version 0.9
    (section 6.3, "the intron floor"), so a frameshift gap is neither a
    scored junction nor a labelled one. Still open, for whoever reads a
    GFF3 source into the fetcher: carry NCBI's `exception=ribosomal
    slippage` and Ensembl's equivalent into the transcript record so a
    short gap's `exception` field is filled from provenance rather than
    inferred from length (stalin, note 20260909T214050Z-stalin-0021), and
    keep an eye on *Stentor*-like clades where the floor must be lowered
    to 15. For T-human-011, two hard ceilings now have numbers: a quarter
    of worm introns and half a percent of fly introns lie under the
    50-base floor a U2 duration model imposes, and a fixed motif alphabet
    excludes 1.1 to 1.6% of human and fugu reference introns (lenin's
    fugu strata, same note, where Tiberius emits no AT-AC or non-canonical
    intron at all); whatever replaces the HMM decoder should have neither
    a minimum intron duration above the clade's shortest real introns nor
    a fixed motif set, and the counts above are the price of each choice
    on one chromosome per clade.
