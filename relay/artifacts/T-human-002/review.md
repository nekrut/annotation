# Independent review of eukaryotic gene prediction — `lenin`

Task: T-human-002 (slot 1 of 5 — a fifth slot, T-human-012, was added by
`human` decision 20260909T015500Z). Status: **draft, run 2 of an expected 3.**
Author: `lenin` (Claude Code / Opus 5). Last updated: 2026-09-09 (run 2).

This draft was written blind. I have not opened any other agent's review
artifact and will not until this task is in `review`, per charter
§Independence. One caveat, recorded rather than hidden: `trotsky` posted a
broadcast `note` (20260909T014817Z) summarising its own T-human-012
architectural conclusions, which my inbox delivered before I could avoid it.
I did not open its artifact, and nothing in §6 below was changed in response
— §6 was written and pushed in run 1, before that message existed. The
synthesis (T-human-006) should treat my §6 as pre-dating it.

Everything below that carries a number carries a citation or names the
artifact that produced it (`repos.tsv`, or the search log in §1). Where I
have a claim I could not verify in this run, it is marked **[unverified]**.
Nothing paywalled was retrieved; only abstracts and metadata were pulled for
non-OA items.

---

## 1. Search log

All searches run 2026-09-09 (UTC), scripted, no browser. Runs 1, 2 and 3
were the same day; rows 5–8 are run 2, rows 9–13 are run 3.

| # | Source | Interface | Queries | Notes |
|---|--------|-----------|---------|-------|
| 1 | Europe PMC | REST `/search`, `resultType=core`, sort by citations | 20 broad topical queries (see below) | Ranking by citation count pulled in irrelevant high-citation papers for short/ambiguous names (`CONTRAST`, `Exonerate`, `Liftoff`, `TWINSCAN`); those needed targeted re-queries. |
| 2 | Europe PMC | REST `/search`, `TITLE:` field queries | 25 targeted queries on named tools | Recovered GENSCAN, SNAP, GlimmerHMM, GeneID, miniprot, Splign, TOGA, CAT, SpliceAI, minisplice, Evo 2, AlphaGenome, phyloP, Cactus, OMArk, ANNEVO. |
| 3 | Europe PMC | REST `/search` by DOI | 19 DOIs | Pulled verified abstract, journal, PMCID, OA flag and licence for each headline method. |
| 4 | GitHub | REST API v3 (`/repos/{full}`, `/commits?since=`) | 32 repositories | Last push, commits in trailing 12 months, open issues, stars, language, SPDX licence, archived flag. Output: `repos.tsv`. Rate limit reached on the last 2. |
| 5 | GitHub | REST API v3 `/search/repositories` | 8 name searches | Resolved the repositories run 1 could not find: MAKER, GlimmerHMM, ensembl-anno, ANNEVO, minisplice, funannotate, OpenSpliceAI, GeneMark. Separate quota from `/repos`, which was still rate-limited. |
| 6 | Crossref | REST `/works/{doi}` | 4 DOIs | Verified the two DOIs run 1 marked `TO VERIFY` (TWINSCAN 2003, CONTRAST 2007) and completed the ANNEVO and GALBA author lists. |
| 7 | raw.githubusercontent.com | README/LICENSE fetches | egapx, BRAKER, Helixer, GeneMark-ETP | Primary-source verification of hardware floors, clade exclusions and licence terms — see §3 and §5. |
| 8 | local shell | actual installs, fresh venvs, Python 3.13.9, gcc 15.2, RTX 5080 | Tiberius, egapx, miniprot, minisplice, Helixer (probe) | Filled the `install_tested` column. Results in §3.1. |
| 9 | OpenAlex | REST `/works?search=`, `from_publication_date:2022-01-01` | 8 topical queries, 25 hits each | 163 unique works. Relevance ranking is noisy (it returns high-citation genome papers that merely *used* a predictor), but it independently re-surfaced everything already in §2 and added the Helixer *Nature Methods* version and [djossou2025overview]. Useful mainly as a **negative** control on run-1 coverage. |
| 10 | Europe PMC | REST `/search` with `SRC:PPR` and date filters | 6 preprint queries | **The productive pass of this run.** 281 + 65 + 131 + 57 hits; yielded the entire 2025–2026 deep-learning annotation cohort that runs 1 and 2 missed (§2.1). |
| 11 | Europe PMC | REST `/search` by DOI, `resultType=core` | 17 DOIs | Verified title, authors, licence and abstract for every row added in §2.1, plus the ANNEVO preprint. |
| 12 | GitHub | REST API v3 via authenticated `gh api` | 13 `/repos` + 7 `/search/repositories` | Backfilled the five run-2 rows whose counts were `?` (unauthenticated quota had been exhausted) and resolved the run-3 repositories. |
| 13 | raw.githubusercontent.com + PyPI | README and package metadata | Vipsania, OrionGeno, GeneCAD, TOGA2 | Primary-source licence and hardware terms for the new cohort — see §3.2. |
| — | arXiv | REST `/api/query` | 5 queries, 3 retries each, custom User-Agent | **Failed: HTTP 429 on every attempt from this host, including a single query in isolation.** Covered indirectly through OpenAlex `type:preprint`, which indexes arXiv; recorded as a gap rather than claimed as done. |

Broad queries in pass 1 (verbatim): GENSCAN; AUGUSTUS ab initio eukaryotic;
BRAKER pipeline; GeneMark-ES self-training; SNAP Korf; Gnomon NCBI; EGAPx;
MAKER; TWINSCAN/N-SCAN; CONTRAST CRF; Helixer; Tiberius; DNA language model
annotation; SpliceAI; Ka/Ks comparative human-mouse; Zoonomia; TOGA; GALBA;
StringTie; benchmark comparison of gene predictors.

Hit counts for the useful ones: `TITLE:"TWINSCAN"` 4, `TITLE:"N-SCAN"` 2,
`TITLE:"Gnomon"` **0**, `TITLE:"EGAPx"` **0**, `TITLE:"ANNEVO"` 2,
`TITLE:"miniprot"` 3, `Helixer` 44, `Tiberius` 6, `EGAPx` (free text) 41 —
all 41 being genome papers that *used* it, not method papers.

**Finding from the search log itself, worth stating plainly:** the two tools
the charter singles out as production reality — NCBI's **Gnomon** and
**EGAPx** — have no methods paper. `TITLE:"Gnomon"` and `TITLE:"EGAPx"` each
return zero hits in Europe PMC. They are describable only through NCBI
database papers [kuhn2013refseq; sayers2026ncbi] and through the repository
documentation. Any benchmark that includes them is benchmarking software, not
a published method, and that asymmetry should be recorded in T-human-007.

### 1.1 What is still not covered (updated after run 2)

Listed so the next run and the synthesis (T-human-006) know the shape of the
hole, rather than assuming coverage:

- ~~No bioRxiv, arXiv, OpenAlex or Semantic Scholar pass yet.~~ **Closed in
  run 3 for bioRxiv (via Europe PMC `SRC:PPR`) and OpenAlex.** This was the
  most consequential gap in the whole review: it produced §2.1, eleven
  methods from 2025–2026 that runs 1 and 2 did not see at all, including two
  that bear directly on the charter's central hypothesis. **arXiv is still
  not covered directly** — its API returns 429 to this host — and Semantic
  Scholar was not attempted. OpenAlex `type:preprint` partly compensates.
- ~~No GitHub *code search* pass.~~ **Partly closed in run 3**: seven name
  searches resolved the run-3 cohort's repositories, and the five run-2 rows
  with `?` counts are now filled. A true dependents/code snowball
  (`/search/code`, `/network/dependents`) was still not run.
- ~~No install was attempted for any repository.~~ **Closed in run 2 for the
  four tools that matter most** — see §3.1. Two build clean, two install
  partially, one (BRAKER) was not attempted and the reason is documented.
  37 of 41 rows read `install_tested=no` after run 3 (2 `yes`, 2 `partial`); most of those are context
  repositories, not baselines.
- Runtime and memory figures below are as-reported by authors. Nothing was
  measured here. Measurement belongs to T-human-009, not to this task, but
  the reported figures are not comparable across papers (different hardware,
  different genomes, different definitions of "annotation") and I flag that
  rather than tabulate them as if they were.

---

## 2. Publications table

Approach classes: **GHMM** = generalized HMM ab initio; **COMP** =
comparative/alignment-informed; **EVID** = evidence-based pipeline;
**DL** = deep learning; **HYB** = hybrid. Accuracy figures are as reported by
the cited authors on their own benchmark; they are *not* mutually comparable.

| Method | Year | DOI | Class | Inputs required | Clades trained/evaluated | Reported accuracy (benchmark, metric) | Runtime / hardware | Cross-species evidence | Code |
|---|---|---|---|---|---|---|---|---|---|
| GENSCAN | 1997 | 10.1006/jmbi.1997.0951 | GHMM | genomic DNA only | human/vertebrate | 75–80% of exons identified exactly, on standardized human/vertebrate sets [burge1997genscan] | not reported | Author reports consistent accuracy across C+G content and vertebrate groups [burge1997genscan] | web only |
| GeneID | 2000 | 10.1101/gr.10.4.511 | GHMM | DNA | Drosophila, human | see paper | very fast, C | per-species parameter files | guigolab/geneid |
| GlimmerHMM / TigrScan | 2004 | 10.1093/bioinformatics/bth315 | GHMM | DNA | plant, human, misc. | see paper | fast | per-species training | see `repos.tsv` note |
| SNAP | 2004 | 10.1186/1471-2105-5-59 | GHMM | DNA + training set | multiple novel genomes | conclusion: "every genome needs a dedicated gene finder" [korf2004snap] | fast | **Explicitly negative**: foreign gene finders are "highly inaccurate"; the nearest phylogenetic neighbour is not necessarily the best donor [korf2004snap] | KorfLab/SNAP |
| GeneMark-ES | 2005 | 10.1093/nar/gki937 | GHMM (self-training) | DNA only | fungi + novel euk. | comparable to or better than supervised training [lomsadze2005genemarkes] | not reported | self-training removes the per-species labelled-data requirement, not the per-species *fit* | non-free (see §3) |
| TWINSCAN / N-SCAN | 2003 / 2007 | 10.1101/gr.830003 / 10.1002/0471250953.bi0408s20 | COMP+GHMM | target DNA + informant genome(s) | human/mouse, rat, worm | rat and C. elegans predictions validated by RT-PCR and sequencing [10.1101/gr.1959604; 10.1101/gr.3329005] | not reported | The historical proof that one informant genome buys real accuracy | unmaintained |
| CONTRAST | 2007 | 10.1186/gb-2007-8-12-r269 | COMP+DL(discriminative) | target + multiple informants | vertebrate | see paper | not reported | phylogeny-*free* multiple-informant design — directly relevant to our geometry question | unmaintained |
| MAKER / MAKER2 | 2008 / 2011 | 10.1101/gr.6743907 / 10.1186/1471-2105-12-491 | EVID | DNA + EST/protein/RNA-seq + ab initio predictors | any, community-driven | outperformed by BRAKER3 [gabriel2024braker3] | heavy, MPI cluster | pipeline, not a model | canonical repo **not found** (see `repos.tsv`) |
| EVidenceModeler | 2008 | 10.1186/gb-2008-9-1-r7 | EVID (combiner) | multiple predictions + evidence | any | see paper | light | combiner | active |
| AUGUSTUS (+hints) | 2006 | 10.1186/1471-2105-7-62 | GHMM (+EVID via hints) | DNA, optional hints | many species, per-species parameter sets | the long-standing reference ab initio; superseded ab initio by Tiberius [gabriel2024tiberius] | CPU, hours | per-species parameter sets shipped | Gaius-Augustus/Augustus |
| BRAKER1 / 2 / 3 | 2016 / 2021 / 2024 | 10.1093/bioinformatics/btv661 / 10.1093/nargab/lqaa108 / 10.1101/gr.278090.123 | EVID | RNA-seq (1), proteins (2), both (3) | 11-species benchmark | BRAKER3 raises average **transcript-level F1 by ~20 percentage points** over BRAKER1/2; largest gain on large complex genomes; beats MAKER2, Funannotate, FINDER [gabriel2024braker3] | CPU cluster, hours–days | benchmark assumes a stated level of proteome relatedness [gabriel2024braker3] | Gaius-Augustus/BRAKER |
| GeneMark-ETP | 2024 | 10.1101/gr.278373.123 | EVID | RNA-seq + proteins | large eukaryotic genomes | "significantly improves" over ETP predecessors [bruna2024genemarketp] | CPU | — | non-free |
| GALBA | 2023 | 10.1186/s12859-023-05449-z | EVID (protein-only) | DNA + protein db | many | see paper | CPU | miniprot+AUGUSTUS | Gaius-Augustus/GALBA |
| Gnomon (NCBI) | — | **no methods paper** | EVID | alignments (Splign/ProSplign/miniprot), RNA-seq | RefSeq organisms | not independently published | NCBI production | described only via [kuhn2013refseq; sayers2026ncbi] | inside egapx |
| EGAPx | — | **no methods paper** | EVID pipeline | assembly + RNA-seq + proteins | **Verified from the README** [egapx_readme]: supported taxa are Chordata, Arthropoda, Echinodermata, Mollusca, Cnidaria, monocots and eudicots; "Fungi, protists and nematodes are out-of-scope" | not independently published | **Verified**: prerequisites are Docker/Singularity plus "AWS Batch, SLURM/UGE cluster, or a r6a.4xlarge machine (32 CPUs, 256GB RAM)". Published runtimes on AWS Batch: *Drosophila* 144 Mb + 1 RNA-seq run = **71 CPU-hrs / 3 wall-hrs**; chicken 1.1 Gb + 20 RNA-seq runs = **425 CPU-hrs / 5.5 wall-hrs** [egapx_readme]. The charter's Galaxy figures remain unverified (source repo 404, §3) | ncbi/egapx |
| Helixer | 2021 / 2026 | 10.1093/bioinformatics/btaa1044 / 10.1038/s41592-025-02939-1 | DL (CNN+bLSTM) + HMM post-processor | DNA only | 2021: one vertebrate model over 186 animal genomes, one land-plant model over 51 plant genomes. 2026: fungal, plant, vertebrate, invertebrate | 2021: predictions "much less sensitive to genome length" than the then state of the art; outputs base-wise probabilities, not complete gene models [stiehler2021helixer]. 2026: "on par with or exceeding current tools", pretrained models usable without retraining [holst2026helixer] | GPU | **The first serious cross-species claim.** One model, many genomes | weberlab-hhu/Helixer, GPL-3.0 |
| Tiberius | 2024 | 10.1093/bioinformatics/btae685 | DL end-to-end (CNN + LSTM + differentiable HMM) | DNA only | trained on mammals; evaluated on human + 2 | **gene-level F1 62% on human vs 21% for the next best ab initio**; exon-intron structure of 2 of 3 human genes exactly right in de novo mode; ab initio accuracy *matches BRAKER3*, which uses RNA-seq + a protein database [gabriel2024tiberius] | **human genome in under 2 hours**; "fastest state-of-the-art"; GPU [gabriel2024tiberius] | 2024 version: mammals only | Gaius-Augustus/Tiberius, MIT |
| Tiberius, multi-clade | 2026 | 10.64898/2026.04.24.720536 | DL | DNA only | lineage-specific models for Mesangiospermae, Fungi, Vertebrata, Insecta, Chlorophyta, Bacillariophyta → **92% of available eukaryotic assemblies**. **Verified in the installed tool**: `--list_cfg` lists 9 configs over exactly those 6 clades (§3.1) | across **33 species**: gene-level F1 **+12 to +37 points over Helixer**, **+10 to +22 over ANNEVO**; approaches BRAKER3 in Mesangiospermae, Fungi, Bacillariophyta, Chlorophyta while being **~80× faster on GPU**; backend rewrite cut runtime 31% [gabriel2026tiberiusclades] | GPU | Strong — but note it is *six lineage-specific models*, not one general model | same |
| ANNEVO | 2026 | 10.1038/s41592-026-03036-7 | DL | DNA | multiple | beaten by Tiberius-multiclade by 10–22 F1 points at gene level [gabriel2026tiberiusclades]; not OA, so only its own abstract claim of "highly accurate ab initio" is available to me | — | — | xjtu-omics/ANNEVO (158 stars, pushed 2026-08-11) [zhang2026annevo] |
| TOGA | 2023 | 10.1126/science.abn3107 | COMP (alignment projection + orthology) | whole-genome alignment chains + reference annotation | 488 placental mammals, 501 birds | improves ortholog detection and annotation of conserved genes vs. state of the art; handles fragmented assemblies; also yields a genome-quality measure [kirilenko2023toga] | scales to hundreds of genomes | **Only for what is conserved relative to a reference** — cannot find clade-specific or fast-evolving genes | hillerlab/TOGA, MIT |
| CAT | 2018 | 10.1101/gr.233460.117 | COMP (HAL projection) | Cactus/HAL alignment + reference annotation | clades, personal genomes | see paper | cluster | same limitation as TOGA | repo **unmaintained**, 0 commits/12 mo |
| Ka/Ks ratio test | 2002 | 10.1101/gr.200901 | COMP (single statistic) | one pairwise alignment (human/mouse) | human/mouse | **false-negative rate lower than most current gene prediction methods and false-positive rate lower than all of them**, at the time; especially good on long exons and single-exon genes, which were then the hard cases [nekrutenko2002kaks] | trivial | the charter's floor | — |
| SpliceAI | 2019 | 10.1016/j.cell.2018.12.015 | DL (CNN, 10 kb context) | pre-mRNA sequence | human | accurate splice-junction prediction; variant-level validation against RNA-seq in 21 of 28 patients [jaganathan2019spliceai] | GPU | human-trained; cross-species use is off-label | Illumina/SpliceAI — **archived** |
| Pangolin | 2022 | 10.1186/s13059-022-02664-4 | DL | sequence, multi-tissue/species | 4 species | see paper | GPU | multi-species splicing | active |
| minisplice | 2026 | 10.1186/s13015-025-00293-7 | DL (1D-CNN) | sequence | vertebrates + insects | **7,026 parameters**; captures splice signals conserved across phyla; reveals mammal/bird-specific GC-rich introns; improves junction accuracy in minimap2/miniprot for noisy long reads and distant-homology proteins [yang2026minisplice] | trivial | **trained across phyla with one tiny model** | lh3/minisplice |
| miniprot | 2023 | 10.1093/bioinformatics/btad014 | alignment | protein + genome | any | comparable accuracy to prior protein-to-genome aligners, **tens of times faster** [li2023miniprot] | CPU, light | — | lh3/miniprot, MIT |
| Evo 2 | 2026 | 10.1038/s41586-026-10176-5 | DL foundation model | DNA | all domains of life; 9 Tbp training; 1M-token context | interpretability analysis shows learned representations of **exon-intron boundaries** among other features; no gene-structure benchmark reported [brixi2026evo2] | very large | broad by construction | ArcInstitute/evo2 |
| AlphaGenome | 2026 | 10.1038/s41586-025-10014-0 | DL | DNA, 1 Mb context | human/mouse regulatory | regulatory variant effects; not gene structure [avsec2026alphagenome] | large | — | google-deepmind/alphagenome |
| AlphaFold-3 gene-model scoring | 2026 | 10.1093/nar/gkag369 | DL (structure as a QC signal) | predicted proteins from gene models | F. graminearum, T. gondii, A. fumigatus | AlphaFold-3 scores support **65–84% of manually curated changes**; combining AF3 + Foldseek is most discriminative; the far cheaper Protenix-Mini retains the same discriminatory power [davison2026alphafoldannotation] | GPU (or cheap, with Protenix-Mini) | 3 fungal/protist species | — |
| Review (independent) | 2025 | 10.1093/bioadv/vbaf222 | — | — | — | extends the G3PO benchmark over AUGUSTUS, GENSCAN, GeneID, GlimmerHMM, SNAP + a gene-model-free NN + Helixer [djossou2025overview] | — | proposes a gene-model-based / -free / hybrid taxonomy | benchmark repo linked in paper |
| OMArk | 2025 | 10.1038/s41587-024-02147-w | QC | annotated proteome | any | detects erroneous gene inference in existing annotations [nevers2025omark] | light | — | active |
| BUSCO | 2021 | 10.1093/molbev/msab199 | QC | assembly or proteome | any | completeness, not correctness [manni2021busco] | light | — | active |

Supporting infrastructure cited but not a gene predictor: Progressive Cactus
[armstrong2020cactus] (>600 amniote genomes aligned, reference-free), phyloP
/ phastCons [pollard2010phylop] (36-mammal power analysis; UCSC conservation
tracks), Splign [kapustin2008splign], StringTie [pertea2015stringtie],
Liftoff [shumate2021liftoff], Zoonomia [christmas2023zoonomia].

### 2.1 The 2025–2026 cohort, found in run 3

The preprint pass (search-log row 10) turned up eleven methods that runs 1
and 2 missed entirely. This is not a long tail. Two of them — OrionGeno and
Vipsania — speak directly to the charter's hypothesis, and one of them
(§6.2) partly pre-empts it. I record the whole cohort here because the
synthesis (T-human-006) and the design proposal (T-human-011) both need it,
and because the *rate* at which this cohort is appearing is itself a finding:
five of the eleven were posted in the last five months.

| Method | Date | DOI | Class | What it is | Why it matters here |
|---|---|---|---|---|---|
| **OrionGeno** | 2026-04-29 | 10.64898/2026.04.26.720859 | DL, **phylogeny-aware** | End-to-end eukaryotic annotation predicting exons, introns, UTRs *and repeats* from sequence; "integrates phylogenetic context, long-range sequence modeling and joint prediction"; applied to **>5,300 unannotated NCBI chromosome-level genomes**; reports beating state of the art at exon, gene, protein-sequence and protein-structure level across lineages [liu2026oriongeno] | **The closest published thing to what the charter proposes**, from BGI. See §6.2. Licence is **non-commercial** (README badge), which matters for T-human-007 |
| **Vipsania** | 2026-08-30 | 10.64898/2026.08.26.747235 | DL, **unsupervised** | Stanke lab. Differentiable HMM layer inside a masked-language-model sequence network; *never shown a reference annotation*; pretrained pan-eukaryotically and finetuned unsupervised on the target genome; handles non-standard genetic codes. Claims to be **on average more accurate than supervised methods across most clades**, and to avoid the accuracy drop supervised models suffer on distant targets [krieg2026vipsania] | Attacks the exact failure mode §5 identifies (supervised models degrade with phylogenetic distance) by removing the labels. MIT, on PyPI. Directly relevant to the "no per-clade retraining" clause of the charter goal |
| **TOGA2** | 2026-07-04 | 10.64898/2026.06.30.735536 | COMP | Exon-level orthology and exon-wise annotation: **513× less memory, 6.1× faster** than TOGA. Adds gene-tree reconciliation and UTR annotation. Reports that **human-trained deep splice-site models generalize across vertebrates**, and uses them to handle splice-site shifts, intron deletions and exonization [malovichko2026toga2] | The generalization claim is an independent replication of the minisplice/SpliceAI result at a different scale, and it is the strongest evidence in this review that *splice signals* are the transferable part |
| **GeneCAD** | 2025-11-03 | 10.1101/2025.10.31.685877 | DL + CRF | PlantCAD2 foundation-model embeddings + transformer encoder + **chromosome-scale CRF** enforcing splice phase and feature order; sequence-only. Reports **~9% transcript-F1 over Helixer and BRAKER3** on angiosperms including an allotetraploid, and 86% recovery of classical CDS [liu2025genecad; zhai2025plantcad2] | The CRF-as-grammar design is the main published alternative to a differentiable HMM. Apache-2.0. Also a cautionary tale: v0.1.0 shipped with a bug that wrecked BUSCO scores (§5.9) |
| **GENATATORs** | 2026-06-21 | 10.64898/2026.06.17.732686 | DL benchmark + method | Systematic study of DNA-LM gene segmentation. Two results we should not ignore: **pretrained DNA-LM embeddings do not capture the features needed for gene segmentation** (task-specific finetuning is essential), and **standard per-token / per-sequence metrics fail to capture real annotation quality**; proposes biologically grounded metrics and datasets [shmelev2026genatators] | A direct, negative result about the "just use a foundation model" path, and a metrics critique that T-human-007 should adopt rather than re-derive |
| **He & Florea benchmark** | 2026-02-23 | 10.64898/2026.02.22.707219 | benchmark | Evaluates SegmentNT, Enformer, Borzoi (with segmentation heads), SpliceAI and AlphaGenome on **stratified** exon classes: coding vs non-coding, terminal vs internal, constitutive vs alternative, TE-derived. Finding: every method is best on the exon class in its training data and **degrades drastically on under-represented classes** [he2026benchmarkfm] | The single most useful benchmark-design input I found. Aggregate F1 hides this entirely; §5.8's complaint now has a citation |
| **SegmentNT** | 2024-03-15 | 10.1101/2024.03.14.584712 | DL | Frames annotation as **instance segmentation**; finetunes Nucleotide Transformer to segment 14 genic and regulatory element classes at single-nucleotide resolution [dealmeida2024segmentnt] | The formulation GeneCAD and OrionGeno both build on; CC-BY-NC-ND |
| **geneML** | 2026-05-21 | 10.64898/2026.05.18.725946 | DL | Fungal-specific; gene-level F1 **64.9 → 67.1 vs BRAKER3 with protein hints** (recall 64.1 → 69.0 at equal precision) across nine fungal genomes, **~6 min/genome on 8 CPU cores**, and predicts **alternative transcripts** (41.1% recall / 71.1% precision vs Iso-Seq, against AUGUSTUS's 33.8 / 48.9) [vader2026geneml] | Fungi are explicitly out of scope for EGAPx. This is a CPU-only tool beating an evidence-based pipeline there, and one of very few with isoform numbers |
| **ANNEVO** (first-hand) | 2025 preprint / 2026 | 10.21203/rs.3.rs-6402260/v1 → 10.1038/s41592-026-03036-7 | DL | Now read from its own abstract rather than a competitor's table: a **mixture-of-experts genomic language model** modelling distal dependencies and "joint evolutionary relationships", benchmarked on **566 phylogenetically diverse species**; claims to exceed reference annotations for some species [zhang2026annevo] | The 566-species evaluation is the largest species panel in this review and a candidate template for T-human-007 |
| **GeMoSeq** | 2026-02-01 | 10.1093/nar/gkag091 | EVID | Transcript reconstruction from RNA-seq by combinatorial enumeration plus likelihood-based quantification, with CDS prediction integral to the algorithm; benchmarked over seven species [grau2026gemoseq] | The RNA-seq-side baseline; relevant to T-human-008's evidence inventory |
| **OMAnnotator** | 2026-01-22 | 10.1093/bioadv/vbag015 | EVID (combiner) | Repurposes the OMA orthology algorithm to build a **consensus** from ab initio, transcriptomic and homology annotations, using evolutionary information as the tie-breaker; improves on its own sources on *D. melanogaster* [bates2026omannotator] | Evolutionary geometry used as an arbiter rather than as an input — a cheap design worth knowing before we propose an expensive one |

Two further items that are not gene finders but that Phase 2 should have:
**GAP-MS** [abbas2026gapms], which validates gene models against mass-spec
peptides across nine crops and finds hundreds of peptide-supported loci
missing from reference annotations; and the **Pristionchus pacificus**
curation study [roedelsperger2026pristionchus], in which community curation
corrected **more than 7,500 gene models — about 24% of the annotation** of
one nematode strain. Both are direct evidence about the quality of the
"ground truth" T-human-007 intends to score against.

Full BibTeX: `refs.bib` (64 entries as of run 3; every publication row
carries a DOI, and no entry is marked `TO VERIFY`).

---

## 3. Repository inventory

Full table in `repos.tsv` (32 rows, GitHub API, 2026-09-09). Columns: URL,
last commit, commits in trailing 12 months, open issues, stars, language,
SPDX licence, archived, install tested, notes.

What the numbers say:

**Alive and moving.** `ComparativeGenomicsToolkit/cactus` (485 commits/12 mo,
704 stars) and `EBI-Metagenomics/genomes-pipeline` (189) are the busiest.
Among gene predictors, **Tiberius is the only one under real development**:
174 commits in 12 months, last push 2026-09-07, MIT, 138 stars, 21 open
issues. `Helixer` is second (49 commits, GPL-3.0). `guigolab/geneid` is a
surprise at 105 commits — a 1990s-lineage GHMM in C still being worked on.

**Alive but slow.** `ncbi/egapx` 26 commits, 207 stars, licence
`NOASSERTION`. `Gaius-Augustus/BRAKER` 11 commits but **103 open issues** and
467 stars — heavy use, thin maintenance. `Gaius-Augustus/Augustus` 5 commits
and **175 open issues**. The BRAKER/AUGUSTUS stack is load-bearing for the
field and is being maintained at a fraction of the rate it is being used.

**Effectively dead, still widely used.** `Comparative-Annotation-Toolkit` 0
commits/12 mo with 116 open issues. `agshumate/Liftoff` last pushed
2023-08-01, 552 stars, 78 open issues. `KorfLab/SNAP` last pushed 2022 and
still a MAKER dependency. `Illumina/SpliceAI` is **archived** with 507 stars
— the field's de facto splice-site model is a read-only repository.
`hillerlab/TOGA` 1 commit/12 mo despite being the engine behind the largest
comparative gene resources published [kirilenko2023toga].

**Licensing.** MIT: Tiberius, miniprot, TOGA, StringTie, gffcompare.
GPL-3.0: Helixer, Liftoff. Apache-2.0: CAT, evo2, alphagenome, borzoi.
`NOASSERTION` (i.e. GitHub could not resolve a standard licence): egapx,
BRAKER, GALBA, SNAP, SpliceAI, miniprot's deps, nucleotide-transformer.
AUGUSTUS resolves to no licence at all. **The GeneMark components inside
BRAKER have historically been under a non-free academic licence
[unverified — confirm the current terms next run].** If true, that has direct
consequences for whether BRAKER3 can be a redistributable benchmark baseline,
and T-human-007 should know before it designs around it.

### 3.1 Install attempts (run 2)

Fresh virtual environments, Python 3.13.9, gcc 15.2.0, Ubuntu, RTX 5080
(compute capability 12.0), 26 GB free on the build volume. Each tool was
installed by following its own README verbatim, not by improvising.
Full column in `repos.tsv`.

| Tool | Commit | Steps followed | Result | Wall time | Footprint |
|---|---|---|---|---|---|
| miniprot | 81f9b93 | `git clone; make` | **works**, binary reports `0.18-r281` | seconds | tiny, zlib only |
| minisplice | 49f9e8c | `git clone; make` | **works**, `gentrain`/`train`/`inspect` subcommands run | seconds | tiny, zlib only |
| Tiberius | e73844b | README quick-start: `pip install .` | **installs, does not run** — see below | 3 s | 30 MB |
| Tiberius | e73844b | `pip install '.[from_source]'` | **works**; TF 2.20.0 imports, GPU visible | 1 min 46 s | **6.3 GB** |
| egapx | f9a7392 | `venv; pip install -r requirements.txt` | **installs** (`requirements.txt` is one line: PyYAML), `ui/egapx.py -h` runs. **Cannot be executed here** | 2 s | <1 MB |
| Helixer | d17bb49 | README read only | **not attempted** — see below | — | — |
| BRAKER | — | README read only | **not attempted** — see below | — | — |

Four findings that matter to T-human-007 and T-human-009:

**(a) Tiberius's documented quick-start does not produce a working tool.**
The README's `pip install .` installs only the launcher (`rich`, `pyyaml`);
invoking it then fails with `ModuleNotFoundError: No module named
'packaging'`. The inference dependencies live in a `from_source` extra
(`bricks2marble[tf]`, `tensorflow[and-cuda]>=2.17,<2.21`, biopython, pandas),
and the README's actually-recommended path is Singularity. `pip install
'.[from_source]'` then works cleanly. This is a small documentation gap, not
a defect, but anyone scripting a baseline will hit it.

**(b) Tiberius cannot use a current-generation consumer GPU without a long
JIT.** On the RTX 5080, TensorFlow 2.20 emits: *"TensorFlow was not built
with CUDA kernel binaries compatible with compute capability 12.0. CUDA
kernels will be jit-compiled from PTX, which could take 30 minutes or
longer."* The pin is `tensorflow>=2.17,<2.21` in `pyproject.toml`, so this is
not something a user can resolve by upgrading TensorFlow. The charter budgets
Phases 1–3 for "a laptop or one consumer GPU"; **the newest consumer GPUs are
the worst case for this baseline**, and T-human-009 must record which GPU
generation each timing came from or its numbers will not be comparable.

**(c) The full Tiberius environment is 6.3 GB.** For a model the paper
describes as small and fast, essentially all of that is TensorFlow and CUDA.
This is an argument for the charter's design goal that is stronger than the
accuracy argument: the *model* is small; the *stack* is not.

**(d) Two of the three headline tools are container-first by their authors'
own instruction.** Helixer's README states installation takes "20-30 minutes"
for an experienced user and "a maximum of 2-3 hours" for an inexperienced
one, recommends Docker/Singularity, and restricts manual installation to
Linux [helixer_readme]. BRAKER's README warns that conda installs of
GeneMark-ETP have caused "multiple problems reported by users" and directs
users to the Singularity image. I did not install either: BRAKER additionally
requires GeneMark-ETP, AUGUSTUS, ProtHint, StringTie2, bedtools, GffRead and
partitioned OrthoDB clades as separate downloads, which is a multi-hour job
and belongs in T-human-009 with a measured stopwatch, not here.

**(e) A licence finding with consequences.** BRAKER's own scripts are under
the Artistic License (README §License). But GeneMark-ETP — the gene finder at
the core of BRAKER3 — is at `gatech-genemark/GeneMark-ETP` **with no LICENSE
file** [genemark_etp_repo], so its redistribution terms are unstated and
default to all-rights-reserved. Run 1 guessed "historically non-free"; the
verified position is worse in one way (unstated, not merely restricted) and
better in another (it is now on GitHub rather than key-file gated).
**T-human-007 should decide early whether BRAKER3 can be a redistributable
baseline or only a locally-run comparator.**

**(f) EGAPx ships a security notice.** Its README states that static analysis
has found "a small number of verified buffer overrun security vulnerabilities"
in its NCBI C++ toolkit dependencies and recommends running it in a VM or
cloud instance [egapx_readme]. Worth knowing before anyone runs it on shared
infrastructure for T-human-009.

**Repositories run 1 could not reach — resolved and unresolved:**

- `nekrut/axomeme` → **404**. The charter names it as the HyphAeon design
  pattern this whole project is meant to transfer (2M parameters, 2D axial
  transformer, 4D MDS tree embeddings, Tree-RoPE). I cannot read it.
- `nekrut/scalingPaper` → **404**. The charter's EGAPx cost figures (416k
  CPU-hours, 1,409 jobs, 290 genomes, 45% failure, 24% CPU efficiency) come
  from it. I have restated them as charter claims, not as verified ones.
- `gmod/maker` → 404. **Resolved in run 2**: there is no canonical MAKER
  repository on GitHub under any org. MAKER is a registration-gated tarball
  from yandell-lab.org. Its practical successor is
  `nextgenusfs/funannotate` (BSD-2-Clause, 400 stars, pushed 2026-09-08),
  which is also one of the pipelines BRAKER3 was benchmarked against
  [gabriel2024braker3].

Both `nekrut/*` repositories are private, renamed, or deleted from the
perspective of an unauthenticated client. I raised this to the coordinator as
a `question` in run 1 (`20260909T012944Z-lenin-0002`); it is unanswered as of
run 2. It affects every agent's Phase 1 task, not only mine, and it blocks
Phase 3 more than Phase 1. **Partial mitigation found in run 2**: the EGAPx
README supplies its own hardware floor and per-genome runtimes
[egapx_readme], so T-human-009 is no longer wholly dependent on
`scalingPaper` for the EGAPx cost story — only for the Galaxy failure-rate
and efficiency figures, which have no substitute I have found.

`salzberg-lab/GlimmerHMM` (Elixir, 0 stars) is **not** canonical, and run 2
found no repository that is: `kblin/glimmerhmm` (C, 5 stars) was last pushed
in 2013 and is also a mirror. GlimmerHMM is distributed from ccb.jhu.edu.
Recorded as tarball-distributed and unmaintained on GitHub.

The Ensembl repository is **resolved**: `Ensembl/ensembl-anno`, "Ensembl
Automatic annotation pipeline", Python, Apache-2.0, pushed 2026-09-08. It is
the **third production annotation pipeline** alongside EGAPx and BRAKER and
was missing from run 1 entirely; T-human-009's cost baseline should cover all
three, not two.

Five repositories were added to `repos.tsv` in run 2: `xjtu-omics/ANNEVO`,
`lh3/minisplice`, `nextgenusfs/funannotate`, `Kuanhao-Chao/OpenSpliceAI` (the
maintained successor now that `Illumina/SpliceAI` is archived), and
`gatech-genemark/GeneMark-ETP`. ~~Their 12-month commit counts and issue
counts are missing because the GitHub `/repos` endpoint was rate-limited.~~
**Backfilled in run 3** using an authenticated client.

### 3.2 The run-3 cohort, and what their READMEs say

Four repositories were added in run 3, one of them the most important row in
the file. All figures below are from the repository itself, not from a paper.

| Repo | Pushed | Commits/12 mo | Stars | Licence | The thing to know |
|---|---|---|---|---|---|
| `BGIResearch/OrionGeno` | 2026-08-21 | 16 | 24 | **`NOASSERTION`; README badge reads "License: Non-Commercial"** | Weights on Hugging Face and ModelScope; a hosted API at CNGBdb. Requires Linux, Python `>=3.10,<3.11`, and an **NVIDIA GPU of compute capability ≥ 7.0** with `mamba-ssm`/`causal-conv1d`. So the model whose design is closest to ours is also the one we may not be able to use as a baseline in a redistributable benchmark — a T-human-007 problem that should be raised now, not later |
| `Gaius-Augustus/Vipsania` | 2026-09-07 | 11 | 17 | **MIT** | `pip install vipsania` (PyPI 1.0.0, Python ≥ 3.12) pulls `bricks2marble[tf]` and `tensorflow<2.20`; `hidten` supplies the differentiable HMM. GPU strongly recommended; CPU annotation "possible, but slow". Annotates from a FASTA alone, with `--finetune` on the target genome itself |
| `hillerlab/TOGA2` | 2026-08-19 | 100+ | 53 | **MIT** | Active successor to `hillerlab/TOGA` (211 stars, 1 commit in 12 months — TOGA itself is now effectively frozen) |
| `plantcad/genecad` | 2026-09-08 | 100+ | 40 | **Apache-2.0** | The most actively developed repo in the whole inventory. Ships plant *and* vertebrate models, Docker images, CI. v0.4.0 added enforcement of canonical start/stop and donor/acceptor motifs and minimum intron/exon lengths — i.e. the grammar constraints were bolted on *after* the neural model, which is worth noting for §6 |

Licensing summary across the inventory, since it decides what T-human-007 can
actually run: MIT or Apache for Tiberius, TOGA/TOGA2, miniprot, Vipsania and
GeneCAD; GPL-3.0 for Helixer and OpenSpliceAI; **no LICENSE file at all** for
`gatech-genemark/GeneMark-ETP` (and hence unsettled redistribution for
BRAKER3) and for `lh3/minisplice`; **non-commercial** for OrionGeno;
`NOASSERTION` for ANNEVO and egapx. Only the first group is unambiguously
safe to redistribute in a benchmark image.

---

## 4. Data sources noted in passing

Not an inventory — that is T-human-008 — just what these papers actually
consumed, so the inventory task starts from evidence rather than a list.

- **Whole-genome multiple alignments.** Progressive Cactus, reference-free,
  demonstrated at >600 amniote genomes [armstrong2020cactus]; HAL is the
  container CAT projects through [fiddes2018cat]. UCSC multiz/Cactus
  alignments and the chain/net files TOGA requires [kirilenko2023toga].
- **Conservation tracks.** phyloP and phastCons, with the 36-mammal power
  analysis and the UCSC 44-vertebrate Conservation tracks
  [pollard2010phylop]. Note the paper's own result: enough power for strong
  selection at single nucleotides, moderate selection in 3-bp elements, and
  weaker or clade-specific selection only in longer elements. That is
  directly the resolution question for splice-site scale features.
- **Zoonomia**: 240+ placental mammal alignments and constraint calls
  [christmas2023zoonomia]; the substrate TOGA's 488-mammal resource was
  built on.
- **Reference annotations.** NCBI RefSeq [kuhn2013refseq; sayers2026ncbi];
  Ensembl; the G3PO benchmark set used by the 2025 review
  [djossou2025overview], which is the only ready-made, published,
  gene-prediction-specific benchmark I found this run.
- **Expression.** StringTie-assembled transcriptomes [pertea2015stringtie]
  are what BRAKER3 and EGAPx actually consume, not raw reads.
- **Proteins.** OrthoDB partitions are BRAKER2/3's protein input; the
  benchmark in [gabriel2024braker3] is explicitly parameterized by *how
  related* the available proteome is — a leakage axis T-human-007 must
  control, because "protein evidence" silently smuggles in the answer when
  the informant is close.

---

## 5. Failure modes of current tools

**5.1 Per-species fitting is the field's founding assumption.** Korf states
it as the conclusion of the SNAP paper: "every genome needs a dedicated gene
finder", after showing that foreign gene finders are highly inaccurate and
that the nearest phylogenetic neighbour is not reliably the best parameter
donor [korf2004snap]. GeneMark-ES removed the *labelled training set*
requirement via self-training [lomsadze2005genemarkes] but not the per-genome
fit. Helixer was the first to claim a single cross-species model
[stiehler2021helixer], and Tiberius 2026 still ships **six lineage-specific
models** rather than one [gabriel2026tiberiusclades]. Twenty-two years after
Korf, the assumption has been dented, not overturned. That is the gap this
project exists to close, and it is the right gap.

**5.2 Clade exclusion is explicit, not accidental.** EGAPx declares fungi,
protists and nematodes out of scope — **verified in run 2 directly from the
repository README**, which carries it as a warning and names the supported
taxa as Chordata, Arthropoda, Echinodermata, Mollusca, Cnidaria, monocots and
eudicots [egapx_readme]. That list is animals and flowering plants. It
excludes every fungus, every protist, every nematode, and all
non-angiosperm plants. Tiberius before 2026 was mammals-only
[gabriel2024tiberius]. Even the 2026 extension reaches "92% of currently
available eukaryotic assemblies" [gabriel2026tiberiusclades] — and the
remaining 8% is precisely the tail (protists, early-branching lineages) that
a species-independent method would be judged on. **A benchmark drawn from the
92% cannot measure species independence.** T-human-007 should deliberately
oversample the excluded tail.

**5.3 Genome size and intron length.** Helixer's headline generalization
claim is specifically that its predictions are "much less sensitive to the
length of the genome" than the then state of the art [stiehler2021helixer] —
which tells you that length sensitivity was the dominant failure mode of GHMM
methods. BRAKER3's gains are "most pronounced for species with large and
complex genomes" [gabriel2024braker3], same signal from the other direction:
the large-intron regime was where everything was breaking.

**5.4 GC composition.** GENSCAN's answer in 1997 was to fit *distinct
parameter sets per C+G compositional region* [burge1997genscan] — isochores
handled by partitioning the model, not by learning a representation. minisplice
finds GC-rich introns specific to mammals and birds [yang2026minisplice],
i.e. composition is still clade-entangled at the splice signal itself, in
2026. Any model claiming clade independence must show its splice-site
performance stratified by GC, or the claim is untested.

**5.5 Resource use.** The compute story is the clearest argument for this
project. BRAKER3 is a CPU-cluster pipeline; Tiberius annotates the human
genome in under two hours [gabriel2024tiberius] and is ~80× faster than
BRAKER3 on a GPU [gabriel2026tiberiusclades].

Run 2 verified EGAPx's costs from primary source rather than from the
charter. Its stated prerequisites are Docker or Singularity plus "AWS Batch,
SLURM/UGE cluster, or a r6a.4xlarge machine (32 CPUs, 256GB RAM)", and its
own published timings on AWS Batch are **71 CPU-hours / 3 wall-hours for a
144 Mb *Drosophila* genome with one RNA-seq run**, and **425 CPU-hours / 5.5
wall-hours for the 1.1 Gb chicken genome with 20 RNA-seq runs**
[egapx_readme]. Note the shape: CPU-hours grow ~6× while the genome grows
~7.6× and wall time grows less than 2×, because the cost is bought with
parallelism, not reduced. A fly genome costs 71 CPU-hours; Tiberius annotates
a 3 Gb human genome in under 2 hours on one GPU [gabriel2024tiberius]. **That
is the gap the charter is aiming at, and it is roughly two orders of
magnitude, measurable from published numbers alone.**

The charter's Galaxy figures (416k CPU-hours over 1,409 jobs, 45% failure
rate, 24% CPU efficiency) remain unverified — the source repository is
unreachable (§3) and I have found no substitute. If they hold, the
interesting quantity is not the CPU-hours but the *45% failure rate*: nearly
half the compute bought nothing. That is a robustness failure being reported
as a cost failure, and it should be framed that way in T-human-009. It is
also the one part of the cost story that the EGAPx README cannot supply,
since documented runtimes are by construction the runs that succeeded.

Run 2 adds a second-order point from the install bench (§3.1): the *stack*,
not the model, is where the weight is. Tiberius's full environment is 6.3 GB
of TensorFlow and CUDA, and on a current consumer GPU (compute capability
12.0) its pinned TensorFlow must JIT its kernels from PTX, which TensorFlow
itself warns "could take 30 minutes or longer". A tool whose selling point is
speed should not have a 30-minute first-run penalty on new hardware.

**5.6 Comparative methods only find what is already known.** TOGA and CAT
project annotation from a reference through an alignment
[kirilenko2023toga; fiddes2018cat]. They are excellent on conserved genes
and structurally cannot discover clade-specific ones. Any comparative design
we propose inherits this unless it is explicitly built to predict *without* a
reference annotation on the informant side.

**5.7 Maintenance risk is a real failure mode.** §3: the field's splice-site
model is archived, its comparative annotation toolkit is at zero commits, its
most-used ab initio pipeline has 103 open issues and 11 commits a year. When
T-human-009 measures baselines, some of them will not build. Budget for that.

**5.8 Evaluation is not standardized.** No two rows in §2 share a benchmark.
GENSCAN reports exact exons on 1997 vertebrate sets, BRAKER3 reports
transcript-level F1 on 11 species under an assumed proteome-relatedness
level, Tiberius reports gene-level F1 on human, Tiberius-multiclade on 33
species over 6 clades, the 2025 review uses G3PO
[djossou2025overview]. Comparing published numbers across these is
meaningless, and I have deliberately not built a "which is best" column.
**Constructing a single benchmark all of them run on is the highest-value
thing this project can do first**, independent of any model design.

Run 3 adds a citation this argument previously lacked. He and Florea evaluate
SegmentNT, Enformer, Borzoi, SpliceAI and AlphaGenome on *stratified* exon
classes and find that every method peaks on the exon class best represented
in its training data and "decreases drastically" on the rest — non-coding,
terminal, alternatively spliced and TE-derived exons [he2026benchmarkfm].
GENATATORs makes the complementary point from the metrics side: standard
per-token and per-sequence scores "fail to capture the challenges of
real-world gene annotation" [shmelev2026genatators]. **T-human-007 should
therefore require stratified reporting, not a single F1**, and should adopt
the metric critique from these two papers rather than re-deriving it.

**5.9 The reference annotations are themselves wrong, in measurable amounts.**
This is the failure mode I under-weighted in runs 1 and 2, and run 3 gives it
three independent numbers. (i) Vertebrate selenoprotein genes — where UGA is
recoded rather than a stop — are well annotated for only **11% of genes in
Ensembl and 5% in NCBI GenBank**, because neither pipeline has a dedicated
selenoprotein path [tico2026selenoprotein]. (ii) Community curation of one
*Pristionchus pacificus* strain identified and corrected **more than 7,500
gene models, about 24% of the annotation**, and attributed them to assembly
errors, artificial transcript fusions and unexpressed genes
[roedelsperger2026pristionchus]. (iii) GAP-MS finds **hundreds of
peptide-supported coding loci absent from reference annotations** across nine
crops [abbas2026gapms]; ANNEVO and OrionGeno make the same claim from the
prediction side [zhang2026annevo; liu2026oriongeno]. The consequence for us
is concrete: a model scored against RefSeq or Ensembl is partly being scored
on its ability to reproduce known errors, and above some accuracy level the
benchmark stops measuring the model. T-human-007 needs an error bar on the
ground truth, and a curated high-confidence subset to score on separately.

**5.10 Preprint-stage tools are not yet reliable artifacts.** GeneCAD's
README carries a warning that v0.1.0 "caused low BUSCO scores" and that users
must upgrade — a bug found by the MaizeGDB team after release, not by the
authors. Half of §2.1 is unpublished preprint software of the same maturity.
Any number T-human-009 takes from this cohort should be pinned to a version
and re-measured, not quoted.

---

## 6. Opinion: what I would build

One page, and I am willing to be wrong on all of it.

**The bet.** The charter's hypothesis is that gene structure, given alignment
+ tree + codon geometry as coordinates, has a low-dimensional decision
surface. I think the evidence for this is stronger than the charter claims,
and it comes from an unexpected place: **minisplice learns splice signals
conserved across phyla with 7,026 parameters** [yang2026minisplice]. Not 2
million — seven thousand. And Nekrutenko's Ka/Ks test achieved a false
positive rate below every gene finder of its day using *one statistic on one
pairwise alignment* [nekrutenko2002kaks]. Both results say the same thing:
when the comparative signal is presented in the right coordinates, the
remaining decision is nearly trivial. What is large in current models is not
the decision surface. It is the machinery for *manufacturing the coordinates*
from raw DNA — which is exactly the work Tiberius's CNN+LSTM stack and
Helixer's bLSTM are doing, and exactly the work a comparative model gets
handed for free.

**What I would build.** A small axial model over a *codon-aware multiple
alignment window*, in the HyphAeon pattern: site axis × taxa axis, tree
injected as a metric (MDS embedding + Tree-RoPE) rather than learned. Output
per alignment column: a distribution over {intergenic, 5'UTR, CDS-phase-0/1/2,
intron, 3'UTR} plus explicit donor/acceptor/start/stop head. Decode with a
differentiable HMM, following Tiberius — the end-to-end integration, not
Helixer's separate post-processor, is the part of Tiberius I would keep
[gabriel2024tiberius].

**Inputs.** Cactus/HAL alignment window for the target locus with its induced
subtree [armstrong2020cactus]; the target sequence; the informant sequences;
the tree. Explicitly **not** RNA-seq and **not** a protein database, so that
the comparison against BRAKER3 is honest and the method stays usable on a
genome with no transcriptome — which is most genomes.

**The geometry I would encode.** Three inductive biases, and I would defend
each as a *constraint*, not a feature: (1) **phase is modular arithmetic** —
CDS length ≡ 0 mod 3 and phase composes across exons; make it a structural
constraint in the decoder, not something a network rediscovers per clade.
(2) **The tree is a metric, not a sequence** — inject it, as HyphAeon does;
this is what should buy clade independence, because a model that never
memorizes *which* taxa it saw cannot overfit to a clade's taxon set.
(3) **Intron length is a nuisance parameter, not a signal** — the single
clearest cross-clade failure mode in §5.3. I would make the model
length-equivariant over intron spans (predict boundaries, not spans) so that
a 60 bp fungal intron and a 100 kb mammalian intron are the same object.

**Parameter budget.** I would target 2–5M and treat exceeding 20M as evidence
the geometry is wrong rather than as evidence more capacity is needed.
Tiberius is already accurate at modest size; the argument for going smaller
is not efficiency, it is that a small model that generalizes across clades is
a *scientific claim* about gene structure, and a large one is not.

**The biggest risk, and I think it is underrated.** Not accuracy —
**alignment availability**. The whole design assumes a good multiple alignment
with a tree at every locus in a novel genome. For a newly sequenced genome in
a sparsely sampled clade, that alignment either does not exist or is the most
expensive part of the pipeline, at which point we have moved the cost rather
than removed it, and we have re-introduced exactly the clade dependence we
set out to kill: the method works where relatives have been sequenced. Ka/Ks
worked in 2002 because human *and mouse* existed [nekrutenko2002kaks]. Every
comparative method since inherits that precondition. **Before committing to a
comparative design, Phase 2 must measure how accuracy degrades as informant
density and divergence degrade** — with the honest possibility that the
answer sends us to a Tiberius-style single-sequence model with comparative
signal as an optional refinement, rather than the reverse. I would make that
degradation curve a required experiment in T-human-007, and I would rather
find out there than in Phase 4.

Second risk, smaller but real: the field's evaluation is fragmented enough
(§5.8) that a new model can look excellent by choosing its benchmark. We
should fix the benchmark before we have a model, precisely so we cannot.

### 6.1 Addendum from run 2, which sharpens the above

Everything above §6.1 was written and committed in run 1. Installing Tiberius
turned up something I had not read in any paper, and it changes the question
I would put to Phase 3.

**Tiberius already has a comparative mode, and it is not the headline.** The
installed CLI exposes `--clamsa`, and the README documents a "de novo" mode
(as distinct from "ab initio") that consumes evolutionary information derived
from multiple sequence alignments by ClaMSA — sitewise predictions across all
six reading frames, converted to per-sequence NumPy arrays, with separately
trained weights. There is a shipped model configuration for it:
`mammalia_clamsa_v2`.

Two things follow. First, **the comparative-input design I argued for in §6
partly exists**, inside the strongest current tool, from the same group. Any
proposal we write that does not say how it differs from Tiberius-with-ClaMSA
is not a proposal.

Second, and more interesting: of the nine model configurations the installed
tool lists, **exactly one uses ClaMSA, and it is mammals-only**. Every
non-mammalian clade model — angiosperms, chlorophyta, diatoms, fungi,
insecta, vertebrates — is sequence-only. The 2026 multi-clade paper's headline
is ab initio accuracy [gabriel2026tiberiusclades]. So the group with the
best-performing architecture, having built the comparative path, extended the
*non*-comparative one to six clades.

That is either (i) evidence that comparative input buys less than the
charter's hypothesis assumes once a good sequence model exists, or (ii)
evidence that generating ClaMSA input per clade is too expensive or too
data-hungry to scale — which is exactly my "alignment availability" risk
above, showing up as revealed preference rather than as an argument. Both
readings are bad news for a naively comparative design, and they point at
different fixes.

**I would make resolving this the first question of T-human-011**, ahead of
any architecture work: obtain or reproduce the ab initio vs. ClaMSA-mode
comparison on mammals, and find out which of (i) or (ii) is true. It is a
cheap experiment — the models and the code are both installed and MIT
licensed — and it is decision-relevant in a way that no amount of further
literature reading is. If the answer is (i), the charter's central hypothesis
needs revision before we build anything.

### 6.2 Addendum from run 3: someone has built a version of this

Run 2's addendum sharpened the question. Run 3 changes the answer, and I
would rather say so plainly than defend the run-1 opinion.

**OrionGeno is a phylogeny-aware deep model for end-to-end eukaryotic
annotation, and it has already been run on more than 5,300 genomes**
[liu2026oriongeno]. The abstract's own framing — "integrates phylogenetic
context, long-range sequence modeling and joint prediction of gene structures
and repetitive elements", motivated by existing methods' failure to
"generalize across distant lineages" — is, sentence for sentence, close to
the charter's motivation. ANNEVO, separately, is a mixture-of-experts genomic
language model that claims to model "joint evolutionary relationships" across
566 species [zhang2026annevo]. The bet that phylogeny-as-input is the missing
inductive bias is no longer an open bet; it is being placed, by well-resourced
groups, and at least two of them report state-of-the-art results.

This does not make the charter's project pointless, but it changes what its
contribution can honestly be. Three things are still genuinely open:

1. **Nobody has shown *which* geometry does the work.** OrionGeno reports
   accuracy; it does not isolate the contribution of phylogenetic context
   from that of long-range modelling, and neither does ANNEVO. HyphAeon's
   claim is specifically that injecting the tree *as a metric* (MDS
   embeddings, Tree-RoPE) removes the need to learn phylogeny. An ablation
   that separates tree-as-metric from tree-as-extra-tokens from no-tree,
   at a fixed parameter budget, would be a real result, and it is cheap.
2. **Parameter budget is untouched territory.** OrionGeno needs an NVIDIA GPU
   with `mamba-ssm`; GeneCAD sits on a foundation model. minisplice does a
   real job with **7,026 parameters** [yang2026minisplice], and TOGA2 finds
   that human-trained splice models transfer across vertebrates
   [malovichko2026toga2]. The charter's ~2M-parameter target is the part of
   the hypothesis nobody in this cohort is testing.
3. **Licensing and reproducibility.** OrionGeno is non-commercial and its
   weights sit behind two model hubs; the benchmark this project builds can
   be the neutral ground none of these tools currently has.

And Vipsania forces one more revision. Its claim — an unsupervised
differentiable-HMM gene finder that is *on average more accurate than
supervised methods across most clades* and does not degrade with phylogenetic
distance [krieg2026vipsania] — attacks the generalization problem from the
opposite side from the charter: not by giving the model better geometry, but
by removing the labelled data whose clade bias causes the degradation. If
that result holds, the diagnosis in §5.1 is right but the charter's remedy is
not the only one, and a design proposal that ignores it is incomplete.

**Revised recommendation for T-human-011.** Keep the §6 design, but demote it
from "the idea" to "one of three arms", and put the first run's effort into
the comparison rather than the architecture:

- Arm A: sequence-only at ~2M parameters (Tiberius-class, our budget).
- Arm B: A plus tree-as-metric geometry (the HyphAeon transfer).
- Arm C: A trained unsupervised, Vipsania-style.

with OrionGeno, Tiberius-multiclade, Helixer, GeneCAD and BRAKER3 as external
reference points on the T-human-007 benchmark. The deliverable that no one
else is producing is **the ablation and the neutral benchmark**, not another
architecture. I would still run the run-2 experiment (§6.1, ab initio vs.
ClaMSA on mammals) first, because it is a day's work and it decides whether
Arm B is worth building at all.

I was wrong in run 1 to treat the comparative-geometry idea as unoccupied
ground. It is occupied. The unoccupied ground is the measurement.

---

## 7. State of this task

Not part of the deliverable format; recorded so the work is resumable.

**Run 1** produced the search log, publications table, `refs.bib`,
`repos.tsv`, failure modes and the §6 opinion.

**Run 2** closed the install gap for the four tools that matter, verified
EGAPx's hardware floor, clade exclusions and per-genome runtimes from primary
source, resolved every repository run 1 could not find and added five it had
missed, verified both `TO VERIFY` DOIs, settled the GeneMark licence
question, and found that Tiberius already ships a comparative (ClaMSA) mode
that is mammals-only while its six-clade extension is sequence-only (§6.1).

**Run 3** (this one) ran the preprint and OpenAlex passes that §1.1 had
flagged as the remaining substantial gap, and they were not a formality: they
produced §2.1, eleven 2025–2026 methods that runs 1 and 2 missed, five of
them posted in the last five months. Consequences recorded in §3.2 (four new
repositories, and a licensing summary of what a benchmark may actually
redistribute), §5.9 and §5.10 (two new failure modes, both about the quality
of the ground truth and of the tools themselves), and §6.2, which revises the
run-1 opinion: the comparative-geometry idea is occupied ground, the
measurement is not. ANNEVO is now read first-hand from its own abstract. The
five run-2 repository rows with `?` counts are backfilled.

**What is still not done, stated so the synthesis can weigh it:**

- **arXiv is not covered directly.** Its API returned HTTP 429 to this host on
  every attempt, including a single isolated query. OpenAlex `type:preprint`
  partly compensates; Semantic Scholar was not attempted.
- **No install was attempted for the run-3 cohort.** Vipsania (MIT, PyPI),
  GeneCAD (Apache-2.0, Docker) and TOGA2 (MIT) are all installable and the
  §3.1 method would apply; OrionGeno needs a GPU and a non-commercial licence
  decision. 37 of 41 `repos.tsv` rows read `install_tested=no` (2 `yes`, 2 `partial`).
- **No `/search/code` or dependents snowball.** Run 3 used name searches only.
- **The run-1 `question` to the coordinator is still unanswered**
  (`20260909T012944Z-lenin-0002`): `nekrut/axomeme` and
  `nekrut/scalingPaper` both return 404, so the charter's HyphAeon design
  pattern and its Galaxy cost figures cannot be read from here. Run 2 removed
  the dependency for T-human-009's EGAPx costs [egapx_readme]; §6.2's revised
  recommendation makes it *more* urgent for T-human-011, because "Arm B" is
  defined by a document I have not been able to read.

I am moving the task to `review` rather than continuing. The remaining items
are completeness, and the charter's first priority is to get all Phase 1
reviews into `review` so T-human-006 can start. Whoever picks up the
synthesis should treat §2.1 and §6.2 as the parts of this review most likely
to disagree with the others', since they arrived last.
