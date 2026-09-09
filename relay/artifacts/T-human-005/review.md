# Independent review of eukaryotic gene prediction: literature and software (T-human-005, marx)

Status: draft, tick 1 of several (2026-09-09). Sections 1 to 3 are populated
from verified metadata; sections 4 to 6 are outlines to be filled in the
next ticks. Every DOI below was resolved through the Europe PMC REST API on
2026-09-09; `refs.bib` holds the full records, `repos.tsv` the repository
survey. This review was written blind: I have not opened any other agent's
`relay/artifacts/T-human-00N/` directory.

## 1. Search log

Runner constraints that shaped the search: this agent runs in a cloud
container. Europe PMC (search and abstracts), bioRxiv API and Crossref were
reachable. OpenAlex answered HTTP 429 with a zero daily budget and Semantic
Scholar answered HTTP 429 on every call, so neither was used. The GitHub
REST API and github.com HTML were blocked by the container's proxy; the
repository survey therefore used anonymous `git clone --filter=blob:none
--shallow-since=2025-09-09`, which yields last-commit dates, commits in the
last 12 months, languages, licences and HEAD hashes but not stars or open
issue counts (recorded as `n/a` in `repos.tsv`; to be filled by an agent
with API access, or by me if access appears).

| # | date | source | query (Europe PMC syntax) | hits | kept |
|---|------|--------|---------------------------|------|------|
| 1 | 2026-09-09 | Europe PMC | 46 title-anchored queries, one per required method (GENSCAN, AUGUSTUS 2003/2006/2016, BRAKER1/2/3, GeneMark-ES/ET/EP+/ETP, SNAP, Gnomon, EGAPx, MAKER/MAKER2, TWINSCAN, N-SCAN, CONTRAST, Helixer, Tiberius, KA/KS, GeneID, GlimmerHMM, Exonerate, GeMoMa, miniprot, Spaln, SpliceAI, Pangolin, foundation models, funannotate, EviAnn, GALBA, Liftoff, TOGA, phastCons, phyloP, EGASP, nGASP, benchmarks, Mikado) | 1 to 386 per query | 44 primary papers |
| 2 | 2026-09-09 | Europe PMC | 41 gap queries: TWINSCAN 2001, N-SCAN 2006, CONTRAST 2007, PhyloCSF, RNAcode, CSF, CESAR, TOGA 2023, SegmentNT, PlantCaduceus, AgroNT, SpliceBERT, OpenSpliceAI, Spliceator, `(gene finding OR gene prediction OR genome annotation) AND (deep learning OR transformer OR language model) AND eukaryot* AND PUB_YEAR:[2023 TO 2026]`, author queries for Stanke, Hoff, Borodovsky, Lomsadze, Denton 2024 to 2026, ProtHint, Splign, intron length, GC isochores, BUSCO, compleasm, OMArk, Ensembl, Cactus, multiz, annotation-quality reviews | 0 to 277 | 26 more |
| 3 | 2026-09-09 | Europe PMC `resultType=core` | abstract and licence retrieval for the 54 DOIs in `refs.bib` | 54/54 resolved | 54 |
| 4 | 2026-09-09 | git (anonymous, shallow) | 36 repositories starting from Gaius-Augustus/Tiberius and ncbi/egapx and following citations and READMEs outward | 30 cloned, 6 non-existent slugs (ncbi/gnomon, Jstacs/GeMoMa, LarsGab/Tiberius, and three guesses) | 30 rows in `repos.tsv` |

Not yet done: bioRxiv full-text for the 2026 preprints (Tiberius clades,
Vipsania, OrionGeno, PlantGeneAnn); README install tests in a fresh venv;
Gnomon documentation from NCBI (no paper and no public repository); GeMoMa
(hosted at jstacs.de, not GitHub).

Observation from the search itself: Europe PMC's title index finds zero
2023 to 2026 papers that combine "gene prediction" with "transformer" or
"language model" in the title, and only six with those terms in the
abstract. The 2025 to 2026 wave (ANNEVO, OrionGeno, PlantGeneAnn, Vipsania,
Tiberius multi-clade, Helixer in Nature Methods, EviAnn) is therefore recent
enough that citation counts are not yet informative; the field moved in the
last 18 months.

## 2. Publications table

Approach classes: GHMM = ab initio (generalized) hidden Markov model;
COMP = comparative or alignment-based; EVID = evidence-based pipeline;
DL = deep learning; HYB = hybrid. "Accuracy" quotes the paper's own headline
figure with its benchmark; figures I have not yet read in full text are
marked (abstract). Runtime and hardware are from the paper or the README at
the HEAD hash in `repos.tsv`.

| Method | Year | DOI | Class | Inputs | Clades trained / evaluated | Reported accuracy (benchmark, metric) | Runtime, hardware | Cross-species evidence | Code |
|---|---|---|---|---|---|---|---|---|---|
| GENSCAN | 1997 | 10.1006/jmbi.1997.0951 | GHMM | genome only | human, vertebrates; isochore-specific parameters | classic Burset-Guigo set: exon Sn/Sp about 0.78/0.81 (to verify against full text; abstract only via Europe PMC) | minutes per Mb, CPU | none; trained on human, applied widely with known over-prediction in other clades | binaries only |
| TWINSCAN | 2001 | 10.1093/bioinformatics/17.suppl_1.s140 | COMP (pairwise informant) | genome + one informant genome | mouse with human informant | "dramatic improvement in exact gene Sn/Sp over GENSCAN" (abstract) | CPU, minutes per Mb | requires informant at suitable distance; distance sensitivity studied in 10.1101/sqb.2003.68.125 | historical |
| N-SCAN (TWINSCAN 3.0) | 2006 | 10.1089/cmb.2006.13.379 | COMP (phylogenetic, multi-informant) | genome + multiple alignment + tree | human, D. melanogaster | "exceeds all previously published whole-genome de novo predictors" in human and fly (abstract) | CPU | first to inject a phylogenetic model of the informants; the closest ancestor of what we want to build | historical |
| CONTRAST | 2007 | 10.1186/gb-2007-8-12-r269 | COMP (discriminative, phylogeny-free) | genome + multiple alignment | human (ENCODE regions) | exact CDS structure for 65% more human genes than N-SCAN; 46% fewer missed exons (abstract; OA) | CPU | uses alignments without a tree: SVM boundary classifiers + CRF; tree-free is the design choice we should argue against | historical |
| KA/KS test | 2002 | 10.1101/gr.200901 | COMP (single signal) | human-mouse aligned windows | human/mouse | about 9.5% FN and 2 to 3% FP for exons (charter, citing PMC155263); abstract: "FN lower than most, FP lower than all current methods" | trivial | pairwise only; needs suitable divergence and exon length | none |
| PhyloCSF | 2011 | 10.1093/bioinformatics/btr209 | COMP (codon substitution frequencies) | multiple alignment + phylogenetic codon model | 12 flies, 29 mammals | strong coding/non-coding discrimination; used to find hundreds of new human exons (10.1101/gr.246462.118) | CPU heavy per window; PhyloCSF++ (10.1093/bioinformatics/btab756) is much faster | needs clade-specific empirical codon models | mlin/PhyloCSF (AGPL), cpockrandt/PhyloCSFpp |
| RNAcode | 2011 | 10.1261/rna.2536111 | COMP (single signal) | multiple alignment | broad | robust coding detection in alignments (abstract) | CPU | designed to be species-independent | ViennaRNA site |
| AUGUSTUS | 2003, 2006 | 10.1093/bioinformatics/btg1080; 10.1093/nar/gkl200 | GHMM (+hints) | genome; optional hints from ESTs, proteins, RNA-seq | human, fly, Arabidopsis, many via training | EGASP 2006 best ab initio and hint-based results (10.1186/gb-2006-7-s1-s11); the 2003 paper introduced an explicit intron length submodel with GC-dependent parameters | CPU, roughly linear; hours per mammalian genome per core | needs per-species training; comparative mode (cgp, 10.1093/bioinformatics/btw494) predicts jointly across a whole-genome alignment | Gaius-Augustus/Augustus |
| SNAP | 2004 | 10.1186/1471-2105-5-59 | GHMM | genome | tested transfer between C. elegans, D. melanogaster, A. thaliana, O. sativa | accuracy collapses when parameters are transferred across species; motivates species-specific training (abstract; OA) | fast, CPU | direct negative evidence on cross-species transfer of HMM parameters | KorfLab/SNAP |
| GlimmerHMM / TigrScan | 2004 | 10.1093/bioinformatics/bth315 | GHMM | genome | Arabidopsis, human (ENCODE via JIGSAW paper 10.1186/gb-2006-7-s1-s9) | comparable to GENSCAN on its test sets (abstract) | fast, CPU | per-species training | JHU site |
| GeneID | 2000, 2007 | 10.1101/gr.10.4.511; 10.1002/0471250953.bi0403s18 | GHMM-like (weighted signals) | genome | Drosophila, many | competitive in the Drosophila GASP | fast, CPU | per-species parameter files | CRG site |
| GeneMark-ES | 2005 | 10.1093/nar/gki937 | GHMM, self-training | genome only | fungi, plants, animals | unsupervised training reaches supervised accuracy in tested species (abstract; OA) | CPU, hours | self-training is the classic answer to cross-species use; struggles on large intron-rich genomes | Borodovsky lab tarball (licence-restricted) |
| GeneMark-ET | 2014 | 10.1093/nar/gku557 | HYB (GHMM + RNA-seq introns) | genome + mapped RNA-seq | fly, Arabidopsis, C. elegans | improved over ES (abstract; OA) | CPU | as ES | same |
| GeneMark-EP+ | 2020 | 10.1093/nargab/lqaa026 | HYB (GHMM + protein hints via ProtHint) | genome + protein database | six eukaryotes at several relatedness levels | accuracy rises with proteome relatedness; details from full text pending (OA) | CPU | explicit study of evidence distance | gatech-genemark/ProtHint |
| GeneMark-ETP | 2024 | 10.1101/gr.278373.123 | EVID/HYB | genome + RNA-seq + protein db | large plant and animal genomes | outperforms MAKER2 and TSEBRA; margin grows with genome size (abstract; OA) | CPU, multi-core, hours to days | claims robustness on large genomes; still per-genome self-training | gatech-genemark/GeneMark-ETP (CC BY-NC-SA) |
| BRAKER1 / BRAKER2 | 2016 / 2021 | 10.1093/bioinformatics/btv661; 10.1093/nargab/lqaa108 | EVID (GeneMark + AUGUSTUS) | RNA-seq (B1) or proteins (B2) | 12 species at several protein-distance levels (B2) | B2: accuracy strongly dependent on protein distance (full text pending) | CPU, hours to days | per-genome training; evidence distance is the driver | Gaius-Augustus/BRAKER |
| BRAKER3 | 2024 | 10.1101/gr.278090.123 | EVID (GeneMark-ETP + AUGUSTUS + TSEBRA) | genome + RNA-seq + protein db | 11 species | best evidence-based accuracy among tested; transcript F1 improvements over B1/B2 (exact numbers pending from OA full text) | CPU, many hours | benchmarked at an assumed proteome relatedness level | Gaius-Augustus/BRAKER |
| GALBA | 2023 | 10.1186/s12859-023-05449-z | EVID (miniprot + AUGUSTUS) | genome + close reference proteins | vertebrates | competitive with BRAKER2 when close proteins exist (abstract; OA) | CPU | designed for the close-reference case | Gaius-Augustus/GALBA |
| MAKER / MAKER2 | 2008 / 2011 | 10.1101/gr.6743907; 10.1186/1471-2105-12-491 | EVID (combiner) | genome + ESTs/RNA-seq + proteins + ab initio predictors | emerging model organisms | annotation edit distance framework (full text pending) | CPU, MPI, days | combiner; accuracy bounded by its predictors | Yandell-Lab/maker (no commits in 12 months) |
| Gnomon / EGAPx | (no paper) | RefSeq: 10.1093/nar/gkae1038; NCBI docs | EVID | genome + RNA-seq + proteins | vertebrates, arthropods, plants; fungi, protists, nematodes out of scope (README at f9a7392) | no published benchmark; RefSeq curation paper describes the pipeline | README: 32 CPUs, 256 GB RAM; charter: about 416,000 CPU-hours for 1,409 Galaxy jobs, 45% failure, 24% CPU efficiency (relay/TASK.md) | clade scope explicitly limited | ncbi/egapx |
| EviAnn | 2026 | 10.1038/s41592-026-03156-0 | EVID (alignment-first, minimal ab initio) | genome + transcripts and/or proteins | eukaryotes incl. mammals | outperforms BRAKER3, MAKER2, FINDER on identical inputs; mammal in under 1 h on one multicore server (abstract; not OA) | CPU, under 1 h claimed | evidence-only; no cross-species model | alekseyzimin/EviAnn |
| TOGA + CESAR | 2023 / 2017 | 10.1126/science.abn3107; 10.1093/bioinformatics/btx527 | COMP (annotation transfer) | genome alignment chains + reference annotation | 488 placental mammals, 501 birds | orthology-based projection; accuracy tied to reference | CPU cluster | transfer, not prediction; fails on lineage-specific genes | hillerlab/TOGA |
| Liftoff | 2021 | 10.1093/bioinformatics/btaa1016 | COMP (lift-over) | reference annotation + minimap2 | same or close species | high within-species fidelity | minutes | not cross-clade | agshumate/Liftoff |
| Helixer | 2021, 2026 | 10.1093/bioinformatics/btaa1044; 10.1038/s41592-025-02939-1 | DL (CNN/biLSTM base labelling) + HMM post-processing | genome only | one vertebrate model over 186 animal genomes, one land-plant model over 51; 2026: fungi, plants, vertebrates, invertebrates | 2021: less sensitive to genome length than AUGUSTUS; 2026: "on par with or exceeding current tools" (abstract; OA) | GPU for training; inference GPU or CPU | explicitly cross-species within a kingdom-scale model; 2026 claims no retraining needed | weberlab-hhu/Helixer |
| Tiberius | 2024 | 10.1093/bioinformatics/btae685 | DL (CNN + biLSTM + differentiable HMM, end-to-end) | genome (+ optional softmasking, ClaMSA evolutionary input) | trained on mammals; evaluated on human and two others | human gene-level F1 62% vs 21% for next best ab initio; matches BRAKER3; human genome under 2 h (abstract; OA) | GPU >= 8 GB recommended (README) | mammal-only models in 2024 | Gaius-Augustus/Tiberius |
| Tiberius multi-clade | 2026 | 10.64898/2026.04.24.720536 | DL | genome | Mesangiospermae, Fungi, Vertebrata, Insecta, Chlorophyta, Bacillariophyta; 33-species benchmark | gene-level F1 12 to 37 points over Helixer, 10 to 22 over ANNEVO; about 80x faster than BRAKER3 on GPU; backend 31% faster (abstract) | GPU | separate model per lineage, applicable to "92% of available assemblies" | same |
| ANNEVO | 2026 | 10.1038/s41592-026-03036-7 | DL (mixture-of-experts genomic language model) | genome | 566 species benchmark | "substantially outperforms existing ab initio methods, comparable to pipelines" (abstract; not OA) | unknown; likely GPU | claims a single evolution-aware model | code link pending |
| OrionGeno | 2026 | 10.64898/2026.04.26.720859 | DL (phylogeny-aware, long-range, joint genes + repeats + UTRs) | genome | "diverse eukaryotic lineages"; ran on >5,300 unannotated NCBI genomes | outperforms state of the art at exon, gene, protein levels (abstract) | unknown | phylogenetic context as input: directly relevant to our design | pending |
| Vipsania | 2026 | 10.64898/2026.08.26.747235 | DL, unsupervised (differentiable HMM inside a deep model) | unannotated genome only | pretrained for "virtually all eukaryotes", fine-tunes without labels | more accurate than supervised methods across most clades; no drop on distant genomes; handles non-standard codes (abstract) | unknown | strongest cross-clade claim in the cohort | Gaius-Augustus/Vipsania |
| PlantGeneAnn | 2026 | 10.64898/2026.06.25.733695 | DL (strand-specific plant foundation model) | genome | fine-tuned on 9 plants; 13-species benchmark | beats four baselines at five levels; data quality beats token volume (abstract) | GPU | plants only | pending |
| Sensor-NN ab initio | 2023 | 10.1093/bioadv/vbad105 | DL (small NN over hand-built sensors) | genome | phylogenetically distant train/test | nucleotide-level accuracy above existing ab initio methods with far less training data (abstract; OA) | CPU | explicit cross-clade design; small model | pending |
| SpliceAI / Pangolin / OpenSpliceAI / Spliceator | 2019 / 2022 / 2025 / 2021 | 10.1016/j.cell.2018.12.015; 10.1186/s13059-022-02664-4; 10.7554/elife.107454; 10.1186/s12859-021-04471-3 | DL (splice-site CNNs) | pre-mRNA sequence | human; four mammals; retrainable; multi-species | splice-site level only; not gene structure | GPU for training | Spliceator and OpenSpliceAI address species transfer | see repos.tsv; SpliceAI licence is PolyForm Strict |
| Nucleotide Transformer / SegmentNT, Evo 2, Caduceus | 2025 / 2026 | 10.1038/s41592-024-02523-z; 10.1038/s41586-026-10176-5 | DL (foundation models) | genome | broad | element-level segmentation, not full gene models; Evo 2 at 7B and 40B parameters | multi-GPU | outside the charter's compute budget; useful as baselines and as evidence about what scale buys | see repos.tsv |
| G3PO benchmark | 2020 | 10.1186/s12864-020-6707-9 | benchmark | curated genes from 147 organisms | broad | independent comparison of ab initio tools; identifies systematic weaknesses (OA) | n/a | directly reusable for our benchmark design | supplementary data |
| EGASP / nGASP | 2006 / 2008 | 10.1186/gb-2006-7-s1-s2; 10.1186/1471-2105-9-549 | benchmark | ENCODE human; C. elegans | human; nematode | community assessments; historical baselines | n/a | n/a | n/a |
| fitild | 2018 | 10.1093/bioinformatics/bty353 | analysis | annotations of >1,000 genomes | broad | intron length distributions carry up to 30% of intron-recognition information in some species; modelled as 1 to 3 Frechet components (OA) | n/a | quantifies exactly the clade-dependence the charter worries about | ogotoh/fitild |

Still to add rows for: Exonerate, Spaln, Splign/ProSplign, GeMoMa, miniprot
(aligners the pipelines depend on), EVidenceModeler, funannotate, Mikado,
phastCons/phyloP, OMArk/BUSCO as evaluation tools, and the 2025 "200 insect
genomes with BRAKER" and 2026 "301 Drosophilidae comparative annotation"
papers as evidence on scaling costs.

## 3. Repository inventory

See `repos.tsv` (30 repositories, HEAD hashes as of 2026-09-09). Summary
of what the activity data say:

- Alive and moving in 2026: Tiberius (174 commits in 12 months), Helixer
  (37), Vipsania (11, new), EviAnn (99), funannotate (279), cactus (480),
  OpenSpliceAI (29), egapx (25), CESAR2.0 (14), BRAKER (11).
- Maintenance or dormant: Augustus (2), miniprot (2), SNAP (0, last commit
  2022), MAKER (0, 2024), GeneMark-ETP (0, 2025), GeneMark-EP (empty tree),
  ProtHint (0, 2023), Liftoff (0, 2023), PhyloCSF (0, 2023), PhyloCSF++
  (0, 2022), EVidenceModeler (0, 2024), Pangolin (0, 2023).
- Licence hazards for a derived tool: GeneMark family (CC BY-NC-SA or
  custom non-commercial), SpliceAI (PolyForm Strict), Nucleotide
  Transformer (CC BY-NC-SA). MIT/BSD/Apache/GPL elsewhere.
- Not on GitHub: Gnomon (NCBI internal), GeMoMa (jstacs.de), GeneID,
  GlimmerHMM. Not found yet: ANNEVO, OrionGeno, PlantGeneAnn repositories
  (need the preprint full text for links).

Install tests have not been run yet; the plan is a fresh `python3 -m venv`
per repository next tick, CPU only, recording success, wall time and the
first failure.

## 4. Data sources noted in passing (outline)

- Whole-genome alignments: UCSC multiz (10.1101/gr.1933104) and Cactus/HAL
  (10.1038/s41586-020-2871-y); Zoonomia 241-mammal alignment; the TOGA
  paper's 488-mammal and 501-bird chains.
- Conservation tracks: phastCons (10.1101/gr.3715005), phyloP
  (10.1101/gr.097857.109), PhyloCSF tracks (10.1101/gr.246462.118).
- Reference annotations used as training or truth: RefSeq
  (10.1093/nar/gkae1038), Ensembl (10.1093/database/baw093), G3PO curated
  genes (10.1186/s12864-020-6707-9), the Tiberius 33-species and ANNEVO
  566-species panels (to extract species lists from full text).
- Expression: VARUS-selected SRA RNA-seq as in BRAKER3 and the 200-insect
  study (10.1101/2025.04.17.649312).
- Quality tools reusable as metrics: BUSCO (10.1093/molbev/msab199),
  compleasm (10.1093/bioinformatics/btad595), OMArk
  (10.1038/s41587-024-02147-w).

## 5. Failure modes of current tools (outline, citations to be expanded)

- Cross-clade transfer of HMM parameters fails (SNAP, 10.1186/1471-2105-5-59);
  supervised deep models also drop on distant genomes (stated as the
  motivation of Vipsania, 10.64898/2026.08.26.747235) and Tiberius ships one
  model per lineage rather than one model (10.64898/2026.04.24.720536).
- Intron length: distributions vary enough across species to be a
  species fingerprint and carry up to 30% of intron-recognition information
  (fitild, 10.1093/bioinformatics/bty353); AUGUSTUS needed an explicit
  length submodel (10.1093/bioinformatics/btg1080).
- GC isochores: GENSCAN and AUGUSTUS use GC-binned parameters; to document
  the accuracy gradient with GC from full text.
- Evidence distance: BRAKER2 and GeneMark-EP+ accuracy tracks proteome
  relatedness (10.1093/nargab/lqaa108, 10.1093/nargab/lqaa026).
- Resource use: EGAPx 32 CPU / 256 GB and its Galaxy failure statistics
  (README f9a7392; relay/TASK.md); Tiberius needs an 8 GB GPU (README
  e73844b); BRAKER3 is about 80x slower than Tiberius on GPU
  (10.64898/2026.04.24.720536).
- Annotation quality control: OMArk finds erroneous gene inference at
  scale (10.1038/s41587-024-02155-w).

## 6. Opinion (to be written after sections 4 and 5; one page)

Provisional direction, to be argued properly next tick: the 2006 N-SCAN
design (multiple alignment plus an explicit phylogenetic model) was right
and was abandoned for lack of alignments and compute, not because it was
wrong. The 2026 cohort is converging back on it from the deep-learning
side (OrionGeno "phylogeny-aware", Tiberius "ClaMSA" input, ANNEVO
"evolutionary relationships"), but each learns the phylogeny implicitly
from sequence. The gap to fill is a model that receives alignment, tree and
codon geometry as coordinates, in the HyphAeon style described in
relay/TASK.md, and only has to learn a low-dimensional decision surface.
The biggest risk is alignment availability and quality in non-coding
regions for clades outside mammals and flies.
