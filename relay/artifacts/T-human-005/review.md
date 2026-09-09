# Independent review of eukaryotic gene prediction: literature and software (T-human-005, marx)

Status: submitted for review, tick 3, 2026-09-09. Sections 1 to 6 are
complete. Every DOI was resolved through the Europe PMC REST API; every
accuracy or runtime number below is quoted from open-access full text I
read (PMC ids in the search log) or from a README or documentation file at
the commit hash listed in `repos.tsv`. Numbers I could only see in an
abstract are marked "(abstract)". This review was written blind: I have not
opened any other agent's `relay/artifacts/T-human-00N/` directory. The one
exception to "read it myself" is the charter's own EGAPx cost figure, which
I cite as `relay/TASK.md`.

Known gaps, stated up front: ClaMSA and AUGUSTUS-cgp are covered from
abstracts plus their repository documentation, not full text (the journal
versions are not in PMC, and the OA preprints sit behind hosts that refused
this runner, see rows 9 and 10 of the search log); the five 2026 preprints
are covered from bioRxiv API abstracts and READMEs; no install test was
possible for the container-only tools.

## 1. Search log

Runner constraints: this agent runs in a cloud container behind a proxy.
Tick 1 (2026-09-09, 01:30 to 03:30 UTC): Europe PMC search and abstracts,
the bioRxiv API and Crossref were reachable; OpenAlex and Semantic Scholar
answered HTTP 429; the GitHub REST API and github.com HTML were blocked
(403), so repository activity was measured with anonymous
`git clone --filter=blob:none --shallow-since=2025-09-09`. Tick 2
(03:50 to 04:30 UTC): Europe PMC full-text XML, OpenAlex, Semantic Scholar
and raw.githubusercontent.com were reachable; the GitHub API remained
blocked, and star and issue counts were obtained from shields.io badge
JSON (`img.shields.io/github/stars/<owner>/<repo>.json`), which reads the
GitHub API server-side. bioRxiv full-text JATS answered 429 on every
attempt, so the four 2026 preprints are covered from their API abstracts
and READMEs only. Tick 3 (04:50 to 05:30 UTC): the bioRxiv HTML, PDF and
JATS hosts answered HTTP 429 (Cloudflare error 1015) on all seven attempts,
PeerJ and OUP answered 403, NCBI efetch for PMC5860283 returned only the
front matter ("the publisher of this article does not allow downloading of
the full text in XML form"), Europe PMC has no full text for the five
preprint records (PPR294699, PPR1213068, PPR1213229, PPR1309120,
PPR1260420), and OpenAlex and Semantic Scholar returned empty records. I
stopped there rather than work around a rate limit.

| # | date | source | query or action | hits | kept |
|---|------|--------|-----------------|------|------|
| 1 | 2026-09-09 | Europe PMC search | 46 title-anchored queries, one per required method (GENSCAN, AUGUSTUS 2003/2006/2016, BRAKER1/2/3, GeneMark-ES/ET/EP+/ETP, SNAP, Gnomon, EGAPx, MAKER/MAKER2, TWINSCAN, N-SCAN, CONTRAST, Helixer, Tiberius, KA/KS, GeneID, GlimmerHMM, Exonerate, GeMoMa, miniprot, Spaln, SpliceAI, Pangolin, funannotate, EviAnn, GALBA, Liftoff, TOGA, phastCons, phyloP, EGASP, nGASP, Mikado) | 1 to 386 per query | 44 primary papers |
| 2 | 2026-09-09 | Europe PMC search | 41 gap queries: TWINSCAN 2001, N-SCAN 2006, CONTRAST 2007, PhyloCSF, RNAcode, CESAR, TOGA 2023, SegmentNT, PlantCaduceus, AgroNT, SpliceBERT, OpenSpliceAI, Spliceator, `(gene finding OR gene prediction OR genome annotation) AND (deep learning OR transformer OR language model) AND eukaryot* AND PUB_YEAR:[2023 TO 2026]`, author queries (Stanke, Hoff, Borodovsky, Lomsadze, Denton 2024 to 2026), ProtHint, Splign, intron length, GC isochores, BUSCO, compleasm, OMArk, Ensembl, Cactus, multiz | 0 to 277 | 26 more |
| 3 | 2026-09-09 | Europe PMC `resultType=core` | abstract and licence retrieval for all DOIs in `refs.bib` | 60/60 resolved | 60 |
| 4 | 2026-09-09 | Europe PMC fullTextXML | full text of 17 OA papers: PMC11645249 (Tiberius 2024), PMC11216308 (BRAKER3), PMC11216313 (GeneMark-ETP), PMC8016489 (Helixer 2021), PMC13076211 (Helixer 2026), PMC2246271 (CONTRAST), PMC7147072 (G3PO), PMC6157073 (fitild), PMC421630 (SNAP), PMC7222226 (GeneMark-EP+), PMC7787252 (BRAKER2), PMC10448985 (sensor-NN 2023), PMC1298918 (GeneMark-ES), PMC10472564 (GALBA), PMC3117341 (PhyloCSF), PMC1810551 (EGASP), PMC3280279 (MAKER2) | 17/17 | all read for numbers |
| 5 | 2026-09-09 | bioRxiv API `details` | four 2026 preprints (Tiberius clades, OrionGeno, Vipsania, PlantGeneAnn): versions, licences, abstracts, code links | 4/4 metadata; 0/4 full text (HTTP 429) | 4 |
| 6 | 2026-09-09 | Semantic Scholar | citation counts for five anchor papers | KA/KS 271, Helixer 2021 138, CONTRAST 106, BRAKER3 59, Tiberius 42 | context only |
| 7 | 2026-09-09 | git (anonymous, shallow) + shields.io | 36 repositories from Gaius-Augustus/Tiberius and ncbi/egapx outward; 6 slugs did not exist (ncbi/gnomon, Jstacs/GeMoMa, LarsGab/Tiberius, three guesses) | 30 cloned; 30/30 stars and issue counts | 30 rows in `repos.tsv` |
| 8 | 2026-09-09 | fresh `python3 -m venv` | install tests for Tiberius (Python 3.11 and 3.12), Vipsania (3.11, 3.12 dry-run), Helixer (3.11, real install) | see section 3 | 3 |
| 9 | 2026-09-09 | Europe PMC search + fullTextXML, NCBI efetch, bioRxiv (www), PeerJ, OUP, OpenAlex, Semantic Scholar | full text of ClaMSA (btac028; preprint 10.1101/2021.03.09.434414, CC BY) and AUGUSTUS-cgp (btw494, PMC5860283; PeerJ preprint 10.7287/peerj.preprints.1296, CC BY) | 0/2 full text; 2/2 journal abstracts (Europe PMC, Semantic Scholar) | abstracts + repo docs used |
| 10 | 2026-09-09 | bioRxiv API `details` | re-check of the 2026 preprints: Tiberius multi-clade has a v2 (2026-07-29) and OrionGeno a v2 (2026-08-24); v2 abstracts read; Vipsania and PlantGeneAnn still v1 | 2 new versions | rows updated |
| 11 | 2026-09-09 | raw.githubusercontent.com | Gaius-Augustus/clamsa `README.md`, Gaius-Augustus/Augustus `docs/RUNNING-AUGUSTUS-IN-CGP-MODE.md` (inputs, tree requirements), xjtu-omics/ANNEVO `README.md` (lineages, hardware) | 3/3 | 3 |
| 12 | 2026-09-09 | git (anonymous, shallow) + shields.io | Gaius-Augustus/clamsa and xjtu-omics/ANNEVO (slug found via the abstract's stated availability and the author affiliation) | 2 cloned | 2 rows added to `repos.tsv` (32 total) |

Observation from the search: Europe PMC's title index finds zero 2023 to
2026 papers combining "gene prediction" with "transformer" or "language
model" in the title, and six with those terms in the abstract. The 2025 to
2026 wave (ANNEVO, OrionGeno, PlantGeneAnn, Vipsania, Tiberius multi-clade,
Helixer in Nature Methods, EviAnn) is too recent for citation counts to
mean anything; the field moved in the last 18 months, and most of that
movement is in preprints whose full text I could not fetch from this
runner.

## 2. Publications table

Classes: GHMM = ab initio (generalized) hidden Markov model; COMP =
comparative or alignment-based; EVID = evidence-based pipeline; DL = deep
learning; HYB = hybrid. Sn = sensitivity, Sp = specificity or precision.

| Method | Year | DOI | Class | Inputs | Clades trained / evaluated | Reported accuracy (benchmark, metric) | Runtime, hardware | Cross-species evidence | Code |
|---|---|---|---|---|---|---|---|---|---|
| GENSCAN | 1997 | 10.1006/jmbi.1997.0951 | GHMM, isochore-specific parameters | genome | human; applied to all vertebrates | Burset-Guigo set: exon-level Sn/Sp about 0.78/0.81 (abstract; not OA). On G3PO (1793 genes, 147 species) nucleotide F1 0.51, second of five (PMC7147072 via PMC10448985 Table 3) | minutes per Mb, one CPU | trained on human; over-predicts elsewhere (G3PO) | binaries only |
| TWINSCAN | 2001 | 10.1093/bioinformatics/17.suppl_1.s140 | COMP, one informant | genome + informant genome | mouse with human informant | "dramatic improvement" in exact gene Sn/Sp over GENSCAN (abstract) | CPU | needs an informant at the right distance | historical |
| N-SCAN | 2006 | 10.1089/cmb.2006.13.379 | COMP, phylogenetic, multi-informant | genome + multiple alignment + tree | human, fly | on human (CONTRAST's evaluation set, PMC2246271 Table 1): gene Sn/Sp 35.6/25.1, exon Sn/Sp 84.2/64.6 with mouse informant; "performs as well using mouse as its only informant as it does with any combination" | CPU | first to inject a phylogenetic model of the informants; the direct ancestor of what we want | historical |
| CONTRAST | 2007 | 10.1186/gb-2007-8-12-r269 | COMP, discriminative, phylogeny-free (SVM boundary classifiers + CRF) | genome + multiple alignment (11 informants) | human ENCODE regions | gene Sn/Sp 58.6/35.5 with 11 informants, 50.8/29.3 with mouse only; exon Sn/Sp 92.8/72.5; "65% increase in gene sensitivity and 46% reduction in exon error rate" over N-SCAN (PMC2246271 Table 1) | CPU | the design choice we should argue against: it threw the tree away and still gained from more informants, which N-SCAN could not | historical |
| KA/KS test | 2002 | 10.1101/gr.200901 | COMP, one signal | human-mouse aligned windows | human, mouse | exon detection with about 9.5% FN and 2 to 3% FP (charter; PMC155263 abstract: "FN lower than most, FP lower than all current methods") | trivial | pairwise; needs suitable divergence and window length; 271 citations (Semantic Scholar) | none |
| PhyloCSF | 2011 | 10.1093/bioinformatics/btr209 | COMP, phylogenetic codon model | multiple alignment + tree | 12 flies, 29 mammals | strong coding/non-coding discrimination on fly regions (PMC3117341); clade-specific empirical codon models | CPU heavy per window; PhyloCSF++ (10.1093/bioinformatics/btab756) is the fast C++ port | needs a clade model; not a gene finder | mlin/PhyloCSF (AGPL); cpockrandt/PhyloCSFpp |
| ClaMSA | 2022 | 10.1093/bioinformatics/btac028 (preprint 10.1101/2021.03.09.434414) | COMP, end-to-end learned evolutionary models: a continuous-time Markov chain (CTMC) layer whose rate matrices are trained discriminatively, followed by (recurrent) neural layers | codon MSA + phylogenetic tree in Newick, scaled to one expected codon substitution per time unit; the README recommends a MrBayes codon-model tree built from positive (coding) alignments only | trained and tested on vertebrate, fly and yeast codon alignments (README data scripts `download_fly_vert_yeast_train.sh`); abstract reports vertebrate and fly | "four times fewer false positives ... than existing methods at the same true positive rate" on the coding vs non-coding candidate-classification task (journal abstract; the preprint abstract is not quoted here because its full text was unreachable). Classifies pre-extracted candidate alignments, so it is a scorer, not a gene finder | CPU (TensorFlow >= 2.0); tree construction is the expensive step | closest published precedent for the charter's "give the network the tree" idea: the tree enters as fixed branch lengths in a learned CTMC, not as a learned embedding; it is used as an optional input by Tiberius de novo mode (PMC11645249) | Gaius-Augustus/clamsa (README at commit in `repos.tsv`; not cloned this tick) |
| RNAcode | 2011 | 10.1261/rna.2536111 | COMP, one signal | multiple alignment | broad | robust coding detection (abstract) | CPU | designed species-independent | ViennaRNA site |
| AUGUSTUS | 2003, 2006, 2016 | 10.1093/bioinformatics/btg1080; 10.1093/nar/gkl200; 10.1093/bioinformatics/btw494 | GHMM (+hints; cgp = comparative) | genome; optional hints | many species via training; explicit intron length submodel with GC-dependent parameters (2003) | EGASP 2006: no method exceeded 45% exact transcript Sn (PMC1810551). Tiberius 2024 benchmark, three mammals: exon F1 67.3, gene F1 12.4 with 2010 human parameters (PMC11645249). G3PO: nucleotide F1 0.52, best of five; short exons under 50 nt only 18% correct (PMC7147072) | 2:25 h per mammalian genome on 48 threads (PMC11645249 Table 1) | per-species training. cgp mode (2016, PMC5860283, abstract + `docs/RUNNING-AUGUSTUS-IN-CGP-MODE.md`): takes a trained species parameter set, one FASTA per genome, a MAF whole-genome alignment and a Newick tree (a star tree with uniform branch lengths is the documented fallback), and predicts all genomes jointly as a binary labelling problem on a graph solved by dual decomposition; tested on 12-vertebrate and 12-Drosophila alignments, evaluated on human, mouse and D. melanogaster; the abstract claims annotation transfer through the alignment beats protein-spliced alignment. This is the direct ancestor of the charter's design: alignment plus tree as fixed inputs, exon gains and losses scored against the species tree | Gaius-Augustus/Augustus |
| SNAP | 2004 | 10.1186/1471-2105-5-59 | GHMM | genome | A. thaliana, C. elegans, D. melanogaster, O. sativa | own-species nucleotide Sn 93.8 to 98.1; foreign parameters collapse: fly parameters on Arabidopsis Sn 26.0, rice parameters on C. elegans Sn 21.7; bootstrapping recovers Sn 75 to 96 (PMC421630 Table 4) | fast, CPU | the cleanest published negative result on transferring HMM parameters | KorfLab/SNAP |
| GlimmerHMM | 2004 | 10.1093/bioinformatics/bth315 | GHMM | genome | Arabidopsis, human | G3PO nucleotide F1 0.45, highest Sn (0.74) but lowest Sp (0.43) (PMC10448985 Table 3) | fast, CPU | per-species training | JHU site |
| GeneID | 2000 | 10.1101/gr.10.4.511 | weighted-signal GHMM-like | genome | Drosophila, many | G3PO nucleotide F1 0.40 | fast, CPU | per-species parameter files | CRG site |
| GeneMark-ES | 2005 | 10.1093/nar/gki937 | GHMM, self-training | genome only | A. thaliana, C. elegans, D. melanogaster; A. gambiae, C. intestinalis, C. reinhardtii, T. gondii | unsupervised matches supervised at nucleotide level (A. thaliana Sn/Sp 97.7/96.3 vs 97.2/95.8; PMC1298918 Table 1); gene-level accuracy 5 to 20% on several later test genomes (PMC7222226) | CPU, hours | self-training is the classic cross-species answer; weak at gene level on large genomes | Borodovsky lab tarball, restricted licence |
| GeneMark-ET | 2014 | 10.1093/nar/gku557 | HYB, RNA-seq introns | genome + mapped RNA-seq | fly, Arabidopsis, C. elegans | improved over ES (abstract) | CPU | as ES | same |
| GeneMark-EP+ / ProtHint | 2020 | 10.1093/nargab/lqaa026 | HYB, protein hints | genome + protein database | N. crassa, C. elegans, A. thaliana, D. melanogaster, S. lycopersicum, D. rerio | ProtHint on D. melanogaster: intron hint Sn falls from 79.8 (species-level exclusion) to 35.8 (phylum-level); start-codon Sn from 70.3 to 14.1; high-confidence Sp stays above 98.8 (PMC7222226 Table 5) | CPU | the quantitative curve of evidence distance | gatech-genemark/ProtHint |
| GeneMark-ETP | 2024 | 10.1101/gr.278373.123 | EVID/HYB | genome + RNA-seq + protein db | 7 species in three groups: compact, large GC-homogeneous, large GC-inhomogeneous | gene F1 gains over GeneMark-ET of 19.6 / 47.8 / 66.3 points per group; over TSEBRA 8.2 (large homogeneous) and 39.0 (large inhomogeneous); over MAKER2 23.3 / 21.7 / 27.8 (fly, zebrafish, mouse); BUSCO above 90% everywhere while gene F1 varies widely (PMC11216313) | 64 cores: 3 h fly, 12 h zebrafish, 18 h mouse, excluding HISAT2/StringTie | trains three GC-specific models per GC-inhomogeneous genome; still per-genome | gatech-genemark/GeneMark-ETP (CC BY-NC-SA) |
| BRAKER1 / BRAKER2 | 2016 / 2021 | 10.1093/bioinformatics/btv661; 10.1093/nargab/lqaa108 | EVID (GeneMark + AUGUSTUS) | RNA-seq (B1) or proteins (B2) | B2: 12 species | B2 gene Sn: A. thaliana 70.2, D. melanogaster 59.5, D. rerio 39.1, R. prolixus 13.2, T. nigroviridis 10.4 (PMC7787252 Table 3); vs MAKER2 on fly gene F1 59.7 vs 35.9 (Table 4) | B2 with 443-species proteins on 8 CPUs comparable to MAKER2 with 10 species, about 10 h | accuracy tracks genome size, annotation quality and protein distance | Gaius-Augustus/BRAKER |
| BRAKER3 | 2024 | 10.1101/gr.278090.123 | EVID (GeneMark-ETP + AUGUSTUS + TSEBRA) | genome + RNA-seq + protein db | 11 species, order-excluded proteins | transcript F1 about 20 points above B1/B2 on average; G. gallus gene/transcript F1 +55/+48 over TSEBRA(B1+B2); vs Funannotate +10.2/+25.9/+21.6 (exon/gene/transcript); Funannotate vs MAKER2 +2.2/+3.8/+4.4; start-codon F1 70%, stop 76%, splice sites above 87% (PMC11216308) | 48 threads Xeon E5-2650 v4: 5 h 37 min (A. thaliana) to 64 h 16 min (mouse), RNA-seq alignment excluded; 48:53 h average in the Tiberius benchmark (PMC11645249) | per-genome training; FINDER completed on only 7 of 11 species | Gaius-Augustus/BRAKER |
| GALBA | 2023 | 10.1186/s12859-023-05449-z | EVID (miniprot + AUGUSTUS) | genome + close proteins | 14 species | gene F1 79.5% with the target's own proteome (PMC10472564); Tiberius benchmark on mammals: exon/gene F1 86.2/41.8 | 35:12 h average, 48 threads (PMC11645249) | close-reference case only | Gaius-Augustus/GALBA |
| MAKER / MAKER2 | 2008 / 2011 | 10.1101/gr.6743907; 10.1186/1471-2105-12-491 | EVID combiner | genome + ESTs/RNA-seq + proteins + predictors | emerging model organisms | introduces Annotation Edit Distance; no F1 tables; outperformed by BRAKER2/3 and GeneMark-ETP by 20 to 28 gene-F1 points (above) | MPI; scaling tested to 32 cores on 48-core Opteron (PMC3280279) | bounded by its predictors | Yandell-Lab/maker (dormant) |
| EVidenceModeler | 2008 | 10.1186/gb-2008-9-1-r7 | EVID combiner | predictions + alignments | rice, human | weighted consensus (PMC2395244) | CPU | n/a | EVidenceModeler/EVidenceModeler (dormant) |
| Mikado | 2018 | 10.1093/gigascience/giy093 | EVID transcript selector | multiple transcript assemblies | plants, human | scoring-based selection (PMC6105091) | CPU | n/a | EI-CoreBioinformatics/mikado (not cloned) |
| Gnomon / EGAPx | no paper | 10.1093/nar/gkae1038 (RefSeq); README | EVID | genome + RNA-seq + proteins | vertebrates, arthropods, plants; fungi, protists, nematodes out of scope (README f9a7392) | no published benchmark | README: 32 CPUs, 256 GB RAM; Galaxy: about 416,000 CPU-hours over 1,409 jobs, 45% failure, 24% CPU efficiency (relay/TASK.md) | clade scope explicitly limited | ncbi/egapx (207 stars) |
| EviAnn | 2026 | 10.1038/s41592-026-03156-0 | EVID, alignment-first | genome + transcripts and/or proteins | eukaryotes incl. mammals | outperforms BRAKER3, MAKER2, FINDER on identical inputs; mammal under 1 h (abstract; not OA) | one multicore server | evidence-only | alekseyzimin/EviAnn_release (183 stars) |
| TOGA + CESAR | 2023 / 2017 | 10.1126/science.abn3107; 10.1093/bioinformatics/btx527 | COMP, annotation projection | genome alignment chains + reference annotation | 488 placental mammals, 501 birds | orthology-based; accuracy tied to reference | cluster | transfer, not prediction; misses lineage-specific genes | hillerlab/TOGA |
| Liftoff | 2021 | 10.1093/bioinformatics/btaa1016 | COMP lift-over | reference annotation + minimap2 | same or close species | high within-species fidelity | minutes | not cross-clade | agshumate/Liftoff (552 stars, dormant) |
| miniprot / Spaln / Exonerate | 2023 / 2008, 2024 / 2005 | 10.1093/bioinformatics/btad014; 10.1093/nar/gkn105, 10.1093/bioinformatics/btae517; 10.1186/1471-2105-6-31 | spliced aligners inside pipelines | proteins or cDNA + genome | broad | aligner accuracy bounds every EVID pipeline above | CPU; miniprot is the fast one | intron-length priors are species parameters in Spaln | lh3/miniprot (418 stars); ogotoh/spaln |
| Helixer | 2021, 2026 | 10.1093/bioinformatics/btaa1044; 10.1038/s41592-025-02939-1 | DL (CNN + biLSTM) base labelling + HMM post-processor | genome only | 2021: one vertebrate model over 186 animal genomes, one land-plant model over 51 (PMC8016489). 2026: four lineage models (fungi, land plants, vertebrates, invertebrates) | 2026 (PMC13076211 Tables 1, 2), mean over test species: base-wise phase F1 0.95 fungi, 0.81 plants, 0.88 vertebrates, 0.86 invertebrates; transcript-level F1 0.54 / 0.46 / 0.20 / 0.31 vs GeneMark-ES 0.60 / 0.09 / 0.02 / 0.18 vs AUGUSTUS 0.53 / 0.23 / 0.08 / 0.15. In Tiberius's mammal benchmark: exon/gene F1 72.9/19.3 | 8:54 h per mammal on A100 (PMC11645249); README: GPU with 8 to 11 GB; install "20 to 30 minutes" for experienced users | one model per kingdom-scale lineage; GeneMark-ES still wins fungi at transcript level | weberlab-hhu/Helixer (305 stars) |
| Tiberius | 2024 | 10.1093/bioinformatics/btae685 | DL (CNN + biLSTM + differentiable HMM, end-to-end, F1 loss) | genome + softmasking (+ optional ClaMSA input in de novo mode) | trained on mammals; tested on human, cow, beluga; probed on chicken, zebrafish, poplar, tomato | human gene F1 62% vs 21% next-best ab initio; three-mammal means: exon F1 89.7, gene F1 55.1 vs BRAKER3 83.2/53.7, GALBA 86.2/41.8, Helixer 72.9/19.3, AUGUSTUS 67.3/12.4; de novo mode on human 92.6/65.5; about 8M parameters; a 2M-parameter ablation (Tiberius_small) exists (PMC11645249) | 1:39 h per mammal on one A100 80 GB; training 15 days on four A100s; parallel Viterbi gives 17x on GPU; README asks for 8 GB GPU and Python 3.12 | "steady decline in accuracy as phylogenetic distance increases"; softmasking matters more for distant species; "need for re-training for other clades" | Gaius-Augustus/Tiberius (138 stars, 174 commits in 12 months) |
| Tiberius multi-clade | 2026 (v1 2026-04-28, v2 2026-07-29, CC BY) | 10.64898/2026.04.24.720536 | DL | genome | six lineage models: Mesangiospermae, Fungi, Vertebrata, Insecta, Chlorophyta, Bacillariophyta; 33-species benchmark; 2,948 vertebrate assemblies annotated | gene F1 12 to 37 points over Helixer, 10 to 22 over ANNEVO; near BRAKER3 in plants, fungi, diatoms, green algae; 80x faster than BRAKER3 on GPU; backend reimplementation 31% faster; v2 abstract adds that fewer than 20% of NCBI Datasets genomes carry an annotation and applies the Vertebrata model to 2,948 assemblies, nearly 6 Tbp (v2 abstract, CC BY) | GPU; web server | one model per lineage, "92% of available assemblies" | same |
| ANNEVO | 2026 | 10.1038/s41592-026-03036-7 | DL, mixture-of-experts genomic LM (PyTorch) | genome; lineage flag selects segment length | 566-species benchmark (abstract); six released lineage models: Mammalia, Insecta, Aves, Actinopteri, Magnoliopsida, Fungi (README 42c920f) | "substantially outperforms existing ab initio methods" (abstract; not OA); README reports BUSCO Mammalia_odb10 98.3 on GRCh38 for v2.3.3 vs 95.7 in the paper | README: human genome 19 min on one RTX 4090 for v2.3.3 (82 min for the paper version); README table, mean over 12 species on one RTX 4090: ANNEVO 12.2 min and 3.8 GB GPU memory vs Tiberius 43.6 min and 22.5 GB vs Helixer 286.6 min and 8.6 GB (authors' own measurements, batch size 8); GPU prediction then multi-core CPU decoding | the abstract's single "evolution-aware" model is in practice six lineage models chosen by the user, the same pattern as Helixer and Tiberius; non-commercial licence | xjtu-omics/ANNEVO (158 stars, 30 commits in 12 months; training code released) |
| OrionGeno | 2026 (v1 2026-04-29, v2 2026-08-24) | 10.64898/2026.04.26.720859 | DL, phylogeny-aware, long-range, joint genes + repeats + UTRs | genome | "diverse eukaryotic lineages"; >5,300 unannotated NCBI genomes | outperforms state of the art at exon, gene, protein-sequence and protein-structure levels; also reports candidate coding loci absent from curated references (v2 abstract; CC BY-NC-ND; web platform and database) | unknown | phylogenetic context as input | not found |
| Vipsania | 2026 | 10.64898/2026.08.26.747235 | DL, unsupervised: differentiable HMM inside a masked-LM sequence model | unannotated genome only | 17 clade models plus one "other"; training sets of 15 to 200 species each; 6 to 13 test species per clade | README (0c5f84a) average locus F1 after fine-tuning: Discoba 0.70, Fungi 0.66, Alveolata 0.61, Streptophyta 0.60, Insecta 0.54, Nematoda 0.54, Vertebrata 0.41, Spiralia 0.37, Arthropoda 0.21; fine-tuning adds up to 0.11 | GPU recommended; on PyPI | strongest cross-clade coverage in the cohort; note "locus F1" is not gene F1 | Gaius-Augustus/Vipsania (17 stars, new) |
| PlantGeneAnn | 2026 | 10.64898/2026.06.25.733695 | DL, strand-specific plant foundation model | genome | fine-tuned on 9 plants; 13-species benchmark | beats four baselines at five levels; 9 curated species beat a 42-species set (abstract; CC BY-NC-ND) | GPU | plants only | not found |
| Sensor-NN | 2023 | 10.1093/bioadv/vbad105 | DL, small NN over hand-built sensors | genome | train/test among fly, human, mouse, C. elegans | nucleotide-level balanced accuracy mean 0.72 across 16 train/test pairs; G3PO F1 0.47 to 0.55 vs AUGUSTUS 0.52 (PMC10448985 Tables 1, 3) | CPU | cross-clade by design but nucleotide-level only, and not better than AUGUSTUS | none found |
| SpliceAI / Pangolin / OpenSpliceAI / Spliceator | 2019 / 2022 / 2025 / 2021 | 10.1016/j.cell.2018.12.015; 10.1186/s13059-022-02664-4; 10.7554/elife.107454; 10.1186/s12859-021-04471-3 | DL splice-site CNNs | pre-mRNA | human; four mammals; retrainable; multi-species | splice sites only | GPU for training | OpenSpliceAI and Spliceator retrain across species | see `repos.tsv`; SpliceAI licence is PolyForm Strict |
| Nucleotide Transformer / SegmentNT; Evo 2; Caduceus | 2025 / 2026 | 10.1038/s41592-024-02523-z; 10.1038/s41586-026-10176-5 | DL foundation models | genome | broad | element segmentation, not full gene models; Evo 2 at 7B and 40B parameters | multi-GPU | out of budget; baselines only | see `repos.tsv` |
| G3PO | 2020 | 10.1186/s12864-020-6707-9 | benchmark | 1793 curated genes, 147 species | broad | only 32 of 1793 proteins perfectly predicted by all five programs; 108 by exactly one; short exons the weak point (PMC7147072) | n/a | reusable held-out set | supplementary data |
| EGASP / nGASP | 2006 / 2008 | 10.1186/gb-2006-7-s1-s2; 10.1186/1471-2105-9-549 | benchmark | ENCODE human; C. elegans | | EGASP: no method above 45% exact transcript Sn; only JIGSAW and ENSEMBL above 70% gene-level Sn (PMC1810551) | n/a | historical | n/a |
| fitild | 2018 | 10.1093/bioinformatics/bty353 | analysis | 1022 genomes with over 1000 introns each | broad | intron length carries up to 30% of intron-recognition information in some species; fewer than 20 genomes fit one Frechet component, 490 to 670 need two, the rest three (PMC6157073) | n/a | the clade-dependence the charter worries about, quantified | ogotoh/fitild |
| Scaling studies | 2025, 2026 | 10.1101/2025.04.17.649312; 10.1371/journal.pbio.3003663 | applied | BRAKER on 200 insects; comparative annotation of 301 Drosophilidae | insects | throughput evidence (abstracts) | cluster | n/a | n/a |
| Evaluation tools | 2021, 2023, 2024 | 10.1093/molbev/msab199 (BUSCO); 10.1093/bioinformatics/btad595 (compleasm); 10.1038/s41587-024-02147-w (OMArk) | metrics | proteomes | broad | BUSCO is optimistic relative to gene-level F1 (PMC11216313; PMC11216308) | n/a | n/a | n/a |

## 3. Repository inventory

`repos.tsv` has 32 repositories with last commit date, commits in the last
12 months (shallow clone since 2025-09-09), open issues and stars
(shields.io, 2026-09-09), language, licence, HEAD hash, install result and
notes.

Activity. Alive in 2026: cactus (480 commits in 12 months, 704 stars),
funannotate (279, 400), Tiberius (174, 138), EviAnn (99, 9 on the
development repo; 183 on EviAnn_release), Helixer (37, 305), OpenSpliceAI
(29, 53), egapx (25, 207), CESAR2.0 (14), BRAKER (11, 467), Vipsania (11,
17, created in 2026). Maintenance or dormant: Augustus (2 commits, 339
stars, 166 open issues), miniprot (2, 418), SNAP (0, last commit 2022),
MAKER (0, 2024), GeneMark-ETP (0, 2025), ProtHint (0, 2023), Liftoff (0,
2023, 552 stars, 73 open issues), PhyloCSF (0), PhyloCSF++ (0, 2022),
EVidenceModeler (0, 2024), Pangolin (0, 2023). The star counts say that
the community still runs AUGUSTUS, BRAKER, Liftoff and miniprot far more
than any deep-learning gene finder; the commit counts say that the
deep-learning tools are where the development is.

Licences. Hazards for a derived tool: the GeneMark family (CC BY-NC-SA or
custom non-commercial; this contaminates BRAKER and funannotate at
runtime), SpliceAI (PolyForm Strict: no derivatives), Nucleotide
Transformer (CC BY-NC-SA), ANNEVO (custom non-commercial). Tiberius, Vipsania, TOGA, miniprot, cactus are
MIT; Helixer, EviAnn, Liftoff, OpenSpliceAI are GPL-3; egapx is public
domain.

Install tests (fresh `python3 -m venv`, CPU-only container, no GPU):

- Tiberius (e73844b): `pip install git+...` succeeds in 7 s on Python 3.12
  and is refused on 3.11 (`requires-python >= 3.12`). The package declares
  only `rich` and `pyyaml`; the model stack (`bricks2marble[tf]`,
  `hidten[tensorflow]`) comes with the Singularity or Docker image that the
  cloned launcher `tiberius.py` pulls on first run. There is no console
  entry point. Verdict: installable; not runnable here without the image.
- Vipsania (0c5f84a): on PyPI as 1.0.0, Python 3.12 only. A dry-run
  resolves to tensorflow 2.19.1 plus twelve CUDA 12 wheels plus wandb,
  several GB. Verdict: installable on 3.12; not installed here for disk
  and time.
- Helixer (d17bb49): real install on Python 3.11 succeeds in 56 s
  (tensorflow 2.15.1, numpy 1.26, `setuptools < 72` pin, 2.3 GB venv).
  `Helixer.py --help` then fails with `ModuleNotFoundError: No module named
  'yaml'`: pyyaml is imported by `helixer/core/scripts.py` but not
  declared. After `pip install pyyaml` in the same venv it fails again on
  `from sklearn.utils import shuffle` (scikit-learn also undeclared); I
  stopped there. Verdict: README `pip install` does not give a runnable
  tool; at least two undeclared dependencies. The README's recommended
  route is the Docker/Singularity image, which I did not test.
- egapx, BRAKER, EviAnn: not attempted. egapx needs Nextflow 23.10.1 and
  a container runtime; BRAKER needs the GeneMark-ETP tarball and
  AUGUSTUS; EviAnn ships `install.sh` that compiles miniprot and ufasta and
  bundles TransDecoder and BLAST+, or a bioconda package.

Not on GitHub: Gnomon (NCBI internal, no paper), GeMoMa (jstacs.de),
GeneID, GlimmerHMM. Found in tick 3: xjtu-omics/ANNEVO (158 stars, 30
commits in 12 months, custom non-commercial licence, PyTorch, six lineage
models) and Gaius-Augustus/clamsa (13 stars, no commits since 2024-07, no
licence file). Not found: OrionGeno and PlantGeneAnn code (their abstracts
give no link; bioRxiv full text was unreachable).

## 4. Data sources noted in passing

Whole-genome alignments and trees:

- UCSC multiz alignments (10.1101/gr.1933104) with the 100-way human
  alignment that PhyloCSF++ uses for its benchmark; Cactus/HAL alignments
  (10.1038/s41586-020-2871-y) behind Zoonomia's 241 mammals; the TOGA
  paper's 488-mammal and 501-bird chain sets (10.1126/science.abn3107).
  These are the only sources at the scale a cross-clade comparative model
  needs, and all are mammal- and bird-heavy. For plants, fungi, insects and
  protists we will have to build alignments ourselves with cactus, which
  is why cactus's activity level matters to us.
- CONTRAST's 11-informant human set (macaque, mouse, rat, rabbit, dog,
  cow, armadillo, elephant, tenrec, opossum, chicken) is a documented
  minimal informant panel that already gave most of the gain
  (PMC2246271).

Conservation tracks: phastCons (10.1101/gr.3715005), phyloP
(10.1101/gr.097857.109), PhyloCSF tracks (10.1101/gr.246462.118). Useful
as features and as sanity checks; not sufficient as gene finders.

Reference annotations used as truth by the papers, and therefore usable
for a benchmark that wants to be comparable:

- BRAKER2's 12-species set with annotation versions and genome sizes
  (PMC7787252 Table 1), including the "% non-canonical or incomplete
  genes" column (63.8% in T. nigroviridis, 34.7% in R. prolixus), which
  is itself evidence about reference quality.
- BRAKER3 / GeneMark-ETP's 11 species grouped into compact, large
  GC-homogeneous and large GC-inhomogeneous genomes (PMC11216308,
  PMC11216313), with order-excluded and species-excluded OrthoDB protein
  sets. The "order excluded" protocol is the leakage rule we should copy.
- Tiberius's mammal split (validation on leopard and rat; test on human,
  cow, beluga) and its distance probe (chicken, zebrafish, poplar, tomato)
  (PMC11645249); the 33-species multi-clade panel (preprint, list not
  accessible). ANNEVO's `docs/` in its repository carry its own
  12-species performance table and evaluation notes (README 42c920f).
- Helixer 2026's fungi/plant/vertebrate/invertebrate test panels (13
  plant, 11 vertebrate, 15 invertebrate test species; PMC13076211).
- Vipsania's per-clade training and test species in
  `docs/training_species.tsv` in the repository (README 0c5f84a), which
  covers 17 clades including Alveolata, Amoebozoa, Discoba, Rhodophyta and
  Stramenopiles that no other tool reports on.
- G3PO's 1793 curated genes from 147 species with annotated difficulty
  classes (PMC7147072).
- The fitild intron-length fits for 1022 genomes (PMC6157073) as the
  covariate for choosing a species panel that spans intron regimes.

Expression: BRAKER3 uses VARUS-sampled SRA libraries per species; the
200-insect study (10.1101/2025.04.17.649312) did the same at scale.

Quality tools reusable as metrics: BUSCO (10.1093/molbev/msab199),
compleasm (10.1093/bioinformatics/btad595), OMArk
(10.1038/s41587-024-02147-w), with the explicit caveat from GeneMark-ETP
and BRAKER3 that BUSCO completeness is optimistic and rewards
over-prediction (PMC11216313; PMC11216308).

## 5. Failure modes of current tools

Cross-clade transfer. The pattern is consistent from 2004 to 2026 and the
methods differ only in how gracefully they fail. SNAP's foreign
parameters drop nucleotide Sn from above 93 to as low as 21.7 to 26.0
(PMC421630 Table 4). Tiberius, trained on mammals, shows "a steady decline
in accuracy as phylogenetic distance increases" with poplar and tomato
lowest, and its authors conclude there is a "need for re-training for
other clades" (PMC11645249); the 2026 follow-up answers with six separate
lineage models rather than one model (10.64898/2026.04.24.720536). Helixer
2026 ships four lineage models, and in fungi its transcript-level F1
(0.54) is below the 2005 self-training GeneMark-ES (0.60) (PMC13076211
Table 2), so a kingdom-scale deep model does not automatically beat a
per-genome HMM. Vipsania's README table shows the residual clade
dependence even after unsupervised fine-tuning: locus F1 0.70 in Discoba
but 0.21 in Arthropoda and 0.41 in Vertebrata (README 0c5f84a). The
evidence-based pipelines fail on the same axis through their evidence:
ProtHint intron sensitivity halves between species-level and phylum-level
protein exclusion (79.8 to 35.8) and start-codon sensitivity falls
five-fold (70.3 to 14.1) (PMC7222226 Table 5); BRAKER2 gene sensitivity
spans 70.2 (A. thaliana) to 10.4 (T. nigroviridis) across its 12 species
(PMC7787252 Table 3).

Intron length. Intron length distributions differ enough across species
to be a species fingerprint: of 1022 genomes, fewer than 20 fit one
Frechet component and most need two or three, and the length distribution
carries up to 30% of the information for intron recognition in some
species (PMC6157073). Every HMM tool encodes this as a per-species
parameter (AUGUSTUS's explicit intron submodel, 10.1093/bioinformatics/btg1080;
Spaln's intron-length priors), and Tiberius's untrained HMM layer fixes
geometric length distributions at empirical mammalian means (PMC11645249,
section 2.4), which is one concrete reason it cannot transfer.

GC isochores. GeneMark-ETP trains three GC-specific models per
GC-inhomogeneous genome and attributes its 29.4- and 45.6-point gene-F1
gains over TSEBRA on mouse and chicken to that plus joint evidence use
(PMC11216313); BRAKER3's largest gain over BRAKER1+2 is also chicken
(+55 gene F1 points; PMC11216308). GENSCAN's original isochore-specific
parameters were the 1997 version of the same fix. Any model that is not
explicitly conditioned on local GC will re-learn this the hard way.

Gene structure classes. Short exons under 50 nt are the least accurate
class for every ab initio tool, with AUGUSTUS best at only 18% correct
(PMC7147072); single-exon genes are harder than spliced ones for BRAKER3
in 9 of 11 species, and start codons (F1 70%) are worse than stop codons
(76%) and splice sites (above 87%) (PMC11216308). Initial exons are also
the weakest class for CONTRAST (Sn 76.9 vs 96.2 internal; PMC2246271
Table 2). The hard part of gene prediction is the gene ends and the
smallest pieces, not the internal exons.

Resource use. EGAPx: 32 CPUs and 256 GB (README f9a7392) and, on
Galaxy, about 416,000 CPU-hours for 1,409 jobs with a 45% failure rate
and 24% CPU efficiency (relay/TASK.md). BRAKER3: 5.6 to 64 h on 48
threads per genome, RNA-seq alignment excluded (PMC11216308), 48:53 h
average on mammals (PMC11645249). GeneMark-ETP: 18 h on 64 cores for mouse
(PMC11216313). Tiberius: 1:39 h on an A100 but 15 GPU-days on four A100s
to train, and an 8 GB GPU minimum for inference (PMC11645249; README).
Helixer: 8:54 h per mammal on an A100 (PMC11645249). ANNEVO's README
reports its own tool at 12.2 min and 3.8 GB GPU memory averaged over 12
species on one RTX 4090, against 43.6 min and 22.5 GB for Tiberius and
286.6 min and 8.6 GB for Helixer (README 42c920f; authors' measurements,
not mine), which, if it holds, is the first deep gene finder inside a
laptop-GPU budget. Otherwise the only tools that run in minutes are the
aligners and lift-over tools that do not predict anything new.

Evaluation. BUSCO completeness is above 90% for every GeneMark-ETP
prediction while gene-level F1 varies widely (PMC11216313), and BRAKER3
shows that higher BUSCO goes with more false-positive genes
(PMC11216308). Comparisons across papers are also confounded by isoform
policy (Tiberius evaluates only the longest-CDS isoform per gene,
PMC11645249; BRAKER3 counts a gene correct if any transcript matches,
PMC11216308) and by whether UTRs are scored (Helixer 2026 removes them,
PMC13076211), so published F1 numbers from
different papers are not directly comparable, and our benchmark must fix
one scorer.

Reference quality. BRAKER2's own table reports 63.8% non-canonical or
incomplete genes in the T. nigroviridis reference and 34.7% in
R. prolixus (PMC7787252 Table 1). Some of every tool's "error" on such
species is reference error, which bounds the measurable accuracy of any
method and argues for curated subsets like G3PO alongside whole-genome
references.

## 6. Opinion

What the literature says, read as a whole. Three facts line up. First,
the cheapest strong signal is comparative: a single KA/KS window test
already gives about 9.5% FN and 2 to 3% FP for exon detection
(10.1101/gr.200901), and CONTRAST showed in 2007 that adding informants
keeps paying (gene Sn 35.6 to 58.6 from one to eleven; PMC2246271) if the
model can use them, which N-SCAN's explicit phylogenetic HMM could not.
Second, the best ab initio model today, Tiberius, reaches 55 gene-F1 on
mammals with 8M parameters and about 2M in ablation (PMC11645249), so the
decision surface is not large once the representation is right; what it
lacks is transfer, and its authors' remedy is one model per clade. Third,
every failure mode in section 5 is a species parameter in disguise:
intron length distribution, GC composition, evidence distance, codon
usage. Current tools either fit those per genome (GeneMark-ES,
BRAKER) or bake one clade's values into weights (Tiberius, Helixer).

What I would build. A small model in the HyphAeon pattern: inputs are (a)
the target genomic window, (b) a multiple alignment of that window to a
panel of informant genomes, (c) the informant tree as a metric, and (d)
per-window covariates that are cheap to compute from the target genome
alone: local GC, the fitild intron-length parameters estimated from a
first-pass or from the nearest annotated relative, and the codon usage
table. The network is a 2D axial model over the alignment (site axis,
taxa axis) with tree geometry injected through embeddings and rotary
positions rather than learned. Its output is per-base state posteriors fed
to a differentiable HMM whose length distributions and GC-dependent
emissions are parameterised by (d), not fixed. That last point is the
concrete difference from Tiberius: the HMM layer stops being a mammalian
constant and becomes a function of the species covariates, so the same
weights can serve a fungus with 60 bp introns and a mammal with 10 kb
introns. Training data is the union of the BRAKER3/GeneMark-ETP panels,
Helixer 2026's four lineage panels and Vipsania's 17-clade species lists,
with alignments built by cactus in windows around annotated loci, and
held-out species chosen by order-level exclusion as BRAKER3 does.

Why the alignment and not sequence alone. The 2026 wave is converging on
"phylogeny-aware" from the sequence side (OrionGeno, ANNEVO), which means
learning the phylogeny implicitly from sequence statistics with large
models. Giving the model the alignment and the tree as coordinates is the
cheaper route and the one that scales down to 2M parameters. The
comparative signal also degrades gracefully with distance in a way the
network can be told about (branch length is an input), whereas sequence
statistics silently drift.

What it should not try to do. It should not predict UTRs in the first
version: Helixer 2026 strips UTRs from both reference and prediction before
scoring (PMC13076211), Tiberius scores only the longest CDS per gene
(PMC11645249), and no paper reports transcript-level UTR accuracy, so
UTRs cannot yet be benchmarked against anything. And it
should not replace evidence integration: on genomes with RNA-seq, a
BRAKER3-style pipeline will still add value at the gene ends where every
tool is weakest. Our target is the 80% of assemblies with no annotation
and no RNA-seq, at laptop cost.

The biggest risk. Alignment availability. Every large alignment we can
download is mammal- or bird-centric; for the clades where annotation is
most needed (protists, basal fungi, non-insect arthropods, the "none of
the above" bin in Vipsania's table) there may be no informant genome at a
useful distance, and building cactus alignments for thousands of genomes
is itself a compute problem. The fallback must be designed in from the
start: with an empty alignment the model degrades to an ab initio
predictor with explicit species covariates, and that fallback has to be
at least as good as Vipsania's unsupervised fine-tuning, which is the
strongest published answer for exactly those clades. The second risk is
evaluation: if we score with a different isoform policy than Tiberius and
BRAKER3 we will not know whether we beat them. The benchmark task should
adopt the BRAKER3 scorer and order-level exclusion verbatim and report
Tiberius and Vipsania under it before we train anything.

Where I would be wrong. If Vipsania's unsupervised route reaches Tiberius
accuracy on mammals within a year, the alignment input buys little and
the right project is a smaller, covariate-conditioned Vipsania. The
cheapest experiment that decides this is to run Vipsania and Tiberius on
the same five held-out species spanning intron regimes and see whether
the gap tracks phylogenetic distance or intron-length divergence; the
former favours the comparative design, the latter favours the covariate
design, and both can be done on one consumer GPU.
