# Eukaryotic gene prediction: synthesized review

Task `T-human-006`. Author: `lenin`. First complete draft, 2026-09-10.

This document merges the five independent Phase 1 reviews into one account the
rest of the project cites. The individual reviews stay where they are and stay
authoritative for their own evidence; this is a merge, not a replacement.

| Slot | Task | Author | Artifact |
|---|---|---|---|
| 1 | T-human-002 | `lenin` | [review.md](../../relay/artifacts/T-human-002/review.md) |
| 2 | T-human-003 | `engels` | [review.md](../../relay/artifacts/T-human-003/review.md) + 22 annexes ([index](../../relay/messages/20260910T012607Z-engels-0025.md)) |
| 3 | T-human-004 | `stalin` | [review.md](../../relay/artifacts/T-human-004/review.md) + 22 annexes ([index](../../relay/messages/20260910T004104Z-stalin-0024.md)) |
| 4 | T-human-005 | `marx` | [review.md](../../relay/artifacts/T-human-005/review.md) |
| 5 | T-human-012 | `trotsky` | [review.md](../../relay/artifacts/T-human-012/review.md) |

Companion files:

- [`disagreements.md`](disagreements.md) — every point where the reviews
  reached different conclusions, with a resolution or an open question.
- [`candidates.md`](candidates.md) — the ranked shortlist of ideas for the new
  model, distilled from the five opinion sections.
- [`../refs/refs.bib`](../refs/refs.bib) — the merged bibliography, 126 works.
- [`repos.tsv`](repos.tsv), [`repo-verification.tsv`](repo-verification.tsv),
  [`doi-verification.tsv`](doi-verification.tsv) — the merged software
  inventory and the checks run over it.
- [`annex-dois.tsv`](annex-dois.tsv), [`annex-refs.bib`](annex-refs.bib),
  [`annex-repos.tsv`](annex-repos.tsv) — the works and repository snapshots
  that appear only in the `engels` and `stalin` annex messages (§1).

---

## 1. How this synthesis was made

The five reviews were written blind of each other, as the charter's
Independence clause requires. Three of them (`lenin`, `engels`, `marx`)
recorded machine-readable search and retrieval logs; `stalin` recorded a
retrieval manifest; `trotsky` recorded query counts only.

Merging was done with the scripts under [`scripts/review/`](../../scripts/review),
so that every count in this document is reproducible from the artifacts:

| Script | What it does | Output |
|---|---|---|
| `merge_refs.py` | deduplicates the five bibliographies plus the annex works by DOI, merging fields and recording which reviews carried each work | `docs/refs/refs.bib` |
| `merge_annex_refs.py` | scans the `engels` and `stalin` annex messages for DOIs their own bibliographies never carried, and fetches Crossref metadata | `docs/review/annex-dois.tsv`, `docs/review/annex-refs.bib` |
| `check_conflicts.py` | finds works cited under more than one DOI | printed report |
| `verify_dois.py` | resolves every DOI against Crossref and compares registered titles | `docs/review/doi-verification.tsv` |
| `merge_repos.py` | unifies the five repository inventories, which used five different schemas, plus the two annex snapshots, and flags fields the reviews disagree on | `docs/review/repos.tsv` |
| `verify_repos.py` | re-queries GitHub for the disputed repository rows | `docs/review/repo-verification.tsv` |

**Evidence rule used throughout.** A number appears in the tables below only
when at least one review recorded where it came from — a PMC identifier, a
full-text read, a pinned commit, or a local run. Where reviews disagree on a
number the table carries the value that more reviews independently retrieved,
and the disagreement is itself recorded in [`disagreements.md`](disagreements.md).
Nothing here was re-measured; measurement is `T-human-009`'s work.

**Two verification passes were run for this synthesis** and are the only new
first-hand evidence in it:

1. **DOI resolution (Crossref, 2026-09-10).** 126 merged entries: 110 resolve
   with a matching title, 4 have no DOI by design (repository documentation
   with URLs), 2 are arXiv DOIs that Crossref does not serve, 1 is a journal
   supplement Crossref stores without a title, and **9 either resolve to an
   unrelated paper or do not resolve at all**. All 9 come from one review; see
   [`disagreements.md`](disagreements.md) §1. The 15 annex works added below
   were resolved through the same Crossref pass; none of them is among the 9.
2. **Repository resolution (GitHub API, 2026-09-10).** The 23 repositories
   whose rows the reviews disagreed about, plus 5 disputed slugs. Results in
   [`repo-verification.tsv`](repo-verification.tsv).

**The annexes.** `engels` and `stalin` each submitted a base review plus 22
annex messages carrying corrections to their own drafts, and their submission
covers state that the annexes take precedence over the base text where they
conflict. The qualifications their covers name are carried in §2 and §5 where
they change a claim; §7 lists them.

Their bibliography and inventory additions are now merged.
[`merge_annex_refs.py`](../../scripts/review/merge_annex_refs.py) scans the 49
messages either agent filed against `T-human-003` or `T-human-004`, extracts
every DOI, and subtracts the DOIs the five artifact bibliographies already
carry. **15 works are cited only in the annexes** — 10 cited by `engels` and
7 by `stalin`, with 2 in common, so exclusively 8 `engels` + 5 `stalin` + 2
shared — and all 15 resolve in Crossref
([`annex-dois.tsv`](annex-dois.tsv)). They enter `refs.bib` with
`reviews = {engels (annex)}` or `{stalin (annex)}` and an `annex` field naming
the messages that cite them. `engels`'s cover reports "eight distinct DOI
sources"; `stalin`'s index lists twelve. Those two counts overlap, and five of
`stalin`'s twelve were already in the merged file from other reviews, so the
measured union is 15, not 20.

Two of the 15 are the preprint versions of works other reviews carry as
journal articles (PhyloCSF++, and ClaMSA's "End-to-end learning of evolutionary
models"). Both are kept as separate entries and reported as preprint/journal
pairs by `check_conflicts.py`, which is exactly `engels`'s annex point that
ClaMSA's recovered tables belong to the preprint version.

Two annex messages also carry a repository snapshot that never reached their
author's own inventory: `abacus-gene/paml`
([stalin 0016](../../relay/messages/20260909T164106Z-stalin-0016.md)) and
`Jstacs/Jstacs`, which hosts GeMoMa
([engels 0018](../../relay/messages/20260909T172646Z-engels-0018.md)). Both are
transcribed into [`annex-repos.tsv`](annex-repos.tsv) with their own observation
time, pin and `not_attempted` install status, and merged as `stalin (annex)` and
`engels (annex)`. `stalin`'s index names four snapshot additions; the other
three (OrionGeno, Vipsania, PhyloCSFpp) were absent from `stalin`'s inventory
but already present in the merged table from `lenin`, `marx` and `engels`.

---

## 2. Publications

Approach classes: **GHMM** generalized-HMM ab initio; **COMP** comparative or
alignment-informed; **EVID** evidence-based pipeline; **DL** deep learning;
**HYB** hybrid; **QC/bench** evaluation infrastructure.

Accuracy figures are as reported by the cited authors on their own benchmark.
**They are not comparable across rows** — see §5.8. The `Reviews` column
records which of the five independent reviews carried the work; a work found by
one review only is not weaker evidence, but it is unreplicated search.

### 2.1 Classical ab initio

| Method | Year | Cite | Inputs | Clades | Reported accuracy | Runtime | Reviews |
|---|---|---|---|---|---|---|---|
| GENSCAN | 1997 | [burge1997prediction] | DNA, GC-isochore-specific parameters | human-trained, vertebrate tests | Burset–Guigó set: exon Sn/Sp 0.78/0.81; 243/570 genes exactly right (0.43). The paper itself warns that short single-gene sequences bias this upward | ~X+5 s for X kb, Sun Sparc10 (1997) | all 5 |
| GeneID | 2000 | [parra2000geneid] | DNA | Drosophila, many | G3PO nucleotide F1 0.40 [scalzitti2020benchmark] | fast, CPU | lenin, marx |
| AUGUSTUS | 2003, 2006 | [stanke2003gene; stanke2006augustus] | DNA (+hints) | per-species parameter sets | human h178 exact gene Sn/prec 46%/45%, exon 80%/80%; fly100 gene 52%/27%. G3PO nucleotide F1 0.52, best of five, but only 18% of exons under 50 nt correct [scalzitti2020benchmark] | ~6 min for 1.6 Mb (2003); 2:25 h per mammal on 48 threads [gabriel2024tiberius] | all 5 |
| SNAP | 2004 | [korf2004gene] | DNA + training set | Arabidopsis, rice, worm, fly | own-species nucleotide Sn 93.8–98.1; **foreign parameters collapse** — fly on Arabidopsis Sn 26.0, rice on worm Sn 21.7; bootstrapping recovers 75–96 (Table 4) | ~30 CPU s and 100 MB per Mb, 1 GHz (2004) | all 5 |
| GlimmerHMM | 2004 | [majoros2004tigrscan] | DNA | plant, human | G3PO nucleotide F1 0.45; highest Sn 0.74, lowest Sp 0.43 | fast, CPU | engels, lenin, marx |
| GeneMark-ES | 2005 | [lomsadze2005gene] | DNA only, self-training | fungi, plants, worm, fly, protists | unsupervised matches supervised at nucleotide level (*A. thaliana* Sn/Sp 97.7/96.3 vs 97.2/95.8); gene-level accuracy 5–20% on several later test genomes [brruna2020genemark] | CPU, hours | engels, lenin, marx, stalin (trotsky cites a DOI that is not this paper) |

### 2.2 Comparative methods and single comparative signals

| Method | Year | Cite | Inputs | Clades | Reported accuracy | Reviews |
|---|---|---|---|---|---|---|
| TWINSCAN | 2001 | [korf2001integrating] | target DNA + one informant | mouse target, human informant | original set 1: exact-gene Sn/prec 24.4%/14.4%, exact-exon 68.3%/66.0%, on selected annotated loci | engels, marx, stalin, trotsky |
| KA/KS ratio test | 2002 | [nekrutenko2002ratio] | one pairwise alignment window + frame | human–mouse | 118/1,244 exons fail (**9.5% FN**); mean **2.6% FP** over 24 simulated length/divergence classes. The paper states the test does not resolve exon boundaries | engels, lenin, marx, stalin (trotsky cites a DOI that does not resolve) |
| N-SCAN | 2006 | [gross2006multiple] | MSA + phylogeny, context-dependent substitution and indel models | human, fly | on human CCDS with a mouse informant, as re-measured by [gross2007contrast] Table 1: gene Sn/prec 35.6/25.1, exon 84.2/64.6. "Performs as well using mouse as its only informant as with any combination" | engels, marx, stalin, trotsky |
| CONTRAST | 2007 | [gross2007contrast] | target + up to 11 informants, **no phylogenetic model** | human ENCODE/CCDS | 11 informants: exact-gene Sn/prec 58.6/35.5, exon 92.8/72.5; mouse only 50.8/29.3. A 65% increase in gene sensitivity over N-SCAN. Initial exons are the weak class (Sn 76.9 vs 96.2 internal) | all 5 |
| PhyloCSF | 2011 | [lin2011phylocsf] | codon MSA + tree, empirical codon models | 12 flies, 29 mammals | minimum average error 8% below dN/dS overall, 11% below on 30–180 nt regions | engels, marx, stalin |
| RNAcode | 2011 | [washietl2011rnacode] | MSA | broad | coding detection designed to be species-independent | marx |
| AUGUSTUS-CGP | 2016 | [konig2016simultaneous] | one FASTA per genome + MAF alignment + Newick tree (star tree documented as a fallback) | 12-vertebrate and 12-fly alignments | joint prediction across all genomes as a labelling problem on a graph, solved by dual decomposition | marx, stalin |
| CESAR 2.0 | 2017 | [sharma2017cesar] | alignment + reference exons | mammals | exon-aware realignment inside TOGA | marx |
| ClaMSA | 2022 | [mertsch2022end] | codon MSA + branch-scaled tree; **learned CTMC layer** | vertebrate, fly, yeast | on 12-way vertebrate candidate alignments at 50% sensitivity, FP rate **0.8% versus PhyloCSF 3.5%**. Classifies pre-extracted candidates: a scorer, not a gene finder | engels, marx, stalin |
| TOGA | 2023 | [kirilenko2023integrating] | alignment chains + reference annotation | 488 placental mammals, 501 birds | best-in-class ortholog detection and projection; **cannot find lineage-specific genes** | lenin, marx |
| TOGA2 | 2026 | [malovichko2026accurate] | as TOGA, plus gene-tree reconciliation and UTRs | vertebrates | **513× less memory, 6.1× faster** than TOGA; reports that human-trained deep splice-site models generalize across vertebrates | lenin |
| Liftoff | 2021 | [shumate2021liftoff] | reference annotation + minimap2 | same or close species | high within-species fidelity; not cross-clade | lenin, marx |

### 2.3 Evidence-based pipelines

| Method | Year | Cite | Inputs | Reported accuracy | Runtime | Reviews |
|---|---|---|---|---|---|---|
| MAKER / MAKER2 | 2008 / 2011 | [cantarel2008maker; holt2011maker] | DNA + EST/RNA-seq + proteins + predictors | gene-overlap Sn/prec 89.81%/91.69% on the worm nGASP region (overlap, not exact); introduces Annotation Edit Distance; later beaten by BRAKER2/3 and GeneMark-ETP by 20–28 gene-F1 points | 4.1 h/Mb single core (2008); MPI | all 5 |
| EVidenceModeler | 2008 | [haas2008automated] | predictions + evidence | weighted consensus combiner | light | lenin, marx |
| GeneMark-ET / EP+ | 2014 / 2020 | [lomsadze2014integration; brruna2020genemark] | + RNA-seq introns / + protein hints | **the evidence-distance curve**: ProtHint intron hint Sn falls 79.8 → 35.8 and start-codon Sn 70.3 → 14.1 between species-level and phylum-level protein exclusion, while high-confidence precision stays above 98.8 | ProtHint+EP+ ≈ 5 h, 8 CPUs, 8 GB (fly) | engels, lenin, marx (+trotsky, wrong DOI) |
| Mikado | 2018 | [venturini2018leveraging] | multiple transcript assemblies | scoring-based transcript selection | CPU | marx |
| BRAKER1 / 2 | 2016 / 2021 | [hoff2016braker; brruna2021braker] | RNA-seq (1) or proteins (2) | B2 gene Sn spans **70.2 (*A. thaliana*) to 10.4 (*T. nigroviridis*)** over 12 species; vs MAKER2 on fly, gene F1 59.7 vs 35.9 | ~10 h, 8 CPUs | all 5 |
| GALBA | 2023 | [bruna2023galba] | DNA + close proteins | gene F1 79.5% with the target's own proteome; on Tiberius's mammal panel exon/gene F1 86.2/41.8 | 35:12 h avg, 48 threads | lenin, marx, trotsky |
| GeneMark-ETP | 2024 | [brruna2024genemark] | DNA + RNA-seq + proteins; GC-specific submodels | gene F1 over TSEBRA **+8.2 points (large GC-homogeneous), +39.0 (large GC-inhomogeneous)**; BUSCO stays above 90% while gene F1 varies widely | 3 h fly, 12 h zebrafish, 18 h mouse on 64 cores, **excluding HISAT2/StringTie2** | engels, lenin, marx, stalin (+trotsky, wrong DOI) |
| BRAKER3 | 2024 | [gabriel2024braker] | DNA + RNA-seq + protein db | ~**+20 points mean transcript F1** over BRAKER1/2 on 11 species; chicken +55/+48 gene/transcript F1 over TSEBRA; start-codon F1 70%, stop 76%, splice sites >87% | 5 h 37 min (*A. thaliana*) to 64 h 16 min (mouse), 48 threads, **RNA-seq alignment excluded** | engels, lenin, marx, stalin (+trotsky, wrong DOI) |
| Gnomon / EGAPx | — | [goldfarb2025ncbi; ncbi2026egapx] | assembly + RNA-seq + proteins | **no methods paper and no published gene-level benchmark.** RefSeq reports mean BUSCO completeness 97.3% for its EGAP annotations, which is neither exact-structure accuracy nor an EGAPx benchmark | README: 32 CPUs / 256 GB; fly 144 Mb = 71 CPU-h / 3 wall-h; chicken 1.1 Gb = 425 CPU-h / 5.5 wall-h. **Supported taxa are Chordata, Arthropoda, Echinodermata, Mollusca, Cnidaria, monocots, eudicots; fungi, protists and nematodes are declared out of scope** | all 5 |
| EviAnn | 2026 | [zimin2026efficient] | DNA + transcripts and/or proteins | reports outperforming BRAKER3, MAKER2 and FINDER on identical inputs; mammal under 1 h (abstract; not OA) | one multicore server | marx |
| OMAnnotator | 2026 | [bates2026omannotator] | ab initio + transcript + homology annotations | uses the OMA orthology algorithm as a consensus tie-breaker; improves on its own sources on fly | lenin |
| GeMoSeq | 2026 | [grau2026improved] | RNA-seq | transcript reconstruction by combinatorial enumeration with integral CDS prediction, over seven species | lenin |

### 2.4 Deep learning

| Method | Year | Cite | Inputs | Clades | Reported accuracy | Runtime | Reviews |
|---|---|---|---|---|---|---|---|
| SpliceAI | 2019 | [jaganathan2019predicting] | pre-mRNA sequence | human | top-*k* accuracy ~0.95 for donors/acceptors; splice sites only | GPU, seconds per locus | lenin, marx, trotsky |
| Helixer | 2021, 2026 | [stiehler2021helixer; holst2026helixer] | DNA only | 2021: one vertebrate model over 186 animal genomes, one land-plant model over 51. 2026: four lineage models | 2026 mean transcript F1: **fungi 0.5386, plants 0.4618, vertebrates 0.1977, invertebrates 0.3066**; base-wise phase F1 0.95/0.81/0.88/0.86. In fungi, GeneMark-ES (0.60) beats it. **Scoring caveat: these are GffCompare transcript F1 after longest-protein selection and UTR removal, without `--strict-match`/`-e`, so identical intron chains match despite differing outer boundaries — they do not establish exact CDS start/stop agreement** and are not comparable to the benchmark's exact-CDS endpoint ([engels's endpoint audit](../../relay/messages/20260909T032715Z-engels-0004.md); §5.8, §7) | 8:54 h per mammal on A100; README asks 8–11 GB GPU | all 5 |
| Pangolin | 2022 | [zeng2022predicting] | sequence, multi-tissue | 4 species | multi-species splicing | GPU | lenin, marx, trotsky |
| SegmentNT | 2024/2025 | [dealmeida2024annotating; dealmeida2025annotating] | DNA | fine-tuned on human + 5 animals; 10 animals and 5 plants held out | plant genic-element mean MCC 0.45 vs 0.34 for human-only fine-tuning; frames annotation as instance segmentation of 14 element classes | 3-kb model: 20 h on 8×H100 | engels, lenin |
| Tiberius | 2024 | [gabriel2024tiberius] | DNA + softmasking (+ optional ClaMSA in de novo mode) | trained on mammals; tested human, cow, beluga | three-mammal mean exon/gene F1 **89.7/55.1** vs BRAKER3 83.2/53.7, GALBA 86.2/41.8, Helixer 72.9/19.3, AUGUSTUS 67.3/12.4; **human gene F1 62% vs 21% for the next best ab initio**; de novo (ClaMSA) mode on human 92.6/65.5, **but that comparative number is not human-label-unseen**: Methods 3.7 says the sitewise ClaMSA input generator was trained using human chromosome 17 RefSeq labels, so the held-out-species split protects the Tiberius network, not the whole comparative pipeline (§7; [stalin's supplement audit](../../relay/messages/20260909T084030Z-stalin-0008.md)). **~8M parameters, with a 2M ablation** | 1:39 h per mammal on one A100; training 15 days on four A100s | engels, lenin, marx, stalin (trotsky cites a preprint DOI that resolves to an unrelated paper) |
| Tiberius multi-clade | 2026 | [gabriel2026accurate] | DNA | six lineage models: Mesangiospermae, Fungi, Vertebrata, Insecta, Chlorophyta, Bacillariophyta → "92% of available eukaryotic assemblies" | gene F1 **+12 to +37 over Helixer, +10 to +22 over ANNEVO** across a 33-species panel; approaches BRAKER3 in plants, fungi, diatoms, algae; BRAKER3 keeps a gene-level advantage overall | mean 26 min vs ANNEVO 30, Helixer 178, BRAKER3 2,170 min; 72 CPU threads, A100 for the neural tools | engels, lenin, marx, stalin |
| ANNEVO | 2026 | [zhang2026highly] | DNA; lineage flag | 566-species benchmark (abstract); six released lineage models | "substantially outperforms existing ab initio methods" (abstract; not OA). README: human 19 min on one RTX 4090 for v2.3.3; 12-species mean 12.2 min / **3.8 GB GPU memory** vs Tiberius 43.6 min / 22.5 GB and Helixer 286.6 min / 8.6 GB (authors' measurements) | consumer GPU | engels, lenin, marx, stalin |
| geneML | 2026 | [vader2026geneml] | DNA | nine fungal genomes | gene F1 **67.1 vs BRAKER3-with-proteins 64.9** (recall 64.1 → 69.0 at equal precision); predicts alternative transcripts (41.1% recall / 71.1% precision vs Iso-Seq) | **~6 min/genome on 8 CPU cores** | engels, lenin, stalin |
| GENATATORs | 2026 | [shmelev2026genatators] | DNA-LM embeddings | mammalian training, broader transfer | two negative results: **pretrained DNA-LM embeddings do not capture the features needed for gene segmentation** without task-specific finetuning; **standard per-token/per-sequence metrics fail to capture real annotation quality** | — | lenin, stalin |
| GeneCAD | 2025 | [liu2025genecad] | PlantCAD2 embeddings + transformer + chromosome-scale **CRF** | angiosperms | ~**+9% transcript F1 over Helixer and BRAKER3**, including an allotetraploid; 86% recovery of classical CDS. v0.1.0 shipped a bug that wrecked BUSCO scores | GPU | lenin |
| OrionGeno | 2026 | [liu2026scaling] | DNA, **phylogeny-aware**, long-range | ">5,300 unannotated NCBI genomes" | reports beating the state of the art at exon, gene, protein-sequence and protein-structure level across lineages; also reports coding loci absent from curated references | NVIDIA GPU, compute capability ≥7.0; **non-commercial licence** | lenin, marx |
| Vipsania | 2026 | [krieg2026vipsania] | unannotated genome only, **unsupervised** | 17 clade models plus "other"; 6–13 test species per clade | README locus F1 after fine-tuning: **Discoba 0.70, Fungi 0.66, Alveolata 0.61, Streptophyta 0.60, Insecta 0.54, Nematoda 0.54, Vertebrata 0.41, Spiralia 0.37, Arthropoda 0.21**; fine-tuning adds up to 0.11. Claims to avoid the accuracy drop supervised models suffer with distance | GPU recommended; MIT, on PyPI | lenin, marx |
| PlantGeneAnn | 2026 | [zhang2026plantgeneann] | DNA | 9 plants fine-tuned, 13-species benchmark | beats four baselines at five levels; 9 curated species beat a 42-species set | GPU | marx |
| minisplice | 2026 | [yang2026improving] | sequence | vertebrates + insects | **7,026 parameters**; captures splice signals conserved across phyla; finds mammal/bird-specific GC-rich introns; improves junction accuracy in minimap2/miniprot | trivial | lenin |
| Evo 2 / foundation models | 2026 | [brixi2026genome; dallatorre2025nucleotide; nguyen2023hyenadna; schiff2024caduceus] | DNA, ≤1M-token context | all domains | interpretability shows learned exon–intron boundary representations; **no gene-structure benchmark reported** | very large | lenin, marx, trotsky |
| AlphaGenome | 2026 | [avsec2026advancing] | DNA, 1 Mb | human/mouse regulatory | regulatory variant effects; not gene structure | large | lenin |
| AlphaFold-3 as annotation QC | 2026 | [davison2026promise] | proteins from gene models | 3 fungal/protist species | AF3 scores support **65–84% of manually curated changes**; AF3+Foldseek most discriminative; cheap Protenix-Mini retains the discriminatory power | GPU or cheap | lenin |

### 2.5 Alignment tools, benchmarks and quality control

| Item | Year | Cite | Why it is here | Reviews |
|---|---|---|---|---|
| Exonerate / Spaln / Splign / miniprot | 2005–2024 | [slater2005automated; gotoh2008space; kapustin2008splign; li2023protein; gotoh2024spaln] | spliced aligners bound the accuracy of every EVID pipeline; miniprot is "tens of times faster" than its predecessors; **Spaln's intron-length priors are species parameters** | lenin, marx, trotsky |
| StringTie | 2015 | [pertea2015stringtie] | what BRAKER3 and EGAPx actually consume; not raw reads | lenin |
| EGASP / nGASP | 2006 / 2008 | [guigo2006egasp; coghlan2008ngasp] | EGASP: **no method exceeded 45% exact transcript sensitivity** on human ENCODE | marx |
| G3PO | 2020 | [scalzitti2020benchmark] | 1,793 curated genes from 147 species with difficulty classes. **Only 32 of 1,793 proteins were predicted perfectly by all five programs; 108 by exactly one.** Short exons are the weak point | lenin, marx |
| fitild | 2018 | [gotoh2018modeling] | intron length across 1,022 genomes: **fewer than 20 fit one Fréchet component, 490–670 need two, the rest three**; length carries up to 30% of intron-recognition information in some species | marx |
| BUSCO / compleasm / OMArk | 2021–2025 | [manni2021busco; nevers2025quality] | completeness, not correctness; **BUSCO is optimistic and rewards over-prediction** | lenin, marx |
| Cactus / multiz / phyloP / phastCons | 2005–2020 | [armstrong2020progressive; siepel2005evolutionarily; pollard2010detection] | the alignment and conservation substrate a comparative model needs | lenin, marx |
| Zoonomia | 2023 | [christmas2023evolutionary] | 240+ placental mammal alignments and constraint calls | lenin |
| GAP-MS | 2026 | [abbas2026gap] | mass-spec validation finds **hundreds of peptide-supported loci missing from reference annotations** across nine crops | lenin |
| *Pristionchus* curation | 2026 | [rodelsperger2026lessons] | community curation corrected **>7,500 gene models, ~24% of one nematode annotation** | lenin |
| Selenoprotein annotation | 2026 | [tico2026overcoming] | vertebrate selenoproteins are well annotated for **11% of genes in Ensembl and 5% in NCBI GenBank** | lenin |
| Scaling studies | 2025–2026 | [saenko2025annotation; dhakad2026comparative] | BRAKER on 200 insects; comparative annotation across 301 Drosophilidae | marx |
| Review of the field | 2025 | [djossou2025overview] | extends G3PO over five classical tools plus a gene-model-free NN and Helixer | lenin |

### 2.6 Works reached only by the annexes

These 15 are cited in an `engels` or `stalin` annex message and in no base
review bibliography (§1, [`annex-dois.tsv`](annex-dois.tsv)). They are listed
here rather than folded into 2.1–2.5 because none of them appears in any base
review bibliography, and their substantive integration is still incomplete
(§7). The annexes were filed during the independent review tasks, so their
absence from the frozen base bibliographies is not evidence that blindness was
lost. Two of the 15 (`rogic2001evaluation`, `wei2006ests`) are cited by both
agents; shared citation does not by itself establish that both verified the
same result. The `annex` field in `refs.bib` names the message for each.

| Item | Year | Cite | Why the annex reached for it | Annex of |
|---|---|---|---|---|
| Gene-finder evaluation on mammals | 2001 | [rogic2001evaluation] | the reference population behind several historical accuracy figures both agents had to qualify | engels, stalin |
| Begin at the beginning (5′ UTR prediction) | 2005 | [brown2005begin] | UTR prediction is in the charter's eventual scope and this is where the classical treatment sits | engels |
| AUGUSTUS at EGASP | 2006 | [stanke2006augustusb] | the incomplete-reference scoring that qualifies AUGUSTUS's EGASP numbers | stalin |
| Several pair-wise informants (MARS, a TWINSCAN extension) | 2006 | [flicek2006several] | alternative-transcript prediction from several pairwise informants — a comparative-input design point | engels |
| ESTs improve de novo prediction (N-SCAN_EST) | 2006 | [wei2006ests] | the EST-evidence denominators behind N-SCAN_EST's reported gain | engels, stalin |
| Iterative prediction and pseudogene removal | 2006 | [vanbaren2006iterative] | pseudogene masking changes what counts as a false positive | engels |
| GeneMark-ES for novel fungal genomes | 2008 | [terhovhannisyan2008gene] | the unsupervised-training extension whose transfer claims stalin's annex qualifies | stalin |
| U12-type intron database | 2020 | [moyer2020comprehensive] | U12 introns can have GT-AG termini too, so a terminal-dinucleotide stratum is not a spliceosome class | engels |
| GeMoMa | 2016 | [keilwagen2016intron] | homology-based prediction from intron position conservation; not covered by any base review | engels |
| GeMoMa + RNA-seq | 2018 | [keilwagen2018combining] | the reciprocal-hit evaluation and the reference-conditioned intron limit | engels |
| ClaMSA preprint | 2021 | [mertsch2021end] | the recovered tables belong to this version, not the journal article | engels |
| PhyloCSF++ preprint | 2021 | [pockrandt2021phylocsf] | same version distinction | engels |
| *Stentor* macronuclear genome | 2017 | [slabodnick2017macronuclear] | genuine 15-base introns: a real counterexample to short-intron floors | stalin |
| *Stentor* tiny-intron splicing | 2022 | [nuadthaisong2022insights] | the mechanism, so the counterexample is not an annotation artefact | stalin |
| Genetic codes with no dedicated stop codon | 2016 | [swart2016genetic] | context-dependent termination breaks a decoder that assumes the standard code | stalin |

---

## 3. Software inventory

Merged table: [`repos.tsv`](repos.tsv) — 61 repositories, of which 28 were
recorded by more than one review. Because the five inventories used five
different schemas and two different commit-counting methods (GitHub REST
`since` filters versus `git clone --shallow-since`), the merged table keeps
every value the reviews recorded rather than picking one, and flags the 23
repositories where they differ. Re-queried values for those 23 are in
[`repo-verification.tsv`](repo-verification.tsv) (GitHub API, 2026-09-10).

**Where the development is.** Among gene predictors, `Gaius-Augustus/Tiberius`
is the only one under heavy development: 174 commits in 12 months (four
reviews agree; `trotsky` records 100), MIT, 138 stars, last push 2026-09-07.
`usadellab/Helixer` is second (37–49 commits, GPL-3.0, 305 stars).
`plantcad/genecad` and `hillerlab/TOGA2` are both new and very active.
`ComparativeGenomicsToolkit/cactus` (≈480–485 commits, 704 stars) and
`nextgenusfs/funannotate` (279 commits, 400 stars) are the busiest
infrastructure repositories.

**Where the usage is.** `agshumate/Liftoff` 552 stars and 73–78 open issues
with no commits since 2023; `Gaius-Augustus/BRAKER` 467 stars, **103 open
issues, 11 commits**; `Gaius-Augustus/Augustus` 339 stars, **175 open issues,
2–5 commits**; `lh3/miniprot` 418 stars, 2 commits. `Illumina/SpliceAI` — the
field's de facto splice model — is **archived** with 507 stars.
`Comparative-Annotation-Toolkit` is at zero commits with 116 open issues. The
load-bearing stack of the field is being maintained at a small fraction of the
rate it is being used, and `T-human-009` should expect some baselines not to
build.

**Licensing, which decides what a benchmark may redistribute.** Resolved
against the repositories themselves rather than GitHub's detector, which
several reviews correctly note is not authoritative:

| Terms | Repositories |
|---|---|
| MIT / Apache-2.0 / BSD | Tiberius, Vipsania, TOGA, TOGA2, miniprot, GeneCAD, funannotate, BrentLab/Twinscan |
| GPL-3.0 / AGPL | Helixer, Liftoff, OpenSpliceAI, EviAnn, PhyloCSF |
| Artistic (per repo docs; GitHub reports NOASSERTION or none) | AUGUSTUS, BRAKER, GALBA |
| Public domain + bundled third-party terms | `ncbi/egapx` |
| **Non-commercial** | GeneMark-ETP and the GeneMark family (CC BY-NC-SA 4.0), ProtHint, ANNEVO, OrionGeno, Nucleotide Transformer, SegmentNT |
| **No derivatives** | SpliceAI (PolyForm Strict 1.0) |
| Academic/commercial split | MAKER |
| No LICENSE file found | `Gaius-Augustus/clamsa`, `Gaius-Augustus/TSEBRA` |

The consequence for `T-human-007`: **BRAKER3 cannot be redistributed in a
benchmark image**, because GeneMark-ETP at its core is non-commercial and
ProtHint adds eligibility restrictions; it can only be a locally-run
comparator. The same applies to ANNEVO and OrionGeno — the two published
models closest to this project's design. Only the MIT/Apache/GPL group is
unambiguously redistributable.

**Installation, from four independent attempts.** All four reviews that tried
converge on one finding: *the documented quick-start of the leading tools does
not produce a runnable tool.*

| Tool | Result | Who tried |
|---|---|---|
| miniprot, minisplice | `git clone; make` works in seconds | lenin |
| SNAP | fresh source build plus both packaged examples succeed | engels, stalin |
| Tiberius | `pip install .` installs a launcher that then fails on a missing module; `pip install '.[from_source]'` works but is **6.3 GB**; requires Python ≥3.12 (this project targets 3.11); on a compute-capability-12.0 GPU its pinned TensorFlow must JIT from PTX, which TensorFlow warns "could take 30 minutes or longer" | lenin, engels, marx, stalin |
| Helixer | README `pip install` succeeds in 56 s then fails on **two undeclared dependencies** (pyyaml, scikit-learn); README's own recommended route is the container, and it states install takes 20–30 min for an experienced user, up to 2–3 h otherwise | lenin, marx |
| TWINSCAN | documented `make linux` fails: the i686 target conflicts with an x86-64 host compiler | stalin |
| BRAKER, EGAPx, EviAnn | not attempted; BRAKER needs the GeneMark-ETP tarball plus AUGUSTUS, ProtHint, StringTie2, bedtools, GffRead and partitioned OrthoDB clades; EGAPx needs Nextflow and a container runtime, and its README warns of verified buffer-overrun vulnerabilities in its NCBI C++ toolkit dependencies and recommends running it in a VM | lenin, marx |

The second-order point, made independently by `lenin` and `marx`: the *model*
is small, the *stack* is not. Tiberius is ~8M parameters inside a 6.3 GB
environment. That is a stronger argument for this project's design goal than
the accuracy argument is.

**Repositories that do not exist or are not canonical**, resolved against the
GitHub API on 2026-09-10:

- `weberlab-hhu/Helixer` **redirects** to `usadellab/Helixer` (305 stars,
  GPL-3.0, pushed 2026-07-22) — that is the canonical repository.
- `gglyptodon/Helixer` exists but is a **personal fork with 0 stars, last
  pushed 2024-03-24**. `trotsky`'s conclusion that Helixer development "has
  stalled" is drawn from this fork; see [`disagreements.md`](disagreements.md) §2.
- `washut/contrast` (given as CONTRAST's code link) returns **404**.
- `NBISweden/GAAS` (given as MAKER2's code link) exists but is an unrelated
  annotation toolkit; MAKER is at `Yandell-Lab/maker`.
- `BrentLab/Twinscan` **does exist** (C, MIT, 2 stars, last pushed 2023-06-23),
  which resolves an open item in two reviews.
- `nekrut/axomeme` and `nekrut/scalingPaper` both return 404 to an
  unauthenticated client. These are the charter's sources for the HyphAeon
  design pattern and for the Galaxy EGAPx cost figures. A `question` to the
  coordinator about this has been open since 2026-09-09 (§7).
- There is no canonical MAKER repository beyond `Yandell-Lab/maker` (47 stars,
  last pushed 2024-08-20); MAKER's practical successor is `nextgenusfs/funannotate`.
- GlimmerHMM, GeneID and Gnomon are not on GitHub at all. GeMoMa is, but not
  under its own name: it is a module inside `Jstacs/Jstacs`, which none of the
  five base inventories recorded. `engels`'s annex 0018 pinned it, so the
  merged table now carries it (§1). Repository-wide activity there — 59 commits
  in the window against 12 touching `projects/gemoma` — overstates work on the
  predictor.

---

## 4. Data sources

Not an inventory — that is `T-human-008`, already `done`. This records what
the reviewed papers actually consumed, so the inventory and benchmark tasks
start from evidence.

**Alignments and trees.** UCSC multiz (the 100-way human alignment PhyloCSF++
benchmarks on); Cactus/HAL, demonstrated at >600 amniote genomes
[armstrong2020progressive]; Zoonomia's 241 mammals
[christmas2023evolutionary]; TOGA's 488-mammal and 501-bird chain sets
[kirilenko2023integrating]. **The large alignments the reviewed papers used
are mammal- and bird-heavy**, but that is a statement about this literature,
not about what UCSC serves. The accepted data inventory (`docs/data-sources.md`,
T-human-008, `done`) records deep public alignments on three non-vertebrate
reference assemblies: `dm6/multiz124way` (124 insects, 5.4 GB of MAF plus a
Newick tree), `ce11/multiz135way` (135 nematodes) and `sacCer3/multiz7way`
(7 yeasts), and the KA/KS baseline (T-human-010, `done`) ran on the fly
124-way. So insects, nematodes and yeasts are already covered on a few
reference assemblies; **plants, most fungi and the protists have no
comparable public alignment**, and those must be built with Cactus, which is
why Cactus's activity level matters to this project. The residual risk for
the covered clades is not absence but coverage: one reference assembly per
clade and sparse informant coverage away from it. Distinguish two things that
an earlier draft of this section ran together. Re-anchoring on a genome that is
*already in* the alignment is documented: HAL supports queries relative to an
arbitrary reference or subtree and MAF export with a selectable reference
([README](https://github.com/ComparativeGenomicsToolkit/hal#readme)), so for
HAL-backed resources this is an extraction step with a resource-specific cost.
UCSC's multiz MAF files are reference-anchored and need that conversion or a
re-projection first. A genome that is *absent* from the alignment is the actual
limitation: no choice of reference recovers sequence the alignment never
contained, so a novel assembly must be aligned in. And neither operation is by
itself evidence of good held-out-species prediction. That distinction changes
the cost premise for §3 of `candidates.md`: for insects, nematodes and yeasts
the comparative candidates cost a download, not a Cactus run.
CONTRAST's 11-informant human panel (macaque, mouse, rat, rabbit, dog, cow,
armadillo, elephant, tenrec, opossum, chicken) is a documented minimal
informant set that already captured most of the available gain.

**Conservation.** phyloP and phastCons [pollard2010detection;
siepel2005evolutionarily], with the caveat from the phyloP power analysis:
enough power for strong selection at single nucleotides, moderate selection in
3-bp elements, and weaker or clade-specific selection only in longer elements.
That is exactly the resolution question for splice-site-scale features.

**Reference annotations already used as truth by these papers**, and therefore
the set a comparable benchmark can draw on: BRAKER2's 12 species with
annotation versions and a "% non-canonical or incomplete genes" column
(**63.8% in *T. nigroviridis*, 34.7% in *R. prolixus***); BRAKER3 and
GeneMark-ETP's 11 species grouped into compact / large GC-homogeneous / large
GC-inhomogeneous, with **order-excluded** OrthoDB protein sets; Tiberius's
mammal split (validate on leopard and rat, test on human, cow, beluga) plus its
distance probe (chicken, zebrafish, poplar, tomato); Helixer 2026's four
lineage panels; Vipsania's 17-clade species lists, which cover Alveolata,
Amoebozoa, Discoba, Rhodophyta and Stramenopiles that no other tool reports on;
G3PO's 1,793 curated genes across 147 species with difficulty classes; MANE for
human and mouse; the fitild intron-length fits for 1,022 genomes as the
covariate for choosing a panel that spans intron regimes.

**Expression and proteins.** StringTie-assembled transcriptomes, not raw reads;
VARUS-sampled SRA libraries as BRAKER3 uses them; OrthoDB partitions as
BRAKER2/3's protein input. Both are leakage surfaces: BRAKER3's benchmark is
explicitly parameterized by *how related* the available proteome is, so
"protein evidence" smuggles in the answer when the informant is close.

---

## 5. Failure modes

Consolidated from all five reviews. The bracketed count is how many reviews
independently reached the finding.

**5.1 Per-species fitting is the field's founding assumption, and it has been
dented rather than overturned. [5/5]** Korf stated it as the conclusion of the
SNAP paper — "every genome needs a dedicated gene finder" — after showing
foreign parameters collapse and that the nearest phylogenetic neighbour is not
reliably the best donor [korf2004gene]. GeneMark-ES removed the labelled
training set but not the per-genome fit [lomsadze2005gene]. Helixer was first
to claim a single cross-species model [stiehler2021helixer]; it now ships four
lineage models. Tiberius 2026 ships **six** [gabriel2026accurate]. ANNEVO ships
six [zhang2026highly]. Vipsania ships seventeen plus "other"
[krieg2026vipsania]. Twenty-two years after Korf, "one model" has become "a few
models". **This is the gap the charter exists to close, and the reviews agree it
is the right gap.**

**5.2 Clade exclusion is explicit, not accidental. [4/5]** EGAPx names its
supported taxa and declares fungi, protists and nematodes out of scope
[ncbi2026egapx]. Tiberius's six clade models reach "92% of available eukaryotic
assemblies" [gabriel2026accurate] — and the remaining 8% is precisely the tail a
species-independent method would be judged on. **A benchmark drawn from the 92%
cannot measure species independence**, so `T-human-007` should oversample the
excluded tail deliberately.

**5.3 Intron length is a species fingerprint. [5/5]** Of 1,022 genomes, fewer
than 20 fit a single Fréchet component and most need two or three, and intron
length carries up to 30% of the information used for intron recognition in some
species [gotoh2018modeling]. Every HMM tool encodes this as a per-species
parameter, and **Tiberius's HMM layer fixes geometric length distributions at
empirical mammalian means** [gabriel2024tiberius §2.4] — one concrete,
identified reason it does not transfer. Helixer's headline generalization claim
is specifically that its predictions are less sensitive to genome length
[stiehler2021helixer], which tells you length sensitivity was the dominant
GHMM failure. BRAKER3's gains are largest on large complex genomes
[gabriel2024braker]: same signal from the other side.

**5.4 GC composition is still entangled with clade at the splice signal
itself. [4/5]** GENSCAN's 1997 answer was distinct parameter sets per C+G
region [burge1997prediction] — isochores handled by partitioning the model.
GeneMark-ETP trains **three GC-specific models per GC-inhomogeneous genome** and
attributes its largest gains to that [brruna2024genemark]. In 2026, minisplice
still finds GC-rich introns specific to mammals and birds
[yang2026improving]. Any claim of clade independence must show splice-site
performance stratified by GC, or it is untested.

**5.5 The compute gap is roughly two orders of magnitude, measurable from
published numbers alone. [5/5]** EGAPx: 32 CPUs and 256 GB, 71 CPU-hours for a
144 Mb fly genome and 425 CPU-hours for 1.1 Gb of chicken [ncbi2026egapx].
BRAKER3: 5.6 to 64 h on 48 threads per genome, RNA-seq alignment excluded
[gabriel2024braker]. Tiberius: 1:39 h for a mammal on one A100
[gabriel2024tiberius], ~80× faster than BRAKER3 [gabriel2026accurate]. ANNEVO's
README claims 12.2 min and **3.8 GB of GPU memory** averaged over 12 species on
one RTX 4090, which if it holds is the first deep gene finder inside a
laptop-GPU budget. Note the shape of the EGAPx numbers: CPU-hours grow ~6×
while the genome grows ~7.6× and wall time grows less than 2×, because the cost
is bought with parallelism, not reduced.

The charter's Galaxy figures (416k CPU-hours over 1,409 jobs, 45% failure, 24%
CPU efficiency) are cited by three reviews and **verified by none**; the source
repository is unreachable. If they hold, the interesting quantity is not the
CPU-hours but the **45% failure rate** — nearly half the compute bought nothing.
That is a robustness failure being reported as a cost failure, and documented
runtimes are by construction the runs that succeeded, so no README can supply
it. See §7 and [`disagreements.md`](disagreements.md) §4.

**5.6 Comparative projection methods only find what is already known. [3/5]**
TOGA and CAT project annotation from a reference through an alignment
[kirilenko2023integrating]. They are excellent on conserved genes and
structurally cannot discover clade-specific ones. Any comparative design
inherits this unless it is explicitly built to predict *without* a reference
annotation on the informant side.

**5.7 A coding-region classifier leaves a structured decision problem. [3/5]**
The KA/KS paper states its test does not determine exon or intron boundaries
[nekrutenko2002ratio]; ClaMSA classifies candidate alignments whose strand and
frame are already given [mertsch2022end]. Both `engels` and `stalin` draw the
same methodological conclusion: **a scorer must be evaluated jointly with its
candidate generator**, not credited with a classifier's conditional accuracy,
and candidate recall must be measured before end-to-end accuracy — because if
the correct boundary is pruned, a perfect downstream model cannot recover it.

**5.8 Evaluation is not standardized, and aggregate F1 hides the failures that
matter. [5/5]** No two rows in §2 share a benchmark. Comparisons across papers
are confounded by isoform policy (Tiberius scores the longest-CDS isoform per
gene; BRAKER3 counts a gene correct if any transcript matches), by whether UTRs
are scored (Helixer 2026 strips them — and its GffCompare invocation carries no
`--strict-match`/`-e`, so transcripts match on identical intron chains
regardless of terminal position: stripping UTRs does *not* make the comparison
exact at coding termini, `engels`'s
[endpoint audit](../../relay/messages/20260909T032715Z-engels-0004.md)), and by
reference population (GeneMark-ETP's
sensitivity and precision use different reference populations, so its harmonic
mean need not be a single-confusion-matrix F1 — `stalin`'s annex). Two 2026
papers make the point directly: He & Florea find every method peaks on the exon
class best represented in its training data and "decreases drastically" on the
rest — non-coding, terminal, alternatively spliced, TE-derived
[he2026benchmarking]; GENATATORs finds standard per-token and per-sequence
metrics "fail to capture the challenges of real-world gene annotation"
[shmelev2026genatators]. And the hard classes are known and consistent: **short
exons under 50 nt are the worst class for every ab initio tool, with AUGUSTUS
best at 18% correct** [scalzitti2020benchmark]; initial exons are CONTRAST's
weakest class; start codons (F1 70%) are worse than stops (76%) and splice
sites (>87%) for BRAKER3.

**`T-human-007` should require stratified reporting rather than a single F1, and
should adopt the metric critique from these two papers rather than re-deriving
it.**

**5.9 BUSCO measures the wrong thing for this project. [3/5]** BUSCO
completeness is above 90% for every GeneMark-ETP prediction while gene-level F1
varies widely [brruna2024genemark], and BRAKER3 shows higher BUSCO going with
more false-positive genes [gabriel2024braker]. It is completeness, not
correctness, and it rewards over-prediction. Use it as supporting QC, never as
the endpoint.

**5.10 The reference annotations are themselves wrong, in measurable amounts.
[3/5]** Vertebrate selenoproteins are well annotated for **11% of genes in
Ensembl and 5% in NCBI GenBank** [tico2026overcoming]. Community curation of one
*Pristionchus pacificus* strain corrected **>7,500 gene models, ~24% of the
annotation** [rodelsperger2026lessons]. GAP-MS finds hundreds of
peptide-supported loci absent from references across nine crops
[abbas2026gap]; ANNEVO and OrionGeno claim the same from the prediction side.
BRAKER2's own table reports **63.8% non-canonical or incomplete genes** in the
*T. nigroviridis* reference [brruna2021braker]. A model scored against RefSeq or
Ensembl is partly being scored on its ability to reproduce known errors.
`T-human-007` needs an error bar on the ground truth and a curated
high-confidence subset scored separately.

**5.11 Circular ground truth. [2/5]** `engels` and `stalin` both raise the
sharper version of 5.10: many reference annotations already inherit
homology-based calls, and alignment availability preferentially favours
conserved loci. Apparent improvement can therefore mean reproducing the
reference pipeline's preferences. The mitigation both propose is the same —
hold out taxonomic groups *and* homologous families, keep an
independently-evidenced subset, and report whole-genome numbers alongside the
alignable subset, so that novel or poorly conserved genes cannot vanish from
the denominator merely because the comparative channel is unavailable.

**5.12 Maintenance and version drift are real failure modes. [4/5]** §3. Beyond
"some baselines will not build": ANNEVO's own update history records
architecture, model and decoding changes after release, and its README
benchmarks differ from its paper's; GeneCAD v0.1.0 shipped a BUSCO-wrecking
bug found by users after release. **Pin weights and code before interpreting
any comparison**, and re-measure rather than quoting.

**5.13 Self-training can converge on repeats. [1/5]** In the BRAKER2 study,
unmasked GC-rich tandem repeats in *Xenopus tropicalis* distorted GeneMark-ES
training and needed extra masking [brruna2021braker]. A learned initialization
does not remove this risk.

**5.14 Defaults embody narrow biological assumptions. [1/5]** geneML's README
documents a **maximum intron length default of 400 bases** [vader2026geneml].
Nuclear eukaryotes also include genuine very short introns — *Stentor
coeruleus* has 15-base introns, with a splicing mechanism to match
[slabodnick2017macronuclear; nuadthaisong2022insights] — and ciliate nuclear
genomes with no dedicated stop codon, where termination is context dependent
[swart2016genetic] (`stalin`'s scope annexes). Any cross-clade experiment must
record and test such constraints rather than inherit them. The mirror-image
error is conflating a terminal-dinucleotide stratum with a spliceosome class:
U12-type introns also occur with GT-AG termini, so zero recovered AT-AC introns
does not establish zero U12 recovery [moyer2020comprehensive] (`engels`'s annex
0021). Report motif strata and independently assigned U2/U12 labels separately.

---

## 6. What to build: the joint position

Full ranked shortlist in [`candidates.md`](candidates.md). The five opinion
sections agree on more than they disagree on, and the agreements are the load
bearing part.

**Unanimous (5/5):**

1. **The decision surface is small once the representation is right.** The
   evidence is not rhetorical. KA/KS achieves 9.5% FN / 2.6% FP with one
   statistic on one pairwise alignment [nekrutenko2002ratio]. minisplice
   captures splice signals conserved across phyla with **7,026 parameters**
   [yang2026improving]. Tiberius reaches 55 gene-F1 on mammals with ~8M
   parameters and has a 2M ablation [gabriel2024tiberius]. What is large in
   current models is not the decision surface; it is the machinery for
   manufacturing coordinates from raw DNA.
2. **The decoder must enforce the grammar** — strand, phase, start/stop
   compatibility, splice-site consensus, length ≡ 0 mod 3 — as a structural
   constraint rather than something a network rediscovers per clade.
3. **Intron length must stop being a baked-in constant.** This is the single
   most-cited concrete defect in the strongest current tool.
4. **The benchmark must be fixed before the model exists**, with stratified
   reporting and one scorer, precisely so that a new model cannot be made to
   look good by choosing its benchmark.
5. **Alignment availability is the biggest risk to the whole comparative
   programme.** Every review names it, independently, as the top risk. The
   design assumes a good multiple alignment with a tree at every locus in a
   novel genome; for a newly sequenced genome in a sparsely sampled clade that
   alignment either does not exist or is the most expensive part of the
   pipeline — at which point the cost has moved rather than gone, and the clade
   dependence is back in a new form: *the method works where relatives have
   been sequenced.* Every comparative method since KA/KS has inherited this
   precondition. **A single-sequence fallback path must be designed in from the
   start, not bolted on**, and the degradation curve as informant density and
   divergence degrade must be a required experiment.

**Strong majority (4/5):** the input should be the target window, the aligned
informant rows, the tree with branch lengths, and **explicit masks** separating
"unknown base", "absent taxon" and "alignment gap"; aligned non-coding sequence
must stay in the input, because a model given only coding alignments receives
part of the answer; reading frame is a hypothesis to score, not an input.

**Where the reviews genuinely differ** — decoder family, whether to fix the
parameter budget before proving the benefit, whether the tree helps at all, and
whether UTRs belong in version 1 — is the subject of
[`disagreements.md`](disagreements.md). Those four questions, not the
architecture, are what `T-human-011` has to resolve.

**The one thing nobody else is doing.** Three of the five reviews arrive at the
same conclusion from different directions: the comparative/phylogeny-aware
design space is now occupied — OrionGeno is explicitly phylogeny-aware and has
been run on >5,300 genomes [liu2026scaling]; ANNEVO claims to model joint
evolutionary relationships across 566 species [zhang2026highly]; Tiberius
already ships a comparative (ClaMSA) mode. What is *not* occupied is the
measurement: **nobody has isolated which geometry does the work**, at a fixed
parameter budget, on a neutral benchmark. An ablation separating
tree-as-metric from tree-as-extra-tokens from no-tree, run through one fixed
decoder on one benchmark, would be a real result and is cheap. That, plus the
benchmark, is the contribution this project can honestly claim.

---

## 7. Limits of this synthesis, and open items for the coordinator

1. **The annexes are merged as sources, not as prose.** `engels` and `stalin`
   submitted 44 annex messages between them, and both covers state the annexes
   take precedence over the base drafts. Their covers name the material
   qualifications — Helixer's "exact CDS" wording is too strong for the
   documented GffCompare procedure; ClaMSA's recovered tables belong to the
   preprint version and mixed-clade classification does not establish
   unseen-clade transfer; Tiberius's comparative input generator carries
   supervised human exposure that must accompany its human comparative result;
   GeneMark-ETP's sensitivity and precision use different reference
   populations; nuclear-eukaryote scope includes genuine very short introns and
   alternative genetic codes. Those qualifications are reflected in §2 and §5
   where they change a claim.

   The bibliography and inventory gap is now closed: the 15 annex-only works
   are in `refs.bib` and listed in §2.6, and the two annex repository snapshots
   are in `repos.tsv` (§1). **What is still not done is reading the 44 annexes
   as evidence.** Each is a page or two of source-level auditing — pinned line
   numbers in `codeml.c`, `GeMoMa.java`, PAML's `RemoveIndel` — and this pass
   used them only to (a) find the citations, (b) fix claims their covers
   flagged, and (c) correct the GeMoMa/GitHub statement in §3. An annex
   finding that neither cover surfaced is not in this document. The most
   consequential known example is `stalin` 0016 on PAML: `codeml` overrides
   `cleandata` for pairwise runmodes, so feeding a full multispecies alignment
   to a pairwise fit deletes codon columns that were clean in the target and
   informant pair. `T-human-010` is already `done` and its accepted
   [baseline](../../baselines/kaks/README.md) uses the Nei–Gojobori counting
   estimator rather than `codeml`, so nothing currently executes that path; the
   warning stands for any future `codeml` adapter, and it reaches the record
   through this note rather than through §5.
2. **The `trotsky` review needs a decision.** 9 of its 26 bibliography entries
   carry a DOI that resolves to an unrelated paper or does not resolve at all,
   and several of its table numbers and repository links cannot be traced to a
   source. Nothing from it is used as sole evidence in this document. What to
   do about `T-human-012` is the coordinator's call; see
   [`disagreements.md`](disagreements.md) §1–2.
3. **Two charter sources remain unreadable.** `nekrut/axomeme` (the HyphAeon
   design pattern this project is meant to transfer) and `nekrut/scalingPaper`
   (the Galaxy cost figures) both return 404. The question has been open since
   [20260909T012944Z-lenin-0002](../../relay/messages/20260909T012944Z-lenin-0002.md).
   This blocks `T-human-011` more than it blocked Phase 1: the recommended
   comparison arm is defined by a document no agent has been able to read.
4. **arXiv was not searched directly** by any review; its API returned HTTP 429
   to two runners. OpenAlex `type:preprint` partly compensates. Semantic
   Scholar was used only for citation counts.
5. **Nothing here was measured.** Every runtime and accuracy figure is
   as-reported. `T-human-009` should treat this document as a list of claims to
   check, not as a source of baselines. Exactly one measured row exists to
   check the tables against, and it is not from this task: AUGUSTUS 3.5.0 over
   the whole *S. pombe* genome, one core, **2,076 s wall clock, 403 MB peak
   RSS, 165 s per Mb**, scored with T-human-007's validation procedure
   ([20260910T032948Z-marx-0026](../../relay/messages/20260910T032948Z-marx-0026.md)).
6. **What this document says about alignment availability comes from the
   inventory, not from the reviewed papers.** The literature reviewed here used
   mammal- and bird-centric alignments; UCSC additionally serves deep insect,
   nematode and yeast alignments that none of these papers used (§4). Where the
   two differ, `docs/data-sources.md` (T-human-008, `done`) is authoritative.
