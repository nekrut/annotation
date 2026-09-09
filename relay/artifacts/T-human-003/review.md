# Independent review of eukaryotic gene prediction — engels

**Status: first-run working draft, 2026-09-09 UTC. Task T-human-003 remains
in progress.** This is an independent review; no other agent's review artifact
has been read. Numerical results below are authors' reports, not measurements
made in this run. Publication tables distinguish extracted results from pending
checks. Tiberius's launcher was installed and checked in a fresh virtual
environment; inference and training have not been attempted.

The most promising starting point is a small comparative coding scorer coupled
to explicit gene-structure decoding. The claim that this can work across all
eukaryotes remains a hypothesis. Existing successes do not establish that a
single model will deliver accurate complete gene structures across all clades,
or that good exon classification is sufficient for correct gene reconstruction.

## 1. Search log

Searches were run on 2026-09-09 UTC. The reproducible Europe PMC query URLs,
retrieval timestamps, hit counts, and returned publication metadata are in
[literature-search.json](literature-search.json); web discovery queries are in
[web-search-log.json](web-search-log.json). Results were inspected by
title and primary-source relevance; this is an initial scoping search, not a
completed systematic review. Repeated preprint and journal records describe
the same work and are not independent evidence.

| Search route | Queries / coverage | Counts and limitations |
| --- | --- | --- |
| Europe PMC REST, exact titles and method names | GENSCAN; AUGUSTUS; BRAKER3; SNAP; GeneMark; MAKER; TWINSCAN; N-SCAN; CONTRAST; KA/KS; PhyloCSF; ClaMSA; EGAPx; Tiberius; Helixer; SegmentNT; ANNEVO; geneML | Exact queries and all API counts in the JSON log. Ambiguous searches such as `TITLE:CONTRAST AND gene` were replaced with exact titles. Failed and zero-hit searches are retained. |
| Europe PMC discovery | `(TITLE:"gene prediction" OR TITLE:"gene annotation") AND ("deep learning" OR transformer) AND FIRST_PDATE:[2023-01-01 TO 2026-09-09]` | 36 hits; all returned metadata inspected after an initial page limited to 8 results. Many refer to essential genes, cancer drivers, functional annotation, bacteria, or viruses and are outside the charter. The query misses relevant papers with different title wording. [Log](literature-search.json). |
| Web search, primary sources followed | Tiberius gene prediction differentiable HMM; Helixer cross-species annotation; gene prediction deep learning 2025/2026; exact KA/KS, CONTRAST, N-SCAN, ClaMSA, ANNEVO, and geneML titles | Discovery only; engine results do not supply a stable exhaustive hit count. Publisher pages, PubMed, open manuscripts, and author repositories supply the evidence. |
| GitHub REST | Started with `Gaius-Augustus/Tiberius` and `ncbi/egapx`, then followed named methods and repository links | Inventory and request checksums in [repo-snapshot.json](repo-snapshot.json); executable refresh procedure in [fetch_repo_inventory.py](fetch_repo_inventory.py). This is a selected inventory, not all repositories matching gene prediction. |

Open full text was retrieved through the Europe PMC XML endpoint where
available; response checksums and reading scope are in [reading-log.json](reading-log.json).
PMC web pages sometimes returned browser checks; the XML endpoint
worked for several modern papers. KA/KS XML returned HTTP 404, so its primary
PMC indexed text and abstract were used initially. ANNEVO's published article
is subscription content; only its public abstract, author repository, and
available preprint routes are eligible for this review. No paywalled text or
private repository content is included in these artifacts.

Reading coverage is uneven and explicitly unfinished: Tiberius's original and
clade-extension main texts, SNAP, and the final Helixer main text were read;
PhyloCSF and GeneMark-ETP methods/results received a closer pass. BRAKER3
main-text gaps from truncation were revisited. Other rows currently range from primary
abstracts and documentation to selected main-text passages. Supplementary
tables, exact training splits, and old software distributions require another
pass. Bibliographic metadata is in [refs.bib](refs.bib).

## 2. Publications table

“Pending” means not yet extracted or verified, not that the paper omits the
information. Gene sensitivity, gene precision, nucleotide scores, BUSCO
completeness, and exact transcript F1 are different endpoints. Comparisons
within a paper retain that paper's evidence restrictions and reference set.

### Classical and evidence-based methods

| Method; publication | Class and input | Training / evaluation and transfer evidence | Reported accuracy and benchmark | Runtime / hardware | Code |
| --- | --- | --- | --- | --- | --- |
| GENSCAN, 1997; [DOI](https://doi.org/10.1006/jmbi.1997.0951) | GHMM; genomic DNA; composition-specific parameters and duration models | Human-derived model, evaluated on vertebrate gene sets; no universal cross-clade claim | On Burset–Guigó's gene set: exon sensitivity 0.78, exon precision 0.81; 243/570 genes exact (0.43). The paper cautions that small, simple genes bias this figure upward. | Pending full-text extraction | Historical distribution not yet verified |
| AUGUSTUS, 2003; [DOI](https://doi.org/10.1093/bioinformatics/btg1080) | HMM / duration modeling; DNA; later modes use hints | Human and Drosophila; species-trained models, intron lengths and GC-dependent parameters | Primary abstract reports gains on longer genomic sequences; exact benchmark extraction pending | Pending | [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) |
| SNAP, 2004; [DOI](https://doi.org/10.1186/1471-2105-5-59) | GHMM; DNA; species parameters or bootstrapped training | Arabidopsis, worm, fly, rice; within-species cross-validation and explicit cross-species tests | Arabidopsis gene sensitivity/precision 54.3%/46.8% in five-fold validation. Training on a different species can severely degrade performance (Table 3). | Approximately 30 CPU seconds and 100 MB per Mb on a 1 GHz machine; historical report, not a current measurement. | [SNAP](https://github.com/KorfLab/SNAP) |
| GeneMark.hmm-E / ES, 2011 usage paper; [DOI](https://doi.org/10.1002/0471250953.bi0406s35) | Statistical gene finder; DNA; ES iteratively self-trains | Fungal, plant, and animal uses described. Self-training on the target is adaptation, not fixed-weight transfer. | Detailed results and original development papers pending | Pending | GeneMark distribution and component licenses pending |
| GeneMark-EP+, 2020; [DOI](https://doi.org/10.1093/nargab/lqaa026) | Self-training plus ProtHint protein-derived splice/start/stop hints; EP+ enforces reliable hints during prediction | Detailed species and exclusion rules pending | Primary abstract reports gains over ES/ET, especially in large genomes; numerical extraction pending | Pending | Component source / binary boundary pending |
| GeneMark-ETP, 2024; [DOI](https://doi.org/10.1101/gr.278373.123) | Iterative GHMM; genome, RNA-seq, proteins; high-confidence genes seed GC-specific training | Seven plant/animal genomes; target-specific retraining. Large-genome evaluation uses reference intersections/unions, requiring care with denominators. | Gene F1 gains over TSEBRA: 8.2 points for large GC-homogeneous genomes, 39.0 for large GC-inhomogeneous genomes; group averages in this study | Fly 3 h, zebrafish 12 h, mouse 18 h on 64 CPU cores; **HISAT2 and StringTie2 excluded**; order-excluded proteins | [GeneMark-ETP](https://github.com/gatech-genemark/GeneMark-ETP) |
| BRAKER2, 2021; [DOI](https://doi.org/10.1093/nargab/lqaa108) | GeneMark-EP+ and AUGUSTUS; genome and proteins | Automated target-specific training; no RNA-seq requirement | Original benchmark extraction pending | Pending | [BRAKER](https://github.com/Gaius-Augustus/BRAKER) |
| BRAKER3, 2024; [DOI](https://doi.org/10.1101/gr.278090.123) | GeneMark-ETP, AUGUSTUS, TSEBRA; genome, short-read RNA-seq, protein database | 11 reference species; species-excluded and order-excluded protein databases; additional novel genomes assessed separately | About 20 percentage points higher mean transcript F1 than BRAKER1/2; this is a within-study result, not a universal margin. | 5 h 37 min (Arabidopsis) to 64 h 16 min (mouse), 48 threads, Xeon E5-2650 v4; **RNA alignment excluded**, target-model training included. | [BRAKER](https://github.com/Gaius-Augustus/BRAKER) |
| MAKER, 2008 issue / 2007 online; [DOI](https://doi.org/10.1101/gr.6743907) | Repeat finding, EST/protein alignment, ab initio prediction and evidence integration | Planarian Schmidtea mediterranea proof of principle; retraining and benchmark details pending | Original numerical results pending. Do not assign later MAKER2 figures to this release. | Pending | [MAKER](https://github.com/Yandell-Lab/maker) |
| Gnomon / NCBI pipeline / EGAPx; current [documentation](https://github.com/ncbi/egapx/blob/f9a7392b6f60d0d24f28db5c4eb69b0e1984616b/README.md) | Alignment chaining and HMM completion; genome, taxid, RNA-seq; automatically selected proteins and HMMs | Supported animal/plant groups; fungi, protists, and nematodes explicitly excluded | Documentation is not a controlled accuracy paper. Dedicated DOI/preprint for this row not located yet. | Example fly workload: 71 CPU h / 3 wall h on mixed AWS Batch instances; evidence load is specified in the README. Prerequisite guidance lists 32 CPUs / 256 GB. | [EGAPx](https://github.com/ncbi/egapx); caller repository does not by itself describe every bundled component |

### Comparative methods and simple signals

| Method; publication | Class and input | Training / evaluation and transfer evidence | Reported accuracy and benchmark | Runtime / hardware | Code |
| --- | --- | --- | --- | --- | --- |
| TWINSCAN, 2001; [DOI](https://doi.org/10.1093/bioinformatics/17.suppl_1.s140) | Extends GENSCAN with separate conservation models for exons, introns, splice sites and UTRs | High-throughput mouse genomic sequences with human homology; exact split pending | Primary abstract reports improved exact-gene/exon sensitivity and precision; numerical extraction pending | Alignment and prediction costs separate; numbers pending | Historical distribution pending |
| N-SCAN, 2006; [DOI](https://doi.org/10.1089/cmb.2006.13.379) | Comparative model; multiple genome alignment and phylogenetic relationships, context-dependent substitutions and indels | Human and Drosophila whole-genome applications; this is not evidence of one species-independent parameter set | Primary abstract reports improved whole-genome prediction; exact metrics pending | Pending | Historical distribution pending |
| CONTRAST, 2007; [DOI](https://doi.org/10.1186/gb-2007-8-12-r269) | Discriminative boundary classifiers plus global structure model; multiple informants without explicit phylogenetic model | Human evaluation; multiple-informant ablations require close examination | Abstract reports 65% more exactly reconstructed human coding structures and 46% fewer missed exons than its comparator; relative improvements, not percentage-point F1 gains. | Pending | Historical distribution pending |
| KA/KS test, 2002; [DOI](https://doi.org/10.1101/gr.200901) | Comparative exon classifier; aligned homologous sequence and frame | Selected human–mouse exons; random-sequence negative controls | Exon false-negative rate 9.5%; mean false-positive rate 2.6% over 24 simulated length/divergence classes with 1,000 pairs each. These are not genomic false-discovery or complete-gene error rates. | Pending; alignment is additional | Original automation not yet recovered |
| PhyloCSF, 2011; [DOI](https://doi.org/10.1093/bioinformatics/btr209) | Coding-versus-noncoding empirical codon-model likelihood comparison; alignment and tree | About 50,000 exon-length coding/noncoding regions from 12 fly genomes; four-fold validation | Minimum average error 8% lower than dN/dS overall and 11% lower on 30–180 nt regions. These are relative error reductions in alignment classification. | Per-region branch-scale optimization; exact runtime pending | [PhyloCSF](https://github.com/mlin/PhyloCSF) |
| PhyloCSF++, 2022; [DOI](https://doi.org/10.1093/bioinformatics/btab756) | Comparative coding score implementation and annotation utilities | Compatibility, benchmark, and source snapshot pending | Pending; do not assume implementation speed changes establish new biological accuracy | Pending | Repository follow-up pending |
| ClaMSA, 2022; [DOI](https://doi.org/10.1093/bioinformatics/btac028) | Learned continuous-time Markov chain layer and neural classifier; codon MSA and scaled tree | Vertebrate and fly alignment classification; repository also supplies yeast training data | Abstract reports fourfold fewer false positives at the same true-positive rate than existing methods. This is candidate alignment classification. | Pending | [ClaMSA](https://github.com/Gaius-Augustus/clamsa) |

### Deep learning and hybrid models

| Method; publication | Class and input | Training / evaluation and transfer evidence | Reported accuracy and benchmark | Runtime / hardware | Code |
| --- | --- | --- | --- | --- | --- |
| Helixer prototype, 2021; [DOI](https://doi.org/10.1093/bioinformatics/btaa1044) | Neural base classification from DNA | Separate vertebrate and land-plant models; evaluates transfer within broad groups | Base-wise/subgenic endpoints; full transcript reconstruction was not the prototype's output. Exact table extraction pending. | Pending | [Helixer](https://github.com/usadellab/Helixer) |
| Helixer with HelixerPost, 2025 online / 2026 issue; [DOI](https://doi.org/10.1038/s41592-025-02939-1); [2023 preprint](https://doi.org/10.1101/2023.02.06.527280) | CNN/recurrent model plus separate HMM decoding; unmasked DNA | Separate clade models; 45 test species. Exact primary CDS comparison removes UTRs. | Table 2: mean transcript F1 0.5386 fungi, 0.4618 plants, 0.1977 vertebrates, 0.3066 invertebrates. Comparator species subsets differ for AUGUSTUS. | Human just under 8.5 h; Oryza brachyantha 27 min; single-threaded pipeline with GPU, workstation specification pending | [Helixer](https://github.com/usadellab/Helixer) |
| Tiberius, 2024; [DOI](https://doi.org/10.1093/bioinformatics/btae685) | CNN/LSTM/differentiable HMM; DNA, optional masking or ClaMSA | Mammal training; human/cow/beluga tests; taxonomic exclusions | Mean exact gene/exon F1 55.1%/89.7%; human comparative gene F1 65.5% | Mean inference 1 h 39 min, A100 80 GB/48 CPU threads; training 15 days/four A100s. Comparative preprocessing separate. | [Tiberius](https://github.com/Gaius-Augustus/Tiberius) |
| Tiberius clade extension, 2026 **preprint**; [DOI](https://doi.org/10.64898/2026.04.24.720536) | Same core hybrid, updated implementation and separately trained lineage models | Added plant, fungal, vertebrate, insect, chlorophyte and diatom models; comparisons across supported clades | All-tools comparison uses 33 species; gene score counts a locus correct if any predicted isoform matches, unlike the earlier longest-CDS-only comparison. BRAKER3 retains a gene-level advantage overall. | Mean 26 min Tiberius versus 2170 min BRAKER3; 72 CPU threads for all, additional A100 80 GB for neural tools. This is not a hardware-normalized speedup. | [Tiberius](https://github.com/Gaius-Augustus/Tiberius) |
| ANNEVO, 2026; [DOI](https://doi.org/10.1038/s41592-026-03036-7); [open preprint](https://doi.org/10.21203/rs.3.rs-6402260/v1) | Mixture-of-experts genomic model and structure decoding; DNA | Public abstract describes 566 species; complete split audit pending | Published abstract and newer repository benchmarks are distinct evidence. Exact published structure metrics pending. | Current repository supplies hardware-specific timings; extraction and version matching pending | [ANNEVO](https://github.com/xjtu-omics/ANNEVO); noncommercial license |
| SegmentNT, 2025; [DOI](https://doi.org/10.1038/s41592-025-02881-2) | Pretrained NT-v2 500M encoder plus 63M-parameter U-Net; multilabel nucleotide segmentation | Fine-tuned on human plus five animals; ten animals and five plants held out from fine-tuning. Pretraining membership still needs audit. | Plant genic-element mean MCC 0.45 versus 0.34 for human-only fine-tuning; this is not exact-gene F1. | Initial 3-kb model training: 20 h on eight H100s; later length fine-tuning adds cost. Inference runtime pending. | [Nucleotide Transformer / SegmentNT](https://github.com/instadeepai/nucleotide-transformer) |
| geneML, 2026 **preprint**; [DOI](https://doi.org/10.64898/2026.05.18.725946) | Fungal neural gene annotation with alternative-transcript support | Fungal scope; primary full-text and split audit pending | Indexed abstract describes gene/transcript benchmarks; no numerical claim accepted here before checking comparator configuration | Pending | [geneML](https://github.com/hexagonbio/geneML) |

The tables intentionally retain a documentation-only Gnomon/EGAPx row. The
task's DOI/preprint requirement for every row is not yet satisfied, and the
task is not ready for review. A DOI for a component or unrelated RefSeq update
would not establish an EGAPx accuracy benchmark.

## 3. Repository inventory

[repos.tsv](repos.tsv) records public metadata at the UTC snapshot time shown
in each row. Commit counts use the preceding calendar year and default-branch
history, including merge commits. Open issue counts exclude pull requests;
GitHub's combined issue/PR total was adjusted using its open-PR endpoint.
Request URLs and response hashes are retained in [repo-snapshot.json](repo-snapshot.json).
These API calls are not an atomic snapshot and may differ slightly if a
repository changes between requests.

Tiberius's documented `pip install .` and `python tiberius.py --list_cfg`
both passed in a fresh Python 3.13 virtual environment without system site
packages. This verifies **the launcher only**. TensorFlow, containers, weights,
GPU operation and gene prediction were not tested. The isolated source tree and
environment were removed afterwards. Exact commands, dependency versions, exit
codes and timing are in [tiberius-launcher-check.json](tiberius-launcher-check.json),
with a [repeatable check](check_tiberius_launcher.py). Other installation
results remain **not attempted**; finding a README or container is not success.

GitHub's license detector is not authoritative. Inspected statements are
recorded in [repo-audit.json](repo-audit.json); the inventory preserves the
original detector value separately where overridden. Remaining `unknown` and
`NOASSERTION` entries need inspection. EGAPx's [license file](https://github.com/ncbi/egapx/blob/f9a7392b6f60d0d24f28db5c4eb69b0e1984616b/LICENSE)
distinguishes NCBI government-written code from third-party components.
ANNEVO's [author documentation](https://github.com/xjtu-omics/ANNEVO) states
noncommercial restrictions; it should not be called permissively open source.
SNAP's current [license](https://github.com/KorfLab/SNAP/blob/4ad1e957cd8e68b63857cc1cb3380d39a7b518b1/LICENSE)
says MIT, whereas its original article described GPL: the release and current
source must be distinguished. No legal interpretation beyond the documents'
stated terms is attempted.

The current Tiberius [package manifest](https://github.com/Gaius-Augustus/Tiberius/blob/e73844bfc7170665bc6632f1806a05dfe59fae8a/pyproject.toml)
identifies version 2.0.7 and requires Python >=3.12. The charter's Python 3.11
project code should therefore keep any later Tiberius benchmark environment
separate and pinned. This is a compatibility observation, not an installation
failure. Two charter-linked repositories returned anonymous API 404 responses;
no private contents were retrieved or used for this review.

## 4. Data sources noted in passing

These are leads, not an approved benchmark or a downloaded dataset.

| Source | Potential use | Evidence and next check |
| --- | --- | --- |
| NCBI RefSeq assemblies and annotations | Genome/structure pairs | Used by [Tiberius](https://doi.org/10.1093/bioinformatics/btae685); recover accession/release versions and split membership. |
| Zoonomia alignment | Comparative coding evidence | Tiberius uses a 64-species ClaMSA subset; recover preprocessing and coordinate mapping. [Paper](https://doi.org/10.1093/bioinformatics/btae685). |
| OrthoDB and SRA | Protein hints and transcript support | [BRAKER3](https://doi.org/10.1101/gr.278090.123) documents excluded taxonomic groups and RNA-seq selection. A protein's presence in a database can convey target-derived annotation even after the target species is removed. |
| ClaMSA distributed train/test sets | A candidate-alignment baseline before full gene decoding | Author [repository](https://github.com/Gaius-Augustus/clamsa) provides fly, vertebrate and yeast data routes, with branch-length scaling requirements. Verify split independence and origin of negative examples. |
| ANNEVO reference/data release lists | Additional cross-clade benchmark candidates | Publisher [data-availability statement](https://www.nature.com/articles/s41592-026-03036-7) points to RefSeq and Ensembl releases. Recover membership rather than treating all listed species as held out. |
| UCSC multiz/Cactus, phyloP/phastCons | Whole-genome aligned noncoding sequence and conservation controls | Charter-suggested leads; availability, licensing, download sizes, and exact locus-fetch procedures remain unverified here. |

## 5. Failure modes and implications

**Transfer is a biological distribution problem.** SNAP's cross-species tests
show that GC composition can matter more than the nearest phylogenetic
neighbor when choosing foreign parameters. Its bootstrap results are useful,
but compact genomes and filtered gene sets limit the inference. A fixed
multi-clade model needs explicit evaluation of shifts in codon usage and
splice signals. [SNAP](https://doi.org/10.1186/1471-2105-5-59).

**A coding-region classifier leaves a structured decision problem.** KA/KS
requires informative homologous comparisons, and its exon-level results do
not resolve boundaries, gene chaining, splice alternatives, or alignment
absence. Its false-negative rate must not be presented as universal gene
error. [KA/KS](https://doi.org/10.1101/gr.200901). ClaMSA and PhyloCSF should
therefore be evaluated both as candidate scorers and as inputs to a shared
decoder. This last sentence is a proposed experimental design.

**Intron duration and exact boundaries matter.** Original Tiberius uses geometric
lengths and a boundary-sensitive objective. It has limitations for spliced
start codons and junction-spanning stops. Audit these cases as well as exact
boundaries; high base scores do not ensure correct genes.
[Tiberius](https://doi.org/10.1093/bioinformatics/btae685).

**Evidence is uneven across genes.** BRAKER3's performance depends on expression
support and protein availability. Its published runtime excludes RNA-seq
alignment. Consequently, an evidence-based benchmark must record input
coverage and preprocessing cost; it cannot treat the produced BAM files as
free. [BRAKER3](https://doi.org/10.1101/gr.278090.123).

**The metric can reverse the interpretation.** BUSCO measures recovery of a
conserved gene set and does not penalize all extra genes or exon-boundary
errors. Primary transcript selection and any-isoform matching are different
tasks. The Tiberius clade-extension paper explicitly explains the advantage
multi-isoform BRAKER3 gets from the latter definition. Report both endpoints
and the allowed output isoforms. [Clade-extension preprint](https://doi.org/10.64898/2026.04.24.720536).

**Software and weights change the comparator.** ANNEVO's [update history](https://github.com/xjtu-omics/ANNEVO/blob/main/docs/update_history.md)
records architecture, model, and decoding changes after its original release.
Pin weights and code before interpreting comparisons made at different dates.
Activity counts are maintenance clues, not evidence of prediction accuracy.

## 6. My opinion: what I would build

I would begin with a compact comparative scorer and a constrained decoder,
with a DNA-only path available wherever alignment evidence is missing. The
first question is whether the comparative information helps reconstruct
complete genes on genuinely held-out species. Reducing the parameter count
before establishing that benefit would optimize the wrong thing. This is my
design recommendation, not a result established by the reviewed papers.

The input should include the target sequence, aligned homologous sequences,
the alignment's tree with branch lengths, and explicit masks for missing or
unreliable observations. Repeat annotations and expression evidence should be
optional channels with their provenance retained. Aligned intronic and
intergenic sequence must remain available: a model given only known coding
alignments would receive part of the answer in its input. The pipeline must
also define how candidate alignments and reading frames are constructed
without access to the target's reference annotation.

I would encode strand symmetry, codon phase, start/stop constraints, and splice
compatibility explicitly. The sequence encoder should learn local evidence;
a structured decoder should assemble valid paths. A candidate splice graph
could connect distant exons without requiring dense attention across every
intervening nucleotide. Its danger is candidate pruning: if the correct
boundary is removed, a perfect downstream model cannot recover it. Candidate
recall must therefore be measured before end-to-end accuracy, with difficult
boundaries and long introns examined separately.

Tree information is a promising bias, but an attractive geometric embedding
is not itself a validated evolutionary model. I would compare an explicit
tree likelihood or learned CTMC layer, a distance-aware small encoder, and an
encoder without tree information while keeping the decoder fixed. ClaMSA is
already a closely related comparator. Any embedding should be tested for
distance distortion, changed taxon order, removed taxa, and different branch
length scales. The charter's HyphAeon analogy motivates this experiment; its
particular efficiency claims are not independently verified in this draft.

For introns, I would compare learned duration conditioning with a simple
heavy-tailed or mixture duration prior rather than hard-coding the length
regime of a familiar mammal. The model should be allowed to express uncertainty
where sequence divergence, alignment gaps, or assembly errors eliminate useful
evidence. Species identity must not become a lookup table for hidden
species-specific parameter sets. I would inspect whether a shared model
actually transfers, including to lineages with unusual composition and sparse
reference annotations.

The biggest risk is circular ground truth. Many reference annotations already
inherit homology-based calls, and alignment availability preferentially favors
conserved loci. Apparent improvement could therefore mean reproducing the
reference pipeline's preferences. I would hold out taxonomic groups and
homologous families where feasible, retain an independent evidence-supported
subset, and report the entire genome alongside the alignable subset. Novel or
poorly conserved genes cannot disappear from the denominator merely because
the comparative channel is unavailable.

The decisive experiment would compare KA/KS, PhyloCSF/ClaMSA, a DNA-only neural
baseline, and the proposed scorer through the same decoder and evaluation
rules. Report exact coding exons, transcripts, loci, start/stop sites and
stratified errors, alongside BUSCO as supporting quality control. Account for
alignment retrieval or construction, feature generation, inference, decoding,
and training separately. I would recommend prototyping only after an accuracy
gain survives these checks within the coordinator's compute budget.

## Remaining work before submission

Local format and bibliography checks are recorded in
[verification.json](verification.json). They verify the artifacts' format and
links, not the remaining scientific claims or task completion.

- Complete open-full-text and supplementary-table extraction, especially
  classical methods, GeneMark, Helixer, PhyloCSF, SegmentNT, ANNEVO and geneML.
- Resolve the Gnomon/EGAPx publication requirement without attaching an
  unrelated DOI; expand GeneMark and AUGUSTUS historical coverage.
- Recover historical code links and review licenses directly; add newly found
  active repositories and perform small, pinned fresh-environment checks.
- Reconcile newer benchmarks by model revision, held-out taxa, evidence,
  gene/transcript definitions, and preprocessing costs.
- Bibliography parsing already passes BibTeX; repeat after adding references.
  Announce completion in a future tick only when the task is ready for review.
