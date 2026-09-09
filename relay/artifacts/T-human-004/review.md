# Independent review of eukaryotic gene prediction — stalin

Status: working draft from the first run, 2026-09-09 UTC. Task remains
`in_progress`. This is an independent review; no other agent's review
artifact was read. Evidence extraction and installation checks remain
unfinished. Numerical results below are reported by their cited authors,
not reproduced here, unless explicitly identified as repository observations.

## 1. Search log

The search starts with the charter's methods and follows primary papers,
open manuscripts, software documentation, and repository links. It is a
targeted review, not a systematic review with a complete screening flow.

Exact Europe PMC queries, timestamps, result counts, retained bibliographic
metadata, and API URLs are in [literature-metadata.json](literature-metadata.json).
The inclusion manifest is [included-publications.json](included-publications.json).
Inclusion in the bibliography means a relevant source was identified; it
does not mean its full text or supplements have been completely reviewed.

| Europe PMC search group | Returned hits | Metadata records retrieved | Interpretation |
| --- | ---: | ---: | --- |
| Classical names, initial | 235 | 100 | AUGUSTUS collided with unrelated names and clinical studies; refined below |
| Evidence pipelines | 8 | 8 | BRAKER3, GeneMark-ETP, GeneMark family, MAKER, RefSeq |
| Comparative names, initial | 81,933 | 100 | CONTRAST was too ambiguous; refined below |
| Recent gene prediction/annotation plus deep learning, transformer, or language model; 2023-01-01 through 2026-09-09 | 43 | 43 | Includes out-of-scope functional, essential-gene, and cancer-gene prediction; these are not method evidence here |
| Classical names with gene/genome qualification | 18 | 18 | Still requires relevance screening and preprint deduplication |
| Comparative names with gene-prediction qualification | 17 | 17 | Some unrelated medical matches remain; selected primary methods only |
| Exact DOI/title follow-up for SNAP, N-SCAN, AUGUSTUS, Helixer, and comparative AUGUSTUS | 8 | 8 | Repairs missing historical titles; identifies additional follow-up work |

Counts are API observations, not numbers of included studies. Broad-query
metadata was reduced to gene/genome/PhyloCSF title matches after retrieval;
the original response-size counts were retained. Duplicate preprints and
journal articles are not independent replications.

Web discovery included these query groups on 2026-09-09: Tiberius and EGAPx
official repositories; the original and multi-clade Tiberius titles; the
KA/KS paper's exact title and PMC identifier; ClaMSA; Helixer's journal DOI;
ANNEVO; geneML; GENATATORs; N-SCAN with Gross and Brent; SNAP with Korf;
GENSCAN with Burge; and AUGUSTUS's original intron submodel. A final GENSCAN follow-up used
`Burge Karlin 1997 prediction complete gene structures human genomic DNA pdf GENSCAN`
and located the primary paper on a university course page. Searches for
recent work explicitly included 2025 and 2026. Web result counts are not
treated as stable hit counts. Reviews, search summaries, and third-party
mirrors were discovery aids; claims below use the original research or
maintainers' documentation.

Access notes: some PMC browser pages presented a challenge, while the
Europe PMC full-text XML endpoint worked for several open articles.
KA/KS was read from the publisher's freely accessible PDF. ANNEVO's journal
page exposed an abstract and metadata but identified the article as
subscription content; its open preprint is a separate version requiring
follow-up. Some full-text endpoints returned 404/403. No access restriction
was bypassed and no paper full text is committed. GitHub's unauthenticated
API was rate-limited; the configured `gh` client retrieved public repository
metadata. The old Helixer repository URL redirects to `usadellab/Helixer`.

## 2. Publications table

`Pending` means not extracted or independently checked during this run,
not absent from the paper. An exon classifier, a per-base segmenter, a
complete CDS predictor, and an evidence pipeline solve different tasks.
There is deliberately no ranking made by pooling their reported scores.
Bibliographic details are in [refs.bib](refs.bib).

| Method; year; primary source | Approach and required inputs | Training/evaluation domain and transfer evidence | Reported accuracy and benchmark | Runtime/hardware and limits | Code or implementation source; reading status |
| --- | --- | --- | --- | --- | --- |
| GENSCAN; 1997; [Burge & Karlin](https://doi.org/10.1006/jmbi.1997.0951) | Ab initio GHMM-like gene structure model with GC-dependent parameters; genomic DNA | Human training, historical vertebrate tests; training overlap in the Burset/Guigo set is explicitly acknowledged | Burset/Guigo: 243/570 complete genes correct, exact-exon sensitivity/precision 0.78/0.81; independent GeneParser sets I/II exact-exon sensitivity 0.79/0.76 | Historical typical time for X kb: X + 5 seconds on a Sun Sparc10 | [Author's server](http://hollywood.mit.edu/GENSCAN.html); [primary paper](https://www.cs.rice.edu/~devika/comp470/papers/burge97prediction.pdf) benchmark definitions, overlap caveat and runtime read |
| AUGUSTUS; 2003, 2006; [original model](https://doi.org/10.1093/bioinformatics/btg1080), [alternative transcripts](https://doi.org/10.1093/nar/gkl200) | HMM/GHMM; DNA and species parameters; original intron model uses an explicit short-length distribution and geometric tail | Original human and Drosophila evaluation; species-dependent parameter estimation | Original exact-gene/exon tables pending extraction | Pending; historical timings need their original CPU context | [Augustus](https://github.com/Gaius-Augustus/Augustus); original author manuscript methods inspected; later paper queued |
| SNAP; 2004; [Korf](https://doi.org/10.1186/1471-2105-5-59) | Semi-HMM; DNA; supervised or bootstrapped parameters | Arabidopsis, rice, C. elegans, Drosophila; foreign parameters can be poor even for a nearby species | Exact transfer matrix pending; qualitative failure of nearest-relative parameter selection verified | Paper reports about 30 CPU seconds and 100 MB per 1 Mb on a 1 GHz machine; historical, not a current benchmark | [SNAP](https://github.com/KorfLab/SNAP); abstract and selected methods/results read |
| GeneMark-ES; 2005; [Lomsadze et al.](https://doi.org/10.1093/nar/gki937) | Ab initio GHMM with iterative Viterbi self-training; unannotated genomic DNA | Species adaptation through retraining; not evidence for one frozen cross-clade model | Reported comparisons to supervised methods; per-species metrics pending | Pending | [GeneMark](https://exon.gatech.edu/GeneMark/); abstract and selected results read |
| GeneMark-ETP; 2024; [Brůna et al.](https://doi.org/10.1101/gr.278373.123) | Hybrid pipeline; genome, RNA-seq, homologous proteins; high-confidence genes initialize iterative training | Plant and animal genomes; GC-specific models for heterogeneous genomes | Direct comparison and input-database choices inspected; absolute per-species F1 values pending | Drosophila, zebrafish, mouse: 3, 12, 18 h on 64 CPU cores; HISAT2 and StringTie2 excluded | [GeneMark-ETP](https://github.com/gatech-genemark/GeneMark-ETP); abstract, methods/results passages and runtime read |
| BRAKER1/2/3; 2016/2021/2024; [BRAKER1](https://doi.org/10.1093/bioinformatics/btv661), [BRAKER2](https://doi.org/10.1093/nargab/lqaa108), [BRAKER3](https://doi.org/10.1101/gr.278090.123) | Evidence pipelines: RNA-seq in BRAKER1, proteins in BRAKER2, both in BRAKER3; genome-specific statistical training | BRAKER3 evaluated on 11 species, with protein relatedness controlled | BRAKER3 reports about 20 percentage points mean transcript-F1 improvement over BRAKER1/2; individual panels pending | BRAKER3: 5 h 37 min for Arabidopsis to 64 h 16 min for mouse, 48 Xeon E5-2650 v4 threads; RNA read alignment excluded, parameter training included | [BRAKER](https://github.com/Gaius-Augustus/BRAKER); BRAKER3 selected full-text passages read; earlier versions queued |
| Gnomon / EGAPx; [RefSeq 2025 context](https://doi.org/10.1093/nar/gkae1038) and [current EGAPx documentation](https://github.com/ncbi/egapx) | Evidence pipeline; genome and taxid select protein/HMM resources; RNA-seq recommended; alignment chaining plus HMM prediction | Maintainers exclude fungi, protists and nematodes; this is a documented product boundary | No matched accuracy experiment established in this run | README lists 32 CPUs and 256 GB RAM; a documented configuration, not a measured minimum or per-genome cost | [EGAPx](https://github.com/ncbi/egapx); README inspected; RefSeq article detailed extraction pending |
| MAKER; 2008; [Cantarel et al.](https://doi.org/10.1101/gr.6743907) | Genome, EST/protein alignments, repeat masking and trainable SNAP predictions | C. elegans nGASP test region and Schmidtea annotation project | Table 1: gene-overlap sensitivity/precision 89.81%/91.69%; these are overlap metrics. Exact-transcript metrics use only the confirmed reference subset | Schmidtea: 4.1 h/Mb including all compute stages on a single-core 2 GHz Mac with 2 GB RAM; historical measurement | [MAKER](https://github.com/Yandell-Lab/maker); [author-hosted full paper](https://planaria.stowers.org/wp-content/uploads/2017/04/maker.pdf) benchmark, methods and timing inspected |
| TWINSCAN; 2001; [Korf et al.](https://doi.org/10.1093/bioinformatics/17.suppl_1.s140) | Extends GENSCAN with feature-specific conservation models; target DNA and informant alignments | Mouse targets with human informants; accession-level eightfold cross-validation, not a held-out species | Original set 1: exact-gene sensitivity/precision 24.4%/14.4%, exact-exon 68.3%/66.0%; selected annotated loci | Original wall-time/hardware not extracted | [TWINSCAN/N-SCAN](https://github.com/BrentLab/Twinscan); [original paper](https://web.stanford.edu/class/cs273a/papers.spr07/lecture12/Twinscan.pdf) methods and tables inspected |
| N-SCAN; 2006; [Gross & Brent](https://doi.org/10.1089/cmb.2006.13.379) | Comparative model using genome alignments, phylogeny, context-dependent substitutions and indels | Human and Drosophila whole-genome predictions; not demonstrated universal transfer | Original absolute metrics pending; later CONTRAST comparison is separate evidence | Pending | [TWINSCAN/N-SCAN](https://github.com/BrentLab/Twinscan); public publisher abstract and shared software README read |
| CONTRAST; 2007; [Gross et al.](https://doi.org/10.1186/gb-2007-8-12-r269) | Discriminative boundary classifiers plus global gene model; target DNA and multiple informants; no explicit phylogenetic model | Human cross-validation; tests adding informants | Table 1: with 11 informants, exact-gene sensitivity 58.6%, reported specificity 35.5%; N-SCAN with mouse: 35.6%, 25.1%. These specificity values are prediction precision, not true-negative rates | Pending | Original source availability pending; OA abstract, results passages and tables read |
| KA/KS exon test; 2002; [Nekrutenko et al.](https://doi.org/10.1101/gr.200901) | Comparative coding classifier; aligned human–mouse exon candidates and reading frame; significant purifying-selection test | Empirical known exon pairs plus simulated negatives; no whole-genome structure benchmark | 118/1,244 exons fail the test (9.5%); simulated-pair false-positive mean 2.6%; denominators and data sources differ | codeml/F3×4 analysis; end-to-end runtime not extracted | [Freely accessible publisher PDF](https://genome.cshlp.org/content/12/1/198.full.pdf); results and limitations read |
| PhyloCSF; 2011; [Lin et al.](https://doi.org/10.1093/bioinformatics/btr209) | Comparative coding classifier; MSA, tree and coding/noncoding codon models | Original 12-species Drosophila alignment benchmark | ROC-style coding discrimination; exact operating-point metrics pending; no gene reconstruction score | Pending | [PhyloCSF](https://github.com/mlin/PhyloCSF); abstract and selected method discussion read |
| ClaMSA; 2022; [Mertsch & Stanke](https://doi.org/10.1093/bioinformatics/btac028) | Differentiable phylogenetic CTMC likelihoods plus classifier; candidate codon MSA and branch-length tree | Vertebrate, fly, yeast alignments; mixed-clade and cross-clade experiments | On 12-way vertebrate candidate alignments at 50% sensitivity, false-positive rate 0.8%, versus PhyloCSF 3.5%; not whole-genome gene accuracy | Runtime table pending; alignment construction is additional work | [ClaMSA](https://github.com/Gaius-Augustus/clamsa); open publisher methods/results inspected |
| AUGUSTUS-CGP; 2016; [König et al.](https://doi.org/10.1093/bioinformatics/btw494) | Comparative complete-gene prediction; detailed inputs pending | Multi-genome annotation; training/transfer details pending | Pending | Pending | [Augustus](https://github.com/Gaius-Augustus/Augustus); metadata and official dataset page inspected; full-text endpoint unavailable |
| Helixer; 2021, 2025; [base model](https://doi.org/10.1093/bioinformatics/btaa1044), [gene-model paper](https://doi.org/10.1038/s41592-025-02939-1) | CNN/bidirectional LSTM per-base predictions with HMM postprocessor; genomic DNA | Fungi, plants, vertebrates, invertebrates; clade-specific pretrained weights | 2025 Table 2 mean transcript F1: fungi 0.5386, plants 0.4618, vertebrates 0.1977, invertebrates 0.3066; verify precise gffcompare conventions before comparison to CDS-only metrics | Own-paper runtime extraction pending | [Helixer](https://github.com/usadellab/Helixer); 2025 methods and tables read; 2021 full extraction queued |
| Tiberius; 2024; [Gabriel et al.](https://doi.org/10.1093/bioinformatics/btae685) | CNN/LSTM plus differentiable HMM; genome, optional repeat track; comparative mode adds ClaMSA | Mammalian training; human, cow and beluga tests; longest-CDS isoform evaluation | Mean exact-gene/exon F1 55.1%/89.7%; comparative human gene F1 65.5%; within this paper's setup | Mean inference 1 h 39 min, A100 80 GB with 48 CPU threads; training 15 days on four A100s; input preparation is separate | [Tiberius](https://github.com/Gaius-Augustus/Tiberius); main article inspected, supplements pending |
| Tiberius multi-clade extension; 2026 preprint; [Gabriel et al.](https://doi.org/10.64898/2026.04.24.720536) | Same general hybrid design; separate lineage models | Adds flowering plants, fungi, vertebrates, insects, green algae and diatoms; not one universal model | Comparative panel covers 33 species; exact per-clade tables pending. The paper reports no Tiberius test species in its training/validation sets | Mean runtimes on shared 33-species panel: Tiberius 26 min, ANNEVO 30 min, Helixer 178 min, BRAKER3 2,170 min; 72 CPU threads, neural tools additionally A100 80 GB | [Tiberius](https://github.com/Gaius-Augustus/Tiberius); open preprint methods/results read; version-sensitive comparison |
| ANNEVO; 2026; [Zhang et al.](https://doi.org/10.1038/s41592-026-03036-7); [2025 preprint](https://doi.org/10.21203/rs.3.rs-6402260/v1) | Mixture-of-experts genomic model; DNA input | Journal abstract reports evaluation across 566 species; split and clade-model audit pending | Journal full text not accessed; current README comparisons are separate author-reported experiments | Journal-matched hardware/runtime extraction pending; do not mix current README speed with earlier paper scores | [ANNEVO](https://github.com/xjtu-omics/ANNEVO); public abstract and current README read; custom noncommercial software licence |
| geneML; 2026 preprint; [Vader et al.](https://doi.org/10.64898/2026.05.18.725946) | Deep-learning fungal gene predictor; DNA; README describes candidate transcript decoding | Nine fungal reference genomes; training/test partition audit pending | Author abstract: gene F1 67.1 versus 64.9 for BRAKER3 with protein hints; full benchmark configuration still requires verification | Author abstract: about 6 min/genome on an 8-core CPU; processor model and complete timing scope pending | [geneML](https://github.com/hexagonbio/geneML); author abstract via Europe PMC and README read; intron-length defaults warrant testing |
| GENATATORs; 2026 preprint; [Shmelev et al.](https://doi.org/10.64898/2026.06.17.732686) | Fine-tuned DNA language models; sequence segmentation, splice/CDS filtering | Mammalian training and broader transfer evaluation; crop-selection and split audit pending | Scope includes UTRs and lncRNAs; aggregate gene-segmentation claims cannot stand in for this charter's CDS-only endpoint | Appendix H hardware/cost extraction pending | [Author model collection](https://huggingface.co/collections/AIRI-Institute/genatator); abstract and main-text sections inspected |

## 3. Repository inventory

[repos.tsv](repos.tsv) records thirteen repositories, with source observations
and commit identifiers in [repo-metadata.json](repo-metadata.json). Each row
has its own UTC snapshot and annual window. Commit counts use the default
branch's reachable commits at the recorded SHA and GitHub's `since`/`until`
filters, including merge commits. The last-commit field is a committer
timestamp, not `pushed_at`. Open issue counts exclude pull requests. Stars,
issues, and activity are observations, not quality scores.

SNAP passed a fresh source build and both documented example commands at
the pinned commit; [snap-install-check.json](snap-install-check.json)
records the commands, system toolchain, timestamps and input checksums.
The environment was a temporary checkout with restricted executable search
path, using host GCC/make; it was not a clean operating system or container.
No accuracy or genome-scale performance claim follows from these examples.
TWINSCAN's documented `make linux` failed in another fresh checkout: the
i686 target conflicted with the host compiler's x86-64 mode, and pkg-config
also reported missing GLib metadata. The source was not patched;
[twinscan-install-check.json](twinscan-install-check.json) preserves the
command, compiler and diagnostics. This is a host-specific installation
result, not evidence that the software cannot run elsewhere. The remaining
installation fields say `not_attempted`.
Documentation inspection is not installation success. Core execution,
dependency licences, optional aligners, model downloads and container
requirements still need separate checks. The inventory preserves the API
licence identifier alongside manual notes where available. ANNEVO and
GeneMark-ETP describe noncommercial terms; MAKER's licence file distinguishes
academic and commercial use. These observations link to the inspected
versions in the inventory, and dependency terms need separate review.
`NOASSERTION` is not a finding of unrestricted use.

The current Tiberius README recommends an 8 GB GPU and requires Python
3.12 or newer, whereas this project's code target is Python 3.11. Its
runtime would therefore need a separate environment when benchmarked.
The source and paper version must both be pinned before comparing current
releases. [Current Tiberius documentation](https://github.com/Gaius-Augustus/Tiberius)

## 4. Data sources noted in passing

These are leads, not a completed data inventory or authorization to
download large datasets.

| Source | Potential reuse and unresolved checks |
| --- | --- |
| [Tiberius data/methods](https://doi.org/10.1093/bioinformatics/btae685) | RefSeq genome/annotation accessions, test exclusions and human-referenced comparative inputs; retrieve supplements and separate annotation labels from informant evidence |
| [Helixer methods and supplements](https://doi.org/10.1038/s41592-025-02939-1) | RefSeq and Phytozome genomes; curated accession lists and clade train/validation/test splits; reference quality and reuse terms need checking |
| [ClaMSA datasets](https://doi.org/10.1093/bioinformatics/btac028) | UCSC MULTIZ vertebrate alignments and Cactus fly/yeast alignments; candidate-exon sampling and known frame are benchmark conditioning, not realistic gene discovery inputs |
| [AUGUSTUS datasets](https://bioinf.uni-greifswald.de/augustus/datasets/) | Historical training/test sets and comparative annotations; distinguish prediction-derived labels and transferred annotation from independent truth |
| [ANNEVO data availability](https://doi.org/10.1038/s41592-026-03036-7) | RefSeq and Ensembl accession lists; recover exact release and split information from open supplements or preprint |
| [GENATATOR dataset](https://huggingface.co/datasets/shmelev/genatator-segmentation-dataset) | Segmentation examples and benchmark conventions; assess how selected windows differ from uninterrupted genome scans |

## 5. Failure modes and comparison traps

**Coding discrimination is only part of gene prediction.** The KA/KS paper
explicitly says its test does not determine exon/intron boundaries. Its
simulated negatives do not establish performance on conserved regulatory
sequence, pseudogenes, repeats or a genome-wide scan over candidate frames.
Short exons and unsuitable divergence reduce power. These are reasons to
retain the test as a baseline, with a separate structural endpoint.
[Nekrutenko et al.](https://doi.org/10.1101/gr.200901)

**Conservation alone is ambiguous.** ClaMSA compares evolutionary coding
models and discusses conserved noncoding sequence as a source of false
positives. Its candidates already specify strand and frame. My inference:
a downstream gene decoder must be evaluated jointly with candidate
generation rather than credited with an exon classifier's conditional
accuracy. [Mertsch & Stanke](https://doi.org/10.1093/bioinformatics/btac028)

**Phylogenetic proximity does not fully determine parameter transfer.**
SNAP's cross-species experiments implicate composition as well as
relatedness. The original AUGUSTUS model explicitly handles short intron
lengths and GC-dependent parameter estimation. A proposed universal model
must test these axes, not just average across a phylogenetic tree.
[SNAP](https://doi.org/10.1186/1471-2105-5-59),
[AUGUSTUS](https://doi.org/10.1093/bioinformatics/btg1080)

**Separate clade weights do not demonstrate frozen-model universality.**
The Tiberius extension trains lineage-specific models. Its comparison
also identifies unequal exposure of test species to competitors' training
or validation data. Published ranking should be read together with those
splits and software versions.
[Tiberius extension](https://doi.org/10.64898/2026.04.24.720536)

**Curated locus tests can overstate genome-wide results.** The original
GENSCAN article acknowledges train/test overlap in one benchmark and
reports a separate independent test. It also distinguishes short, usually
single-gene sequences from long genomic contigs. These are longstanding
benchmark hazards, relevant when evaluating modern selected-window tests.
[GENSCAN](https://doi.org/10.1006/jmbi.1997.0951)

**Metrics and output scope change conclusions.** Helixer reports both
base-phase and structure metrics, whose values differ substantially.
GENATATORs evaluates outputs including noncoding transcripts and UTRs,
which exceed this charter's initial coding-gene scope. The original
Tiberius comparison selects one coding isoform per gene. Use separately
defined CDS, transcript and isoform endpoints.
[Helixer](https://doi.org/10.1038/s41592-025-02939-1),
[GENATATORs](https://doi.org/10.64898/2026.06.17.732686),
[Tiberius](https://doi.org/10.1093/bioinformatics/btae685)

**Runtime boundaries are inconsistent.** BRAKER3's cited experiment omits
RNA read alignment; GeneMark-ETP's runtime omits HISAT2 and StringTie2.
GPU inference and CPU pipeline times answer different hardware questions.
Future measurements should report preparation, training/adaptation,
alignment, inference and decoding, plus their end-to-end total.
[BRAKER3](https://doi.org/10.1101/gr.278090.123),
[GeneMark-ETP](https://doi.org/10.1101/gr.278373.123)

**Defaults can embody narrow biological assumptions.** geneML's current
README documents a maximum intron length default of 400 bases. This is a
configuration finding, not a demonstrated error rate. Any cross-clade
experiment must record and test such constraints.
[geneML documentation](https://github.com/hexagonbio/geneML)

## 6. Opinion: what I would build

This is a provisional design hypothesis, not a prototype decision. I would
build a compact comparative feature model coupled to an explicit decoder
for coding gene structure. The first scientific question is whether a
small model can preserve exact exon boundaries and reading frame under
changes in clade, composition and intron length. A good average coding
score is insufficient evidence.

The inputs should be the target DNA window on each strand, locally
available genomic alignment rows, their tree distances, and masks that
distinguish unknown bases, absent taxa and alignment gaps. Coding frame is
a hypothesis to score, not an input supplied by a known annotation. I
would preserve nucleotide coordinates near splice sites while exposing
codon substitution features for each candidate frame. Noncoding aligned
sequence should remain available: removing everything outside an inferred
ORF would remove evidence needed for its boundaries.

I would compare two evolutionary encoders before choosing the more
complex one: a small bank of codon-likelihood features with a learned
readout, and the charter's axial site/taxon attention design. ClaMSA is a
particularly relevant existing comparator because it already trains
phylogenetic likelihoods discriminatively. It makes the novelty question
sharper: improvement must come from better structural prediction,
robustness or total cost, not simply from adding a tree to a classifier.
[ClaMSA](https://doi.org/10.1093/bioinformatics/btac028)

The decoder should enforce strand, phase, start/stop compatibility and
splice transitions, while making intron duration a visible modeling
choice. I would test an explicit duration model or a mixture of durations
against a simple geometric distribution. Candidate pruning must expose
its recall ceiling; otherwise a cheap decoder can look excellent after
discarding difficult genes. The initial parameter target would follow the
charter's small-model ambition, with capacity and input preparation cost
reported separately. [Charter](../../TASK.md)

For the proposed tree geometry, I would test stability to taxon ordering,
missing taxa and rescaling of branch lengths. MDS coordinates are not
unique under rotations or sign changes. If the encoder depends on those
arbitrary coordinates, the representation needs a stable convention or
invariant operations. A low-dimensional embedding is an approximation
whose distortion should be measured, rather than assumed harmless.
These are mathematical design checks I propose, not failures established
for the charter's reference implementation.

The biggest risk is that usable alignments and annotation labels decide
which loci enter the dataset. A method can appear general while seeing
only conserved, well-aligned genes and learning another pipeline's errors.
I would therefore keep an explicit alignment-free input path, evaluate
missing-alignment loci, hold out whole clades, and use independently
supported reference subsets. Comparative preprocessing must count toward
cost. The recommendation would change if a simple boundary model plus
codon likelihoods matches the larger encoder on the same blinded tests;
in that case the simpler model should advance.

## Continuation checklist

1. Finish full-text and supplementary-table extraction for the historical
   methods, MAKER, Gnomon/EGAPx, and newer neural tools. Replace `Pending`
   cells with supported measurements or precise statements that the
   source does not report them.
2. Audit native test splits, reference isoform rules, endpoint definitions,
   software versions, and the availability/cost of alignment inputs.
3. Extend the repository inventory to remaining implementations and
   dependencies; attempt bounded fresh-environment installs and distinguish
   installation from successful model execution.
4. Reconcile bibliography metadata with publisher dates, inspect relevant
   open supplements, and tighten the one-page recommendation. Submit only
   after the task's definition of done is met, with its required relay
   announcement on a later tick.
