# Candidate ideas for the new model, ranked

Task `T-human-006`. Author: `lenin`, 2026-09-10.

Distilled from the five opinion sections. This is input to `T-human-011`, not a
substitute for it: the charter asks that task for two or three candidate model
designs with a recommendation, and the coordinator for a `decision`.

**Ranking criterion.** Decision-relevance per unit of compute, weighted by
whether the ground is occupied. The reviews converge on the observation that
the phylogeny-aware deep-model design space now has well-resourced groups in it
— OrionGeno has been run on >5,300 genomes [liu2026scaling], ANNEVO reports a
566-species evaluation [zhang2026highly], Tiberius already ships a comparative
mode — while the *measurement* that would say which part of the design does the
work has not been done by anyone. Candidates that produce a decision rank above
candidates that produce an architecture.

Provenance column: which reviews proposed or supported the idea.

| # | Candidate | Cost | Decides | From |
|---|---|---|---|---|
| 1 | Neutral benchmark with stratified metrics and a leakage protocol | already scoped as `T-human-007` (`done`) | everything downstream | 5/5 |
| 2 | Informant-degradation curve | small | whether the comparative programme is viable at all | 5/5 |
| 3 | Tiberius ab initio vs ClaMSA mode, replicated on mammals | ~1 day, 1 GPU | whether tree/alignment input buys anything given a good sequence model | lenin |
| 4 | Fixed-decoder encoder ablation | weeks, 1 GPU | which geometry does the work | engels, stalin, lenin |
| 5 | Covariate-conditioned decoder | weeks | whether "one model, many clades" is achievable without per-clade weights | marx |
| 6 | Fixed-encoder decoder ablation (HMM vs CRF vs splice graph) | weeks | isoforms and the intron-length regime | trotsky, engels, stalin |
| 7 | Unsupervised fallback arm (Vipsania-style) | weeks | what to do where no alignment exists | lenin, marx |
| 8 | Tiny-model probe at minisplice scale | days | how low the parameter budget can honestly go | lenin |

---

## 1. The neutral benchmark, with stratified metrics — unanimous, and first

**What.** One held-out species panel, one scorer, one isoform policy, one
leakage rule, run over Tiberius (multi-clade), Helixer, ANNEVO, BRAKER3,
Vipsania, GeneCAD and the KA/KS floor, with results reported **stratified** and
not as a single F1.

**Why it is first.** Every review reaches this independently, and two 2026
papers supply the argument as a citation rather than an opinion: He & Florea
find every method peaks on the exon class best represented in its training data
and degrades drastically on the rest [he2026benchmarking]; GENATATORs finds
standard per-token and per-sequence metrics fail to capture real annotation
quality [shmelev2026genatators]. Without this, no later result means anything —
and the reviews note the converse risk explicitly: a new model can be made to
look excellent by choosing its benchmark, so the benchmark should be fixed
before the model exists, precisely so that it cannot be.

**The strata the reviews name.** Exon length (short exons under 50 nt are the
worst class for every ab initio tool, AUGUSTUS best at 18% correct
[scalzitti2020benchmark]); exon position (initial exons are CONTRAST's weakest
class); single-exon versus spliced genes; GC regime; intron-length regime, using
the fitild fits over 1,022 genomes as the covariate [gotoh2018modeling]; and
alignment coverage. Report both the whole genome and the alignable subset, so
that poorly conserved genes cannot vanish from the denominator when the
comparative channel is unavailable.

**Constraints the reviews already settled.** Adopt BRAKER3's order-level
exclusion protocol verbatim and one fixed scorer, or published numbers will not
be comparable to ours. BUSCO is supporting QC only — it is above 90% for every
GeneMark-ETP prediction while gene-level F1 varies widely
[brruna2024genemark], and higher BUSCO tracks *more* false-positive genes
[gabriel2024braker]. Carry an error bar on the ground truth and a curated
high-confidence subset ([`README.md`](README.md) §5.10). Note that BRAKER3,
ANNEVO and OrionGeno cannot be redistributed in a benchmark image
([`README.md`](README.md) §3), so they are locally-run comparators.

`T-human-007` is already `done`; this section is a checklist against which its
output should be read, not a request to redo it.

## 2. The informant-degradation curve

**What.** Measure how accuracy falls as informant density and divergence fall:
fix a target species, vary the number of informants and their phylogenetic
distance, and plot.

**Why.** All five reviews name alignment availability as the top risk, in
almost the same words. The design assumes a good alignment with a tree at every
locus in a novel genome. For a newly sequenced genome in a sparsely sampled
clade that alignment either does not exist or is the most expensive part of the
pipeline — the cost has moved rather than gone, and the clade dependence is
back in a new form. Every large alignment that can be downloaded is mammal- or
bird-centric; for protists, basal fungi and non-insect arthropods there may be
no informant at a useful distance.

The evidence that this curve is steep already exists on the evidence side:
ProtHint's intron-hint sensitivity falls 79.8 → 35.8 and start-codon sensitivity
70.3 → 14.1 between species-level and phylum-level protein exclusion
[brruna2020genemark]. There is no reason to expect the alignment side to behave
better.

**What it decides.** Whether Phase 4 builds a comparative model with a
single-sequence fallback, or a single-sequence model with comparative
refinement. Those are different projects, and the reviews are honest that the
answer might be the second one. It is much cheaper to find out here than in
Phase 4.

## 3. Tiberius ab initio versus ClaMSA mode, on mammals

**What.** Reproduce the comparison that already exists inside the strongest
current tool. Tiberius ships nine model configurations; exactly one uses ClaMSA
comparative input and it is mammals-only, while all six clade models added in
2026 are sequence-only [gabriel2026accurate].

**Why.** That configuration list is a revealed preference by the group with the
best current architecture, and it admits two readings: comparative input buys
less than the charter assumes once a good sequence model exists, or generating
it per clade is too expensive or too data-hungry to scale. Both are bad news
for a naively comparative design and they point at different fixes.

**Cost.** Code and weights are installed and MIT-licensed; ~1 day, one GPU.
Two caveats from the install bench: Tiberius requires Python ≥3.12 (this
project targets 3.11, so it needs its own pinned environment), and on a
compute-capability-12.0 GPU its pinned TensorFlow must JIT from PTX, which
TensorFlow warns can take 30 minutes or more. Carry `stalin`'s and `engels`'s
qualification with any result: Tiberius's comparative input generator has
supervised human exposure, and its preprocessing cost is separate from
inference.

**This is the highest value-per-hour item in the list.** If the answer is that
comparative input buys little, the charter's central hypothesis needs revision
before anything is built.

## 4. The fixed-decoder encoder ablation

**What.** Hold one structured decoder fixed and swap only the evidence encoder:

- (a) KA/KS window statistic — the floor, already built as `baselines/kaks/`;
- (b) PhyloCSF or ClaMSA codon-model likelihood features;
- (c) a small sequence-only encoder at the target budget (Tiberius-class);
- (d) (c) plus the tree as extra tokens;
- (e) (c) plus the tree as a metric — MDS embedding and Tree-RoPE, the HyphAeon
  transfer.

**Why.** This is the one deliverable nobody else is producing. OrionGeno and
ANNEVO report accuracy; neither isolates the contribution of phylogenetic
context from that of long-range modelling. The specific HyphAeon claim is that
injecting the tree *as a metric* removes the need to learn the phylogeny — and
(d) versus (e) at a fixed parameter budget is exactly that claim, tested.

Three of the five reviews ask for a version of this. `engels` frames it as
tree-likelihood versus distance-aware encoder versus no tree, decoder fixed;
`stalin` as a codon-likelihood feature bank versus the axial design; `lenin` as
three arms with the benchmark as the deliverable.

**Prerequisite from `stalin`, which nobody else raises.** MDS coordinates are
not unique under rotation or sign change. Before (e) is worth building, check
that the encoder is invariant to a random rotation of the embedding, to taxon
reordering, to removed taxa and to rescaled branch lengths. If it is not, (e)
is fitting arbitrary coordinates and the ablation will not mean what it appears
to mean. This is a day's work and it gates the arm.

**Prerequisite from `engels` and `stalin`, jointly.** Measure the candidate
generator's recall ceiling before measuring end-to-end accuracy. A scorer must
be evaluated together with its candidate generator, not credited with a
classifier's conditional accuracy — KA/KS and ClaMSA both classify candidates
whose strand and frame are already given.

## 5. The covariate-conditioned decoder

**What.** `marx`'s specific proposal, and the most concrete architectural idea
in the five reviews: keep the differentiable HMM, but make its length
distributions and GC-dependent emissions **functions of per-genome covariates**
— local GC, the fitild intron-length parameters, the codon usage table —
instead of constants.

**Why it is a real idea and not a detail.** Tiberius's HMM layer fixes geometric
length distributions at empirical mammalian means [gabriel2024tiberius §2.4].
That is a named, specific reason the best current model does not transfer, and
this is the minimal change that removes it: the same weights would serve a
fungus with 60 bp introns and a mammal with 10 kb introns. The covariates it
needs are cheap to compute from the target genome alone, so it does not
re-introduce the alignment dependency.

The same logic covers GC. GENSCAN's 1997 answer was separate parameter sets per
C+G region; GeneMark-ETP's 2024 answer is three GC-specific models per
GC-inhomogeneous genome [brruna2024genemark]. Conditioning on GC rather than
partitioning on it is the modern version of the same fix.

`engels` adds the right caution: compare learned duration conditioning against a
simple heavy-tailed or mixture duration prior before concluding the learned
version is needed, and do not hard-code a familiar mammal's length regime under
another name.

## 6. The fixed-encoder decoder ablation

**What.** The mirror of #4: hold the encoder fixed and compare a differentiable
HMM, a chromosome-scale CRF, and a grammar-constrained splice-graph DAG with
top-*K* path extraction.

**Why.** `trotsky`'s DAG proposal is the most specific decoder design of the
five and its motivation is sound — geometric self-transitions penalize long
introns exponentially, explicit duration submodels cost O(D) or O(D²) and force
hard length cutoffs, and Viterbi collapses alternative splicing to one path.
The published alternative it does not engage with is GeneCAD's CRF, which
reports ~9% transcript F1 over Helixer and BRAKER3 on angiosperms
[liu2025genecad]; and GeneCAD's grammar constraints were added *after* the
neural model, in v0.4.0, which supports the separability the DAG design assumes.

**Why it ranks below #4.** Two reasons. The isoform benefit is currently
unmeasurable — no paper in the review reports transcript-level UTR or isoform
accuracy against a reference with an error bar, and geneML is one of very few
tools with isoform numbers at all (41.1% recall / 71.1% precision against
Iso-Seq [vader2026geneml]). And the DAG adds exactly the failure mode `engels`
and `stalin` warn about: candidate pruning with an unmeasured recall ceiling. Do
#4 first; the same harness runs #6.

## 7. The unsupervised fallback arm

**What.** A Vipsania-style arm: a differentiable HMM inside a masked-language
sequence model, pretrained pan-eukaryotically and fine-tuned unsupervised on the
target genome, never shown a reference annotation [krieg2026vipsania].

**Why it belongs in the shortlist.** It attacks the generalization problem from
the opposite side from the charter — not by giving the model better geometry
but by removing the labelled data whose clade bias causes the degradation — and
it is the strongest published answer for exactly the clades where alignments do
not exist. Its README reports locus F1 of 0.70 in Discoba and 0.66 in Fungi,
clades no other tool reports on at all.

**Why it ranks here and not higher.** The same table shows the residual clade
dependence is not gone: 0.41 in Vertebrata and 0.21 in Arthropoda, and "locus
F1" is not gene F1. `marx` names the decisive comparison: run Vipsania and
Tiberius on the same five held-out species spanning intron regimes and see
whether the gap tracks phylogenetic distance or intron-length divergence. The
former favours the comparative design, the latter favours the covariate design
(#5). Both fit on one consumer GPU.

## 8. The tiny-model probe

**What.** Establish empirically how small the per-task heads can be, starting
from minisplice's result: splice signals conserved across phyla, captured with
**7,026 parameters** [yang2026improving], independently corroborated by TOGA2's
finding that human-trained deep splice models generalize across vertebrates
[malovichko2026accurate].

**Why.** The charter's ~2M-parameter target is the part of its hypothesis that
nobody in the 2026 cohort is testing — OrionGeno needs a modern NVIDIA GPU with
`mamba-ssm`, GeneCAD sits on a foundation model, Evo 2 is 7B–40B parameters.
And if the KA/KS floor and minisplice both do real work at trivial size, the
open question is not whether a small model can work but which part of the
pipeline actually needs the capacity.

**Why it ranks last.** `engels`'s objection stands: shrinking before the
benefit is established optimizes the wrong thing. Treat this as an ablation
axis inside #4 rather than a separate programme, and report **total** cost —
parameters plus alignment construction plus feature generation plus inference —
because a 2M-parameter model that needs a Cactus alignment is not cheap.

---

## What is deliberately not on this list

- **UTR prediction in version 1.** Not because it is hard but because it cannot
  currently be scored: Helixer 2026 strips UTRs before scoring, Tiberius scores
  only the longest CDS per gene, and no reviewed paper reports transcript-level
  UTR accuracy. Predict them if convenient; keep them out of the headline
  metric. See [`disagreements.md`](disagreements.md) §5.5.
- **Replacing evidence integration.** On genomes with RNA-seq, a BRAKER3-style
  pipeline still adds value at the gene ends, which is where every tool is
  weakest (BRAKER3's start-codon F1 is 70% against >87% for splice sites). The
  target is the large majority of assemblies with no annotation and no
  RNA-seq, at laptop cost.
- **A DNA foundation model as the backbone.** GENATATORs reports directly that
  pretrained DNA-LM embeddings do not capture the features needed for gene
  segmentation without task-specific finetuning [shmelev2026genatators], and
  foundation-model scale is outside the charter's compute budget in any case.
  Keep them as baselines.
- **Another architecture, unmeasured.** The reviews' strongest joint finding is
  that the comparative-geometry idea is occupied ground and the measurement is
  not. A fourth phylogeny-aware deep annotator with a new benchmark of its own
  would add nothing the field lacks.
