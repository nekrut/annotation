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
| 3 | Tiberius ab initio vs ClaMSA mode, replicated on mammals | ~1 day, 1 GPU **on prepared ClaMSA features**; a preprocessing bill if they must be regenerated | how much *this supplied comparative pipeline* buys on *these mammalian regions* — prioritization evidence for #4, not a verdict on comparative input in general | lenin |
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
back in a new form. The alignments the *reviewed literature* used are mammal-
and bird-centric, but the accepted data inventory (`docs/data-sources.md`,
T-human-008, `done`) records three deep non-vertebrate alignments UCSC serves
directly: `dm6/multiz124way`, `ce11/multiz135way` and `sacCer3/multiz7way`.
So insects, nematodes and yeasts are a download, and the degradation curve can
be measured on them this quarter at no alignment-construction cost. What is
genuinely missing is plants, most fungi and the protists — and, for every
covered clade, informants away from the one reference assembly the alignment
is built on. Two cases must be kept apart here. A genome that *is* in the
alignment but is not the one it was built around is reachable: HAL documents
queries relative to an arbitrary reference or subtree and MAF export with a
selectable reference ([README](https://github.com/ComparativeGenomicsToolkit/hal#readme)),
so the cost is a resource-specific extraction, not a missing capability — and
the multiz MAF files UCSC serves are reference-anchored, so they need that
conversion or a re-projection first. A genome that is *absent* from the source
alignment is the real gap: no change of reference recovers sequence the
alignment never contained, and a new assembly must be aligned in. Neither case
is evidence about held-out-species accuracy on its own.

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

**Cost.** Code and weights are installed and MIT-licensed; ~1 day, one GPU —
*conditional on the published ClaMSA features for the evaluated regions being
downloadable and compatible with the released checkpoints*. That is the whole
estimate: it is inference on prepared input. If the features have to be
regenerated, the comparison acquires a preprocessing bill that this number does
not contain and that must be reported as a separate line; the only figure the
reviews have for the historical shared feature-generation workload is 492 CPU
node-days with the node core count unstated, which is neither a per-genome
charge nor a measurement anyone here has made
([stalin's supplement audit](../../relay/messages/20260909T084030Z-stalin-0008.md)).
Record the checkpoint identity and the exact species and regions evaluated with
any result. Two further caveats from the install bench: Tiberius requires
Python ≥3.12 (this project targets 3.11, so it needs its own pinned
environment), and on a compute-capability-12.0 GPU its pinned TensorFlow must
JIT from PTX, which TensorFlow warns can take 30 minutes or more.

**What it cannot decide.** The measured quantity is the benefit of *this
supplied comparative pipeline* as a whole. It does not separate tree geometry
from alignment information, and it does not separate either from the
feature generator's supervised human exposure (Methods 3.7: the sitewise
ClaMSA input generator was trained with human chromosome 17 RefSeq labels).
Candidate #4, the fixed-decoder encoder ablation, is still required for the
causal claim.

**What its result licenses.** A gain measured here is a statement about *this
pipeline, these checkpoints and these mammalian regions*, and it is
prioritization evidence, not a gate. A large gain raises the priority of #4; a
small one lowers it, and says that this particular supplied encoder adds little
on this material. Neither outcome decides whether a different comparative
encoder, a different informant set, or a different clade would help, and a
small gain here is not grounds for revising the charter's central hypothesis.
Attribution to the tree, or to the alignment, requires #4 with the decoder and
inputs held fixed — and even then the claim is scoped to the encoders,
informants, species and split #4 evaluates, not to comparative encoders or
clades in general. `engels` and `stalin` both asked for this scope limit
([20260910T042437Z-engels-0027](../../relay/messages/20260910T042437Z-engels-0027.md),
[20260910T043806Z-stalin-0028](../../relay/messages/20260910T043806Z-stalin-0028.md)).

**It is still the highest value-per-hour item in the list**, because it is a
day of inference that reorders the rest of the list.

## 4. The fixed-decoder encoder ablation

**What.** Hold one structured decoder fixed and swap only the evidence encoder:

- (a) KA/KS window statistic over the MSA — the floor, already built as
  `baselines/kaks/`;
- (b) PhyloCSF or ClaMSA codon-model likelihood features over the same MSA;
- (c) **target DNA only**, a small sequence-only encoder at the target budget
  (Tiberius-class) — the no-alignment baseline;
- (d) **MSA, no tree**: (c)'s encoder plus the aligned informant rows, with no
  phylogeny supplied in any form;
- (e) **MSA + tree as tokens**: (d) plus the tree as extra tokens;
- (f) **MSA + tree as a metric**: (d) plus MDS embedding and Tree-RoPE, the
  HyphAeon transfer.

**What each contrast isolates, and what must be held fixed.** (c) vs (d)
measures the value of the alignment; (d) vs (e) and (d) vs (f) measure the
value of the tree *given* the alignment; (e) vs (f) is the HyphAeon claim
itself. That reading only holds if (d), (e) and (f) consume **the same aligned
informants, the same masking, the same candidate support, the same fixed
decoder, the same training and evaluation split, and a controlled parameter
budget** — the tree representation is then the only thing that varies. Without
the (d) arm, "sequence-only versus plus the tree" changes alignment content and
tree content at once and measures neither. `engels` and `stalin` both raised
this
([20260910T042437Z-engels-0027](../../relay/messages/20260910T042437Z-engels-0027.md),
[20260910T043806Z-stalin-0028](../../relay/messages/20260910T043806Z-stalin-0028.md));
specifying the arms is a Phase 3 writing task, not a Phase 4 run.

**Why.** This is the one deliverable nobody else is producing. OrionGeno and
ANNEVO report accuracy; neither isolates the contribution of phylogenetic
context from that of long-range modelling. The specific HyphAeon claim is that
injecting the tree *as a metric* removes the need to learn the phylogeny — and
(e) versus (f) at a fixed parameter budget, over identical alignments, is
exactly that claim, tested.

Three of the five reviews ask for a version of this. `engels` frames it as
tree-likelihood versus distance-aware encoder versus no tree, decoder fixed;
`stalin` as a codon-likelihood feature bank versus the axial design; `lenin` as
three arms with the benchmark as the deliverable.

**Prerequisite from `stalin`, which nobody else raises.** Two different checks,
and an earlier draft of this file collapsed them into one gate; `stalin`'s
opinion asks for *stability tests*, not for a blanket invariance requirement
([review](../../relay/artifacts/T-human-004/review.md), and
[20260910T033840Z-stalin-0027](../../relay/messages/20260910T033840Z-stalin-0027.md)).

1. **Representation symmetries — require invariance or equivariance.** MDS
   coordinates are not unique under rotation or sign change, and taxon order
   is an arbitrary indexing choice. Under a consistent permutation of the
   taxon axis, or an arbitrary orthogonal transformation of the embedding,
   nothing about the evidence has changed, so the prediction should not
   change either. If it does, (f) is fitting arbitrary coordinates and the
   ablation will not mean what it appears to mean. The alternative to an
   invariant operation is a stable canonical-coordinate convention, which is
   the other option `stalin` names.
2. **Changes in evidence — measure degradation, do not require equality.**
   Removing an informant removes evidence; rescaling branch lengths changes
   the modeled evolutionary distances. Predictions *may legitimately* move in
   both cases, so an equality test is the wrong instrument: passing it is not
   evidence that the encoder uses its comparative channel, and failing it is
   not evidence that it does. Nor does the converse hold: predictions that stay
   put on redundant or uninformative evidence do not by themselves show that an
   encoder ignores its comparative input, and unchanged decoded labels do not
   even imply unchanged scores. Neither universal equality nor
   universal change is the test. What to report is robustness, calibration and the shape
   of the degradation curve against the number and distance of the informants
   retained.

§5.4 of [`disagreements.md`](disagreements.md) states this the same way.
Check 1 is a day's work and it gates the arm; check 2 is an evaluation the
arm produces, not a gate on building it.

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
