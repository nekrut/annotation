# Geometric gene prediction: candidate design working draft

Task: [T-human-011](../../tasks/T-human-011.md). Author: stalin.
Written 2026-09-15 UTC against accepted main `09b5386` (claim `9239764`).
Status: first bounded design pass, **in progress**. This is the working
artifact for the eventual `docs/design/proposal.md`, not the submitted
proposal or a Phase 4 implementation. Numerical architecture choices below
are proposed settings; arithmetic is an estimate, not measured performance.

## 1. Recommendation to develop

Develop **B: a small sequence model with optional comparative refinement**,
using **A: the same sequence model and decoder** as its mandatory control.
Keep **C: a splice graph using B's scores** as a later decoder comparison.
This ordering permits a useful whole-genome result with no informants and
isolates the comparative channel before adding candidate-path pruning.
It follows the accepted [encoder/decoder ablation plan][synthesis].

The recommendation is conditional on candidate coverage, geometry checks,
and measured cost. The present arithmetic does **not** establish that B
meets the CPU budget on gene-dense genomes. A fixed compute cap may disable
comparative refinement on such genomes; that must be declared as a distinct
input regime and evaluated, not hidden inside an aggregate.

### Accepted evidence that constrains the design

| Evidence | Consequence |
|---|---|
| [Data inventory][data], sections 2 and 7: usable alignments cover only part of the panel and several resources use different assemblies | DNA-only prediction must be available everywhere; never substitute an older assembly silently. |
| [KA/KS baseline][kaks]: short-exon sensitivity and whole-window precision conflict; its worm sample performs poorly | A fixed KA/KS window must not be a hard gate for any coding prediction. Treat its scores as evidence with an explicit undecidable state. |
| [Benchmark][benchmark], sections 3 and 4: heldout and heldout_paired are different tests; start, stop, phase and complete CDS chains matter | Train only on the train split and decode structures on chromosome coordinates, not independently inside inference windows. |
| [Cost baseline][cost], section 5: 5 M parameter ceiling, 15 CPU-s/Mb on its named runner, 0.5 GPU-s/Mb on RTX 4090 class hardware, 8 GB host and device targets | Apply wide layers after taxon pooling or sequence downsampling; count both strands, every unknown frame and repeated halo tokens. |
| [Geometry disagreement][geometry]: taxon order and MDS axes are arbitrary | Use patristic distances directly as the first tree representation. Coordinate-wise Tree-RoPE requires a separate stability audit. |

The KA/KS results are a seeded gene-window diagnostic, not a held-out
whole-genome accuracy estimate. Some exploratory inputs contain held-out
species, so neither those windows nor their best-performing hyperparameters
become training/development data. Reproduce the frozen baseline on those
same windows for continuity, then report the benchmark evaluation separately.

## 2. Shared input and training contract

**Sequence.** Target nuclear FASTA, assembly accession/checksum, sequence
selection, ambiguity mask, local GC, and a declared nuclear translation
table. Use the same weights on both orientations and map outputs back to
reference coordinates. Start with 4,096-base chunks retaining 3,072 central
bases; the 512-base halos cover the proposed local encoder receptive field.
The chromosome decoder carries state across chunk boundaries. These sizes
are design settings, not limits on gene or intron length.

**Comparative input for B and C.** Initially use the inventory's multiz
sources: mm39/35way, dm6/124way, ce11/135way and sacCer3/7way for training
where source and panel coordinates agree. Human hg38/100way is evaluation
only. Chicken and fugu require an exact-assembly source or a verified,
declared mapping before any comparative score; otherwise run DNA-only on
the panel assembly. Zebrafish remains DNA-only until a compliant alignment
exists on GRCz12ab. Zoonomia stays excluded under the accepted inventory.

Choose up to seven distinct informant species plus the target (K <= 8),
with K = 16 as a cost/accuracy sensitivity experiment. Selection uses
alignment coverage and tree diversity without annotation, a fixed rule
frozen on training development chromosomes. Preserve one deterministic
species selection across a chunk and all compared tree arms. Drop rows
with no observed sequence and mask padding; record actual K per example.
Duplicate assemblies are not independent taxa. No ancestral rows or
multiple paralog copies enter the first comparison. Ambiguous duplicate
mapping is marked unavailable unless a documented syntenic choice exists;
highest identity alone is not an orthology guarantee ([data], section 6.3).

Provide base/gap/unaligned/ambiguous categories, insertion lengths and a
mapping to the supplied tree, including its branch-length units. A gap
inside an alignment and absence of alignment remain different inputs.
Splice-site flanks retain nucleotide resolution. Within candidate exonic
spans, form consecutive genomic triplets in **all three offsets on each
strand**, using weight sharing. These are frame hypotheses, not supplied
CDS phase. Do not concatenate exons using truth to construct model inputs.
Triplets crossing an inferred splice boundary need the decoder's phase
state; until that boundary is known, the nucleotide channel remains valid
and a genomic triplet is only local evidence.

**Extra tracks.** No RNA-seq, protein evidence or shipped conservation track
is required in the first comparison. Motif scores come from the target
sequence and are soft features. An optional later evidence arm must declare
accessions, exclusions and processing cost; it cannot be pooled with the
DNA/comparative result. In particular, removing held-out MAF rows does not
remove their contribution to an existing phyloP or phastCons track.

**Labels and split.** Use the ten `train` species in [panel.tsv][panel];
reserve whole chromosomes from those species for tuning, calibration,
early stopping and candidate thresholds. No held-out annotation, proteins,
projected labels or checkpoint exposure enters fitting. Remove every
held-out species/assembly alias from training multiz rows, including hg38
for mm39 and apiMel4 for dm6, and declare the complete removal manifest.
Jointly inferred alignments require rebuilding, as the benchmark states.
Audit informant membership against all held-out species, not only these
two known row names. Pretraining is initially none.

Sample species with equal weight, with training-derived strata for GC,
exon length and intron regime, plus gene-free/background windows. Generate
candidates from sequence alone over complete development chromosomes;
gene-centered sampling cannot estimate their whole-genome false positives.
Use separate strand labels and real, ordered transcript CDS chains. The
current cutter's union/owner-per-base labels cannot represent every
overlapping transcript. A first chain loss can use its declared
`longest-cds` policy with conflicting same-strand overlaps masked; keep
all real boundary sites for auxiliary losses and count the omitted chains.
An isoform-union mask must never be presented as one valid transcript.
Extending the loader to structured chains is an identified Phase 4 need.

**Covariates.** Local GC and sequence-only composition summaries may
condition a shared decoder. A learned mapping is fitted on train genomes
only; do not calculate a held-out genome's intron distribution or codon
usage from its reference annotation. Start with a pooled duration mixture
and compare sequence-derived conditioning with a taxonomy-conditioned
variant that uses one shared set of weights and an unknown-clade value.
Per-genome parameter fitting is outside this initial training plan.

## 3. Candidate A: compact DNA encoder and duration-mixture CRF

**Encoder.** A narrow nucleotide-resolution convolutional stem and boundary
heads accompany a context path downsampled by 12. The context path has four
local transformer blocks, width 96, with 16 attended positions per token.
It supplies contextual features to the fine-resolution heads; interpolation
does not decide splice coordinates. The four-block context has approximately
442,368 weight-matrix parameters; reserve a total **0.70 M** for the stem,
projections, heads, conditioning and decoder. This reservation is a design
ceiling pending a layer-by-layer inventory, not an instantiated count.

**Encoded structure.** Shared orientation processing avoids separate strand
models. Relative genomic offsets express locality; the nucleotide heads
and decoder expose the three possible reading phases explicitly. A coding
state advances phase only when it emits a CDS base; an intron preserves
it. The model therefore need not discover the modulo-three accounting,
which is distinct from discovering which bases actually code.

**Decoder.** A sparse-transition conditional random field (CRF) has
intergenic, coding-phase, intron-phase and start/stop boundary states.
At intron entry choose a duration component and keep it through the
intron; each component has a geometric tail. A finite mixture avoids one
fixed mammalian mean and leaves positive support for arbitrarily long
introns, but is still asymptotically exponential rather than a true
power-law tail. Compare it with the single-component control. Length
and motif penalties are separate; noncanonical motifs receive finite
scores rather than a zero-probability alphabet mask.

Forward/Viterbi state passes run on the complete selected sequence, using
sparse transitions. Chunk checkpoints and recomputation bound resident
traceback storage; do not retain a genome-sized neural activation tensor.
Restart only at real sequence boundaries. This makes a megabase intron
possible without a megabase attention window or donor-acceptor edge list.
The exact transition/state inventory and traceback resource bound remain
required before the proposal is submitted.

The proposed minimum splice-gap convention follows the scorer's default
20-base distinction between introns and short CDS gaps; record anything
below it separately. It is not a universal biological minimum and must be
configurable for a future benchmark with shorter real introns. Do not
impose a 50-base intron floor or prune short coding exons. The source
inventory documents the recall loss those restrictions cause ([data],
section 9). Complete genes require compatible start/stop and phase;
partial sequence-edge calls are explicitly marked. Translation table 6
uses TAA/TAG as glutamine and TGA as stop ([NCBI genetic codes][codes]);
the table is input metadata, not inferred from held-out CDS labels.
Programmed recoding/frameshifts and same-strand overlapping genes need
explicit exception accounting, not silent forcing into a legal ordinary ORF.

**Fallback and extremes.** A already uses no informants. Short introns are
represented by the short-duration component and fine-resolution boundaries;
long introns retain decoder state across windows. Independent orientations
allow opposite-strand overlaps; a single path per strand cannot represent
every same-strand overlap or alternative transcript.

**Three leading failures.** (1) Narrow, downsampled context loses sequence
signals despite fine-resolution heads. (2) A pooled or conditioned mixture
miscalibrates unseen clades or the long tail, producing fusions/splits.
(3) One chain per strand loses same-strand overlapping loci and isoforms.

## 4. Candidate B: A plus a narrow comparative encoder

**Candidate support.** Cheap sequence scans plus A's frozen DNA scores
nominate coding-like spans and boundary neighborhoods. A still provides
scores at every base; a region omitted from B is predicted by A. Thus the
gate limits comparative improvement rather than deleting all predictions
there. Freeze the support independently of informants/tree for the encoder
comparison. Count CDS-base, exon, donor/acceptor and complete-chain support
recall against development truth, including short/noncanonical cases.
Report both the fraction of genome selected and the actual emitted tokens
after merging, halos, padding, both strands and all frame hypotheses.

**Encoder.** Two axial blocks at width 32 act along local codon positions
(32 attended neighbors) and taxa. Each block has two attention projections
and one expansion-four MLP: approximately 16 d^2 parameters per block,
or 32,768 for the two-block core. Pool along taxa before wider target-side
projections. Reserve **0.90 M total** including A and all extra heads;
parameter reuse across frame hypotheses does not eliminate their FLOPs.
Masked residual comparative scores augment A's emissions; for K = 1
(target only), the comparative residual is identically disabled.
Training includes whole-alignment dropout and taxon/gap degradation, so
missing comparative input is represented during training.

**Tree geometry.** Use supplied patristic distance D(i,j) directly:

`attention_logit(i,j) = dot(q_i, k_j)/sqrt(d_head) + b_head(D(i,j))`.

The scalar distance basis/bias is shared across taxa; there is no species
ID embedding or taxon-index positional encoding. This gives the network
pairwise evolutionary separation without asking it to reconstruct the
tree from bases. The bias does not assert that nearby taxa are independent
observations. Structural attention encoding has a precedent in
[Graphormer][graphormer]; this tree application is a proposed design.

Patristic distance is exactly invariant to arbitrary taxon indexing when
rows, masks, distances and pooling are permuted together. It also avoids
MDS axis choices. For an edge e of length w_e, assign a coordinate sqrt(w_e)
to leaves on one side of that edge and zero to the rest. Then
`||z_i-z_j||^2 = sum(lengths on their path) = D(i,j)`.
This is our explicit construction: square-root tree distances have an
exact Euclidean embedding with enough coordinates; a four-dimensional
projection is generally an approximation. Applying classical MDS to D
and to sqrt(D) are different operations. Report which is used and its
distance distortion if testing a compressed coordinate arm.

Coordinate-wise RoPE does not automatically inherit rotational invariance
of an MDS embedding. [RoFormer][rope] establishes rotary relative-position
encoding, not invariance to arbitrary rotations of tree coordinates.
The charter's [axomeme repository][axomeme] returned HTTP 404 during this
tick, so the actual Tree-RoPE implementation was not audited. Retain it
as a named comparison arm, contingent on an accessible implementation and
either invariant operations or a stable canonical convention. Do not
represent the distance-bias candidate as an implementation of Tree-RoPE.

Caching a distance matrix avoids repeated tree traversal, but the bias
still requires attention-time additions; coordinate rotations also cost
operations. Therefore the [cost baseline][cost]'s statement that a metric
adds no per-token FLOP is not adopted literally. These terms are usually
smaller than dense projections, and still belong in the cost inventory.

**Fallback and extremes.** B becomes A where alignment/refinement is absent.
Unknown/ambiguous alignment supplies no negative coding evidence by default.
Non-coding flanks and fine-resolution boundary scores remain available;
the decoder's phase and duration carry through short and long introns.

**Three leading failures.** (1) Support gating withholds comparative rescue
exactly on weak/short exons. (2) Alignment/paralog errors and sparse training
clades make the residual misleading elsewhere. (3) Dense candidate support
or many informants defeat the CPU budget even with very few parameters.

## 5. Candidate C: B's evidence with a splice graph

Use the identical A/B encoders, inputs, masks and training split. Replace
the chain decoder with a directed acyclic graph on genomic start, stop,
donor and acceptor candidates. Exon/intron edges carry strand, phase,
sequence and duration scores; only compatible paths form transcripts.
An intron is an edge between endpoints, so it can cross chunks without
attending to every intronic base. Permit multiple paths and independent
overlapping loci; use the benchmark's all-emitted-transcript precision
to prevent extra isoforms becoming free recall.

A proposed **1.20 M** total allocation includes B and boundary/edge scoring.
Short introns follow the same configurable gap convention and finite
noncanonical-motif scores as A. Long-distance bins retain candidate edges
through the full chromosome, not just the encoder window. A proposed
8-bin, 4-edge-per-bin pruning rule gives at most 32 retained successor
edges per candidate, but it is an approximation: report complete-reference
path survival, including the length-tail stratum, before scoring the model.
No claim of exact all-path decoding survives this pruning.

**Fallback.** With no informants the graph uses A's DNA evidence. Candidate
graph failure can return A's valid chain output, with a declared fallback
counter. Structural recall then remains limited by A in those regions.

**Three leading failures.** (1) Boundary or edge pruning removes the true
path. (2) Long-distance pairing fuses neighboring genes or expands the
graph on repeats. (3) Top-path output overpredicts isoforms and still has
insufficient labels for overlapping loci. C is third because these risks
are additional to the encoder question and its decoder cost is unbounded
until candidate counts are measured.

## 6. Compute arithmetic and expected operating range

The companion [budget-arithmetic.tsv](budget-arithmetic.tsv) is reproduced
by the formula below using Python standard-library arithmetic. It counts
matrix projections and attention products only. Values divided by the
[cost baseline][cost]'s assumed 3e10 CPU FLOP/s and 1e13 GPU FLOP/s are
**conditional arithmetic times**, not measured runtimes or guarantees.
Packing, stem/heads, nonlinearities, masking, I/O, tree setup, decoder and
device synchronization are omitted and must be added in a complete budget.

Let N = 1,000,000 bases, f be pre-halo comparative support fraction, K total
rows including target, and use a provisional comparative halo/padding factor
1.2. All are scenario assumptions; f is not the genome's annotated CDS
fraction. Dense coding genomes are why f = 0.75 and 1 are included.

```python
nA = 2 * N * (4096 / 3072) / 12       # strands, halo repeats, downsampling
pA = 4 * 12 * 96**2
flopA = 2 * pA * nA + 4 * 4 * nA * 16 * 96
nB = 6 * N * f * K * 1.2 / 3         # every strand/frame, not known CDS
pB = 2 * 16 * 32**2
flopB = 2 * pB * nB + 2 * 4 * nB * 32 * (32 + K)
```

| Candidate/scenario | Counted GFLOP/Mb | CPU s/Mb at assumed rate | GPU s/Mb at assumed rate |
|---|---:|---:|---:|
| A core | 202.069 | 6.74 | 0.020 |
| B, K=8, f=0.05 | 274.814 | 9.16 | 0.027 |
| B, K=8, f=0.75 | 1293.244 | 43.11 | 0.129 |
| B, K=8, f=1 | 1656.969 | 55.23 | 0.166 |
| B, K=16, f=0.05 | 351.491 | 11.72 | 0.035 |
| B, K=16, f=0.75 | 2443.401 | 81.45 | 0.244 |
| B, K=16, f=1 | 3190.511 | 106.35 | 0.319 |

A's initial engineering target is <=15 CPU-s/Mb and <=0.5 GPU-s/Mb.
For B, the same targets are plausible only conditionally: sparse support
leaves some CPU headroom, while dense support already exceeds the CPU
ceiling in counted operations alone at the assumed sustained rate. On a
GPU the arithmetic leaves headroom but does not predict small-kernel and
decoder efficiency. C inherits B's encoder time **plus an as-yet unpriced
graph/traceback cost**; assigning it a precise seconds/Mb forecast now
would invent a measurement. A two-pass strand computation shares weights,
not execution time. All three must include preprocessing other than the
separately displayed alignment construction bill in these targets.

Report cold setup and warm prepared-input inference separately, then the
end-to-end total including retrieval/construction of alignments as required
by [benchmark][benchmark] section 4.7. A cached public alignment is not a
zero-cost input for a newly sequenced species. Memory targets are 8 GB
host and 8 GB device for every candidate; tiling makes these plausible but
does not establish a peak. CPU-only mode and missing-alignment mode are
separate declared runs, using the same frozen model weights.

## 7. Experiments specified for Phase 4

1. Freeze data, candidate policy, decoder and training budget. Compare
   DNA-only; same MSA/no tree; same MSA/tree tokens; same MSA/patristic bias;
   and an audited MDS/Tree-RoPE arm. Keep informants, masks, support,
   decoder, split and effective capacity controlled. Tree-based selection
   must be fixed across arms; no-tree means no tree supplied to the encoder,
   conditional on the same selected alignment, not absence of phylogeny
   from upstream alignment construction. Report token-arm overhead.
2. Compare the KA/KS score and codon-likelihood/ClaMSA features under the
   same structure decoder as separately declared controls. Do not reuse
   a supervised feature generator exposed to held-out labels as if it
   were an untrained statistic. Charge feature generation explicitly.
3. Check taxon permutation and MDS rotation/sign representation symmetries
   at emission-score level before comparing final labels. Measure, rather
   than constrain to equality, changes after removing taxa or rescaling
   branch lengths. Verify K=1 fallback, all-masked attention behavior,
   reverse-complement mapping and shifted chunk seams.
4. Evaluate true-path support before model accuracy for C, and comparative
   support recall before attributing B's misses to its encoder. Report
   A/B results on the entire genome and on aligned/unaligned subsets with
   the same complete-genome denominator retained. No new scorer definition
   is needed; the input regime belongs in run metadata and report grouping.
5. Report the benchmark's exact transcript, exon, splice, start/stop and
   locus metrics, fusions/splits and GC/intron/exon strata. Keep paired
   versus cross-clade aggregates separate; Tetrahymena is reported and
   excluded from ranking as specified. Select checkpoints only on train
   development chromosomes. Repeat-seed uncertainty and per-species
   changes precede any claim that geometry wins.

## 8. Work remaining before a review PR

- Finish the decoder transition/phase convention, duration parameterization,
  handling of partial ends and resource bounds for checkpointed traceback.
  Quantify which overlapping transcripts a single-path control cannot emit.
- Replace parameter reservations with an explicit layer inventory and
  include stem/heads/decoder in the CPU/GPU budget. Specify an observable
  candidate/edge-density experiment to price C without prototype training.
- Decide the exact non-neural support scan, density cap/fallback policy and
  support-recall criteria on train development data. Resolve how selected
  short-exon halos change the provisional 1.2 token multiplier.
- Turn this working artifact into `docs/design/proposal.md` on
  `work/T-human-011-stalin`, open a PR, request reviews from at least two
  other agents through relay, and send the final ranked proposal to human
  for the required Phase 4 decision. No model training was run this tick.

## Sources

Accepted repository documents above are pinned to this draft's main commit.
External sources were read through public pages on 2026-09-15 UTC; only
metadata and our own design notes are stored here.

[synthesis]: ../../../docs/review/candidates.md
[geometry]: ../../../docs/review/disagreements.md#54-is-tree-as-metric-mathematically-well-posed
[data]: ../../../docs/data-sources.md
[benchmark]: ../../../docs/benchmark.md
[panel]: ../../../benchmark/panel.tsv
[cost]: ../../../docs/cost-baseline.md
[kaks]: ../../../baselines/kaks/README.md
[codes]: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG6
[graphormer]: https://arxiv.org/abs/2106.05234v5
[rope]: https://arxiv.org/abs/2104.09864v5
[axomeme]: https://github.com/nekrut/axomeme
