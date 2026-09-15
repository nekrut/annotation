# Geometric gene prediction: candidate design working draft

Task: [T-human-011](../../tasks/T-human-011.md). Author: stalin.
Written 2026-09-15 UTC against accepted main `09b5386` (claim `9239764`);
decoder/resource revision against main `86abfc5`, then layer/tile accounting
against main `f424eb7` later that day.
Status: third bounded design pass, **in progress**. This is the working
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
reference coordinates. Use 4,104-base chunks retaining 3,072 central
bases, with 516-base halos. All three lengths are divisible by 12 so the
context-pooling grid stays fixed on chromosome coordinates (section 3.5).
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

**Encoder.** A nucleotide-resolution stem of width 16 accompanies a context
path downsampled by 12. The context path has four local transformer blocks,
width 96, with 16 attended positions per token. Repeat each contextual
vector over its 12 target positions and fuse it with the fine-resolution
stem before the emission head; this coarse vector does not decide splice
coordinates. The explicit inventory in section 3.5 gives **455,841 total
trainable scalars**, including biases, normalization and the pooled decoder.
This is a count of a specified design, not an instantiated model. Covariate
conditioning beyond the input GC channel is a separately budgeted ablation.

**Encoded structure.** Shared orientation processing avoids separate strand
models. Relative genomic offsets express locality; the nucleotide heads
and decoder expose the three possible reading phases explicitly. A coding
state advances phase only when it emits a CDS base; an intron preserves
it. The model therefore need not discover the modulo-three accounting,
which is distinct from discovering which bases actually code.

**Decoder.** The conditional random field (CRF) below retains unfinished
codon bases as well as reading phase. Phase alone cannot recognize a stop
codon split by an intron. Use a sparse, chromosome-wide recurrence with
shifted geometric duration components, finite noncanonical-motif scores,
and checkpointed traceback. The minimum-length prefix can be evaluated
with delayed entry rather than a separate update for every prefix state
at every base. Sections 3.1 to 3.4 specify the proposed recurrence and its
resource accounting; this is a design, not an implemented decoder.

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

### 3.1 Coding state and boundary contract

Work in oriented, zero-based sequence coordinates. A state at boundary t
describes bases strictly before t; a usual transition consumes base x[t].
The prefix q is the unfinished codon in the **spliced** CDS, not the last
genomic bases. Its length p is the number of CDS bases emitted modulo
three. Maintain these logical state families:

| State | Meaning | Maximum allocation, derived from the four-base alphabet |
|---|---|---:|
| U | Outside a coding chain | 1 |
| E(q) | Ordinary CDS after a completed initiator; q has length 0, 1 or 2 | 1 + 4 + 16 = 21 |
| S(q) | Initiator not yet complete; q has length 1 or 2 | 4 + 16 = 20 |
| I(c, r) | Intron tail holding the complete coding state c in E or S and duration component r | P R, with P = 41 as the unpruned allocation |

Unreachable start prefixes have score minus infinity. Using the supplied
translation table to remove them is a possible exact optimization, not
assumed in the resource bounds. Each consuming coding transition appends
x[t] to q. At length three, S accepts only a declared initiator and moves
to E(empty); E returns to E(empty) for a sense codon or terminates the
chain in U for a stop. There is no continuation through a stop as ordinary
CDS. U either consumes an intergenic base or starts S with the first
initiator base. The start score is attached to that first base; the stop
score is attached to completion of the terminal codon. Split initiators
and terminal stops therefore use the same prefix machinery.

From an E/S state c, a donor can open an intron without changing c.
Closing that intron consumes the **first following CDS base**, applies
its ordinary coding transition, and adds the acceptor score. This ensures
at least one CDS base between interior introns without a separate
zero-length-exon state. The same c, including the nucleotide identities
in q, survives any number of intronic bases. A graph decoder for C needs
this compatibility information too; matching only phases is insufficient.

The required score channels per oriented base are: U (one), CDS by p
(three), intron by p (three), and start/stop/donor/acceptor (four): **11**
channels, by this inventory. Prefix states share the phase emission;
they do not require separate neural output channels. Motifs add finite
boundary scores and never prohibit noncanonical junctions. GFF3 CDS
phase for a segment beginning after p emitted CDS bases is `(3-p) % 3`;
derive it along the traced chain, with strand-aware coordinates. Emit
terminal-stop bases inside CDS and explicit codon features, including
split features when needed, under the [benchmark][benchmark] convention.

The declared table supplies the coding alphabet and allowed initiators,
using the [NCBI genetic codes][codes]. In table 6, TAA/TAG are sense and
TGA is stop; table 1 treats all three as stops. For ambiguous observed
bases, use a normalized observation prior over the permitted bases
(initially uniform), sum the weighted choices in the forward recurrence
and maximize them in Viterbi. This avoids rewarding coding states merely
for having more compatible latent bases. Mark any affected emitted CDS as sequence-uncertain;
a compatible latent base assignment is not evidence of an intact ORF.
No held-out CDS annotation is used to resolve ambiguity or choose a table.

At an actual sequence edge, allow E/I entry or E/S/I exit with an explicit
partial-end prior; recover the partial flag from traceback. An entry in
I is a residual intron, with no invented upstream donor. Discard paths
with no observed CDS. Initialize unfinished ordinary codons with a
normalized prior over phase and prefix so phases with more prefixes do
not gain probability merely by multiplicity. Interior chunk boundaries
have no partial-entry/exit option. Recoding, programmed frameshifts and
same-strand overlapping chains remain declared exceptions to this grammar.

### 3.2 Duration recurrence and exact delayed entry

Proposed settings are minimum intron length m = 20 and R = 3 components;
these are design choices, not biological estimates. The benchmark counts
splice gaps of at least 20 bases ([section 4.3][benchmark]); changing that
scoring convention is separate from changing the decoder's configurable
minimum. For phase p, weights pi[p,r] sum to one and 0 < q[p,r] < 1:

`Pr(length = l | p) = sum_r pi[p,r] * (1-q[p,r]) * q[p,r]**(l-m), l >= m`.

The first comparison learns pooled weights and hazards on train data.
They stay constant through the intron. Local GC can condition donor and
acceptor scores now; conditioning the duration law itself is a separate
arm whose parameters must be fixed at entry or constant for the sequence.
Recomputing a hazard from local GC at each intronic base would define a
different, position-dependent duration model. The mixture supports every
length above m, with an exponential asymptotic tail; it is not a power law.

A direct expansion has m-1 mandatory intron states per coding state,
followed by R tail states. An equivalent delayed-entry recurrence uses
only the tails as active intron states. Let D[t,c] be the coding-state
score, d[t] the donor score at boundary t, and u[p,t] the intron emission.
The score for entering tail r at boundary t is

`D[t-m,c] + d[t-m] + sum(u[p,j] for j in range(t-m,t)) + log(pi[p,r])`.

Combine it with continuation from `I[t-1,c,r] + u[p,t-1] + log(q[p,r])`.
At exit, add `log(1-q[p,r])`, the acceptor score at t, and the normal
coding transition consuming x[t]. At l = m no continuation factor is
paid; at l = m+1 exactly one is paid. A ring buffer of the last m donor
entries and rolling phase-emission sums evaluates the fixed-length jump
without enumerating donor/acceptor pairs or visiting m prefix states per
base. It forbids lengths below m and permits arbitrarily long tails.

This equivalence assumes the mandatory prefix uses the same per-phase
intron emission as the tail and has no age-dependent neural features.
More elaborate short-intron shapes would need another derivation and
budget. Keep the expanded recurrence as the specification oracle for
Phase 4 checks of the optimized form.

The CRF forward pass sums over latent duration components. For a known
training chain its numerator also sums those components; masked or
ambiguous labels require a constrained numerator, not an invented union
transcript. Viterbi chooses the best **joint** chain/component assignment.
It is not marginal-MAP over chains after summing component assignments.
That inference choice must be identical in the A/B encoder comparison.

### 3.3 Checkpoints, scratch and decoding work

The following bounds are our arithmetic for this specification, not
measurements. Use P = 41, R = 3, m = 20, a traceback block of B = 65,536
bases and an illustrative chromosome N = 250,000,000 bases. Both strands
run sequentially. For unambiguous input, each coding state has at most
one coding successor. The expanded state/transition bounds are
`S = 1 + P*(m+R)` and `E = 2 + P*(m+3*R)`; delayed entry gives
`S = 1 + P*(1+R)` and `E = 2 + P*(1+3*R)` transition candidates per base.
The latter includes P R jump entries; donor-buffer updates, emission
arithmetic and output serialization are additional work.

| Quantity, using those assumptions | Expanded minimum-length chain | Delayed entry |
|---|---:|---:|
| Active scores per boundary | 944 | 165 |
| Float64 slots reserved per checkpoint | 944 | 1,069 |
| Transition candidates/base, unambiguous input | <=1,191 | <=412 |
| Candidates/Mb, both strands, forward plus traceback replay | <=4.764 billion | <=1.648 billion |
| Local uint16 traceback storage, B+m positions | 123,769,728 bytes | 21,633,480 bytes |
| All checkpoints for illustrative N, one strand | 28,818,432 bytes | 32,634,432 bytes |

The delayed-entry checkpoint reserves `S + m*P + 4*(m+1)` doubles:
active scores, pending donor scores, phase-emission history/sums and
normalization-offset allowance. The checkpoint count is
`ceil(N/B)+1 = 3,816`. Use `8*slots*count` bytes for checkpoint scores
and `2*S*(B+m)` bytes for local traceback. Float64 rolling scores are
additional. A uint16 incoming-edge code suffices for this specified
transition inventory, including the four possible latent base choices;
the jump code recovers donor boundary t-m and c. A traceback jump may
cross a block boundary: resume from its true predecessor using the prior
checkpoint, not from a fabricated boundary state. Keep numerical offsets
for delayed scores consistent with active-state normalization.

Thus the decoder arrays alone fit comfortably within the host-memory
target in this example. This does **not** price neural activations, input
packing or library overhead. It also does not establish runtime: a
transition reduction is not a dense-matrix FLOP, and a serial chromosome
recurrence may dominate GPU wall time. Implementation and measurement
remain Phase 4 work. Ambiguous bases can add coding-successor candidates;
the canonical-input transition figures are not all-input upper bounds.

Traceback replay requires the same emissions again. There are two
explicit resource regimes:

- **Spool emissions.** Feed forward scores to the decoder while writing
  the 11 float32 channels to scratch; replay reads them. This is 44 bytes
  per oriented base, or 11 GB scratch for illustrative N when processing
  strands sequentially. Writing and reading both strands transfers
  176 MB/Mb, excluding FASTA/alignment I/O. Scratch need not be resident,
  and none of it would be committed to git. Disk bandwidth/time belongs
  in the inference result; limit scratch by sequence and remove it after
  traceback. This retains the one-neural-pass arithmetic in section 6.
- **Regenerate emissions.** Replay the same frozen chunks, masks, support
  and kernels, with dropout disabled and deterministic tie handling.
  This approximately doubles the counted neural work as well as replaying
  the decoder. The companion [TSV](budget-arithmetic.tsv) includes those
  scenarios. With the revised layer/tile inventory, A's counted neural
  work alone rises to 15.43 conditional CPU-s/Mb; B at K=8 and f=0.05
  rises to 21.81, before omitted operations. Regeneration therefore
  exceeds the CPU target even for A under the assumed throughput.
  It cannot be silently treated as free.

Adopt emission spooling as the proposed first implementation regime,
with scratch usage reported alongside RAM. Preserve emission precision
between forward and replay; lossy caching must be declared and evaluated.
Streaming checkpoint files instead of retaining every checkpoint can
further bound RAM for longer chromosomes, at an additional I/O cost.
The memory argument concerns inference, not training through a whole
chromosome. Training can use truth-derived boundary conditions within
train-only blocks; inference must retain unrestricted chromosome state.

### 3.4 Hand-checkable acceptance cases for Phase 4

These are consequences of the specified grammar and required future
implementation checks, not results from an implemented predictor. In
the examples, `|intron|` consumes at least m genomic bases and none of
the displayed CDS bases; concatenate the displayed exonic strings to
verify the codons.

| Example | Required behavior |
|---|---|
| `ATGAAAT |intron| AA` | Table 1 complete CDS `ATG AAA TAA`; retain prefix T across a phase-1 intron. |
| `ATGAAATA |intron| A` | The same stop split after TA; retain the phase-2 prefix. |
| `A |intron| TGAAATAA` | Complete the split initiator through S states; do not require genomic adjacency of its bases. |
| `ATGAAAT |intron| AATGA` | Table 6 permits `ATG AAA TAA TGA`; table 1 must terminate at the earlier TAA. |
| Lengths m-1, m and m+1 | Reject the first; assign pi*(1-q) and pi*q*(1-q) per component to the other two before neural scores. |
| An intron spanning a checkpoint or an encoder chunk seam | Preserve prefix, duration component and score; no new start/end permission at the seam. |

Also require exhaustive tiny-lattice forward/Viterbi agreement between
expanded and delayed entry, brute-force agreement of the CRF partition,
reverse-complement coordinate/phase agreement, ambiguous-base handling,
and shifted-chunk replay checks once Phase 4 is authorized. Quantify the
single-path label ceiling on train development chromosomes by comparing
the retained real chains with all reference chains and their overlapping
spans. No ceiling percentage is asserted without that count.

### 3.5 Explicit layer inventory and receptive field

All dimensions and counts here are proposed settings; the formulas are
the source of the counts. A uses eight scalar channels per base: four
ACGT indicators, ambiguity, soft masking, local GC and real-sequence
availability. GC is the fraction of G/C among unambiguous bases in a
centered 129-base window, with a fixed zero value if none exist. Ambiguity
and availability remain separate channels. These features need no
reference annotation. There are no species embeddings or pretrained
parameters in the first model.

After a kernel-nine convolution from eight to 16 channels, use three
residual blocks, each with one layer normalization, a kernel-nine
depthwise convolution and a 16-to-16 pointwise convolution. Their
dilations are 1, 2 and 4. Use GELU in the residual branches. Mean-pool
groups of 12 stem positions, project to width 96, and apply four
pre-normalized attention/MLP blocks with four heads, attention offsets
-8 through +7, and expansion-four MLPs. Each block has its own learned
relative-offset bias. Repeat the context vector at nucleotide resolution;
concatenate it with the 16-channel stem, project 112 to 32 with GELU,
then 32 to the 11 emission channels. Evaluate these last heads only on
the retained central positions.

| A component | Trainable-scalar formula | Count |
|---|---|---:|
| Stem convolution, including bias | `8*16*9 + 16` | 1,168 |
| Three depthwise/pointwise residual blocks, biases and norms | `3*(16*9 + 16 + 16*16 + 16 + 2*16)` | 1,392 |
| Pooled-context projection and bias | `16*96 + 96` | 1,632 |
| Four attention/MLP blocks, biases, two norms and offset biases | `4*(12*96**2 + 9*96 + 4*96 + 4*16)` | 447,616 |
| Fine/context fusion and bias | `112*32 + 32` | 3,616 |
| Emission projection and bias | `32*11 + 11` | 363 |
| Pooled decoder | `3*3*2 + 2*16 + 4` | 54 |
| **A total** | Sum above | **455,841** |

The decoder count includes phase/component mixture logits and hazard
logits, donor/acceptor dinucleotide score tables, and four partial-entry/
exit family scores (coding versus intron). All start/stop compatibility
and prefix transitions are fixed by the declared genetic code. Learned
covariate-to-duration maps and taxonomy-conditioned variants are absent
from this initial inventory, and must add their own parameters and work.

The stem has nucleotide radius `4 + 4*(1+2+4) = 32`; the GC channel
adds at most 64. Pooling and repetition add at most 11, and four context
blocks add at most `4*8*12 = 384`. The maximum dependency radius is thus
**491 bases**, below the proposed 516-base halo. A full 4,104-base chunk
has exactly **342 context tokens**. Core starts and pooling groups are
anchored to the oriented chromosome origin; masked edge padding adds no
new sequence. An alternative downsampling grid is a distinct computation,
so a chunk-seam check must hold that grid fixed. No full-chunk statistics,
global attention or normalization over spatial positions is allowed in
this radius argument; normalization is over channels at each position.

## 4. Candidate B: A plus a narrow comparative encoder

**Candidate support.** Cheap sequence scans plus A's frozen DNA scores
nominate coding-like spans and boundary neighborhoods. A still provides
scores at every base; a region omitted from B is predicted by A. Thus the
gate limits comparative improvement rather than deleting all predictions
there. Freeze the support independently of informants/tree for the encoder
comparison. Count CDS-base, exon, donor/acceptor and complete-chain support
recall against development truth, including short/noncanonical cases.
Report both the raw selected intervals and the actual emitted tokens
after tile expansion, halos, padding, both strands and all frame hypotheses.
Section 4.1 fixes the accounting geometry; the seed scan and its thresholds
still require a training-development density/recall audit.

**Encoder.** Two axial blocks at width 32 act along local codon positions
(32 attended neighbors) and taxa. Each block has separate site and taxon
attention modules and one expansion-four MLP: 16 d^2 matrix parameters
per block. Each triplet input concatenates its three bases' eight
categories (ACGT, ambiguous, gap, unaligned, padding) and one `log1p`
insertion-length value per base, for 27 inputs projected to 32. The two
axial blocks use four heads, three pre-norms per block and GELU in the MLP.
After taxa mixing, gather the designated target row, with the designation
permuted alongside rows in symmetry checks. At each retained target base,
concatenate the three covering frame-hypothesis vectors with A's local
stem (3*32 + 16 = 112), then project 112 to 32 to 11 residual channels.
The all-frame gather preserves distinct hypotheses without using truth
phase. Section 4.1 counts **39,180 additional scalars**, or **495,021
total including A**. Parameter sharing does not eliminate repeated FLOPs.
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
The charter's [axomeme repository][axomeme] returned HTTP 404 during the initial
design pass, so the actual Tree-RoPE implementation was not audited. Retain it
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

### 4.1 Comparative inventory and an observable token meter

Each axial block has site offsets -16 through +15 and a learned
four-head offset-bias table. For taxon attention, use eight fixed linear
hat basis functions with evenly spaced knots on `D/(1+D)` in [0,1], with
four learned coefficient vectors per block. Supplied distances must be
nonnegative with documented units; a change of branch-length scale is
a change of input evidence. Precompute the bias per selected tree/K and
charge its attention-time additions separately from matrix products.

| B addition | Trainable-scalar formula | Count |
|---|---|---:|
| Triplet projection and bias | `27*32 + 32` | 896 |
| Two axial blocks: matrices, biases, three norms, site and tree biases | `2*(16*32**2 + 13*32 + 6*32 + 4*32 + 4*8)` | 34,304 |
| Frame/fine fusion and bias | `112*32 + 32` | 3,616 |
| Residual emission head and bias | `32*11 + 11` | 363 |
| Shared scalar residual scale | `1` | 1 |
| **B additional / total with A** | Sum above / plus section 3.5 | **39,180 / 495,021** |

The scale is bounded with tanh. The residual is zero at bases with no
observed informant base, as well as for K=1; gap-only or unaligned rows
therefore cannot create negative evidence by themselves. An observed
base includes an explicitly ambiguous observation; its category remains
available to the network. Degradation/dropout must recompute this mask.
Padding is always excluded from taxon attention. All-masked cases bypass
the comparative module rather than normalizing an empty attention set.

Replace the old provisional 1.2 halo multiplier with this explicit first
schedule. Partition each oriented chromosome into **384-base cores**,
anchored at that orientation's coordinate zero. Expand its nominated
intervals to the cores they intersect, and score all three frame
hypotheses on every selected core. Process the strands independently;
support in one orientation does not wait for the other's DNA scores.
Adjacent cores remain separate in the initial meter; do not credit hypothetical
halo reuse. Each core gets 102 flanking bases per side, hence a 588-base
nominal span, and each frame hypothesis allocates **196 triplet tokens**.
Offsets one and two need at most two additional right-flank input bases;
fetch them but do not add tokens. Reverse orientation fetches the mirrored
flank. Two site-attention blocks reach at most 32 triplets, or 96 bases,
and token extent adds at most two bases, within the 102-base halo.
The nucleotide stem fused at the target base is A's already-counted stem.

For S_plus and S_minus selected cores and constant K, the exact allocated
comparative token count is `3*(S_plus+S_minus)*196*K`. Cores have unique
central coordinates within each orientation, so
each has one owner of its residual outputs; its neighbors' halos never
emit duplicate predictions. A 3,072-base A core contains exactly eight
such B cores. For variable K, sum the expression per core using allocated
rows, including any batch padding. Report useful observed rows separately.
Clip output coordinates at sequence edges, but charge padded tokens.

Use `f = 384*(S_plus+S_minus)/(2*N)` for **mean allocated support across
orientations** in the scenarios in section 6. Then the halo factor is
`588/384 = 1.53125`, before the extra
raw input bases needed for the other offsets. This replaces the previous
1.2 assumption; it is not a measured data property. Distinguish raw seed
coverage, real bases in selected cores, allocated central bases and token
count, reporting each orientation and the physical-coordinate union
separately. On a genome of many short sequences the allocated fraction can
exceed one. For example, a one-base isolated seed occupies a whole
384-base core in this schedule; annotated CDS fraction cannot predict
that allocation. Thresholds selected before core expansion would hide it.

On train development chromosomes, report these four counters alongside
CDS-base support, fully covered exons, donor/acceptor neighborhoods and
whole-chain support, stratified by exon length, motif class and clade.
Whole-chain comparative support means all its required coding/boundary
positions are supported, not that the entire intron interior is selected.
These counters are a required density/recall audit, not results obtained
in this tick. The exact seed scan, threshold selection and density policy
remain outstanding. In particular, this draft does not assume a scan can
retain 5% allocated support on every genome with useful recall.

The one-pass arithmetic assumes A's current core scores and stem features
stay buffered while B refines its eight eligible cores in that orientation,
then the combined
11-channel emissions are spooled. A decision based on whole-genome A
scores cannot be known at that time. Any global density-cap policy that
first completes A must separately charge storage/reload of its stem and
scores, or their regeneration. This dependency must be resolved before
claiming a concrete CPU-capped inference schedule.

## 5. Candidate C: B's evidence with a splice graph

Use the identical A/B encoders, inputs, masks and training split. Replace
the chain decoder with a directed acyclic graph on genomic start, stop,
donor and acceptor candidates. Exon/intron edges carry strand, phase,
unfinished codon prefixes, sequence and duration scores; only compatible
paths form transcripts under the coding rules in section 3.1.
An intron is an edge between endpoints, so it can cross chunks without
attending to every intronic base. Permit multiple paths and independent
overlapping loci; use the benchmark's all-emitted-transcript precision
to prevent extra isoforms becoming free recall.

A proposed **1.20 M ceiling** includes B and boundary/edge scoring; this
remains a reservation, unlike the explicit A/B layer counts, until C's
edge scorer and candidate-density audit are specified.
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
by the Python standard-library arithmetic below; the rows are stored to
six decimal places. It now counts all specified A/B matrix projections,
convolutions and attention products, including stem, embeddings, fusion
and output heads. Bias additions, layer normalization, nonlinearities,
softmax, pooling, gathers, residual additions, masks, feature extraction,
I/O, tree setup, decoder and synchronization remain **uncounted work**.
Parameter counts include these layers' trainable scalars; FLOP counts do
not price every scalar operation they perform. Small depthwise kernels
also need not sustain the same throughput as large matrix products.

Division by the [cost baseline][cost]'s assumed 3e10 CPU FLOP/s and
1e13 GPU FLOP/s gives **conditional arithmetic times**, not measured
runtimes or guarantees. Section 3.3 separately counts decoder transitions
and scratch traffic; those cannot be timed using a matrix-throughput rate.
The table assumes one neural evaluation followed by emission replay from
scratch, under the streaming dependency in section 4.1. The TSV also
shows the two-neural-evaluation regime. No model, decoder or GPU inference
was implemented to obtain these numbers.

Let N = 1,000,000 bases and K include the target. Here f is the mean
allocated central support across orientations after 384-base tile
expansion, before halos, as defined
in section 4.1. Both encoder counts are normalized long-sequence rates;
for a real assembly, sum `ceil(sequence_length/3072)` A chunks and count
selected B cores and their actual allocated rows. The A output heads can
skip terminal padding. The scenario f is not raw seed coverage, alignment
coverage, nor the genome's annotated CDS fraction. Use the explicit
1.53125 B halo allocation, not the earlier 1.2 placeholder.

```python
N = 1_000_000
nStem = 2 * N * 4104 / 3072          # both strands and stem halos
nA = 2 * N * 342 / 3072             # 342 context tokens per A chunk
nOut = 2 * N                       # heads only on retained positions
flopA = (
    2 * (8*16*9 + 3*(16*9 + 16*16)) * nStem
    + 2 * (16*96) * nA
    + 2 * (4*12*96**2) * nA
    + 4 * 4 * nA * 16 * 96
    + 2 * (112*32 + 32*11) * nOut
)
rows = []
for K, f in [(0, 0), (8, .05), (8, .75), (8, 1),
             (16, .05), (16, .75), (16, 1)]:
    nB = 6 * N * f * K * 196 / 384  # both strands, all frames, halo tokens
    flopB = (
        2 * (27*32 + 2*16*32**2) * nB
        + 2 * 4 * nB * 32 * (32 + K)
        + 2 * (112*32 + 32*11) * nOut * f
    ) if K else 0
    total = flopA + flopB
    rows.append([
        'B' if K else 'A', K, float(f),
        flopA/1e9, flopB/1e9, total/1e9,
        total/3e10, total/1e13,
        2*total/1e9, 2*total/3e10, 2*total/1e13,
    ])
```

| Candidate/scenario | Counted GFLOP/Mb | CPU s/Mb at assumed rate | GPU s/Mb at assumed rate |
|---|---:|---:|---:|
| A, full listed layers | 231.460 | 7.72 | 0.023 |
| B, K=8, f=0.05 | 327.190 | 10.91 | 0.033 |
| B, K=8, f=0.75 | 1667.405 | 55.58 | 0.167 |
| B, K=8, f=1 | 2146.052 | 71.54 | 0.215 |
| B, K=16, f=0.05 | 427.150 | 14.24 | 0.043 |
| B, K=16, f=0.75 | 3166.805 | 105.56 | 0.317 |
| B, K=16, f=1 | 4145.252 | 138.18 | 0.415 |

For A the components are 12.5685 GFLOP/Mb for the stem convolutions,
0.684 for context projection, 202.464 for context blocks/attention, and
15.744 for fusion/output heads: **231.4605 GFLOP/Mb** in total by the
formulas above. The revised B count includes the triplet projection and
all three frame-hypothesis gathers' shared output head, as well as the
larger token allocation. These values supersede the previous core-only
202.069 and sparse-B 274.814 GFLOP/Mb scenarios.

A's engineering targets remain <=15 CPU-s/Mb and <=0.5 GPU-s/Mb.
B at K=8, f=0.05 leaves only **4.09 conditional CPU-s/Mb** of the 15-second
budget for every omitted operation, including decoding and scratch I/O.
That subtraction is arithmetic, not a measured allowance. At K=16 the
same support leaves only **0.76 seconds**. Before pricing any omitted
work, the formal maximum allocated fractions at the assumed CPU rate are
`(450e9-flopA)/flopB(f=1)` = **0.1141 at K=8** and **0.05584 at K=16**.
These are necessary conditions under the throughput assumption, not
support policies or promises of sufficient recall. Dense support exceeds
the CPU ceiling in the counted work alone. A's regenerated-emission
scenario also exceeds it (15.43 CPU-s/Mb before omitted operations).

On a GPU the matrix arithmetic leaves headroom, but it does not predict
small-kernel, I/O or serial decoder efficiency. C inherits B's encoder
work **plus an unpriced graph/traceback cost**; assigning a precise total
seconds/Mb forecast now would invent evidence. A two-pass strand schedule
shares weights, not execution time. All three must include preprocessing
other than the separately displayed alignment-construction bill.

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

- Review the now-specified prefix-state decoder, duration recurrence,
  partial-end rules and scratch/checkpoint accounting. Quantify which
  overlapping transcripts the single-path training control cannot emit;
  retain the Phase 4 correctness checks in section 3.4 as prerequisites
  to any encoder comparison, without implementing the prototype now.
- A/B layer and token inventories are now explicit (sections 3.5, 4.1
  and 6). Review their fixed-grid/halo assumptions and streaming buffers;
  price non-matrix operations, decoder and I/O by measurement only after
  Phase 4 authorization. Specify C's edge scorer and an observable
  candidate/edge-density audit; its 1.20 M allocation remains a ceiling.
- Decide the exact non-neural support scan, density cap/fallback policy and
  support-recall criteria on train development data. Use the explicit tile
  counters, replacing the old 1.2 multiplier, and resolve any whole-genome
  decision's need to buffer or regenerate A's scores and stem features.
- Turn this working artifact into `docs/design/proposal.md` on
  `work/T-human-011-stalin`, open a PR, request reviews from at least two
  other agents through relay, and send the final ranked proposal to human
  for the required Phase 4 decision. No model training was run this tick.

## Sources

Accepted repository documents above are pinned to this draft's main commit.
External sources were read through public pages on 2026-09-15 UTC; the
NCBI genetic-code source was rechecked for the decoder revision. Only
metadata and our own design notes are stored here. Decoder state counts,
storage quantities, layer/tile inventories and example strings are our
proposed specification and arithmetic, not biological measurements or
implemented-model results.

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
