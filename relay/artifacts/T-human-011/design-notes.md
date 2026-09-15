# Geometric gene prediction: candidate design working draft

Task: [T-human-011](../../tasks/T-human-011.md). Author: stalin.
Written 2026-09-15 UTC against accepted main `09b5386` (claim `9239764`);
decoder/resource revision against main `86abfc5`, then layer/tile accounting
against main `f424eb7` later that day. Support/strand and annotation-topology
audit against main `0899198`, then representative-label and carried-quota
audits against main `9fbe9cd` on the same date.
Status: fifth bounded design pass, **in progress**. This is the working
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
meets the CPU budget on gene-dense genomes. Declare A on CPU, B with capped
support on CPU, and B with unrestricted support on GPU as separate regimes.
The cap withholds comparative refinement from excess tiles; A still predicts
those positions. Neither partial comparative support nor GPU matrix arithmetic
establishes the end-to-end time target. Sections 4.2 and 6.1 make those
distinctions explicit.

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
`longest-cds` policy with conflicting same-strand overlaps masked as
specified and audited in section 4.4; keep all real boundary sites for
separate auxiliary losses and count the omitted chains.
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
Sections 4.1 and 4.2 fix the accounting geometry and proposed seed/cap rule;
threshold fitting and its density/recall audit require a trained A in Phase 4.

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
design pass and again in the fourth pass's public-page check (2026-09-15 UTC).
Lenin independently reported the same result from another runner in
[lenin-0040](../../messages/20260915T030649Z-lenin-0040.md). The actual
Tree-RoPE implementation remains unaudited; the eventual decision must
record this access gap. Retain it
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
The sequence-score counters remain a required Phase 4 density/recall audit.
Section 4.3 now measures annotation-only tile requirements on four train
species; those oracle counts are not a sequence candidate detector. This
draft does not assume useful support recall at 5% allocation on every genome.

The one-pass arithmetic assumes A's scores and stem features stay buffered
until B refines that region, then the combined 11-channel emissions are
spooled. Section 4.2 uses bounded groups to satisfy that dependency. A
decision based on whole-genome A scores would instead require separate
storage/reload of its stem and scores, or their regeneration; it is not
part of the first schedule.

### 4.2 A concrete streaming support rule and cap

This is a proposed Phase 4 schedule, not a trained gate. Begin with A alone.
At each oriented base define the local ranking score
`s(t) = max_p emission_CDS[p,t] - emission_U[t]`. It is a score contrast,
not a CRF marginal or a calibrated coding probability. Nominate bases with
`s(t) >= tau`, dilate their intervals by 12 real bases on each side, then
expand to the fixed 384-base cores. The dilation is a design setting.
A tile's rank is the maximum seed score among seeds whose dilated interval
intersects it; ties use increasing oriented coordinate. No minimum ORF or
exon length, canonical motif, alignment coverage, tree, or KA/KS threshold
is a hard eligibility condition. A supplies the only learned scan, already
charged genome-wide; the subsequent max/threshold/dilation pass is
non-neural work whose time must be counted. Separate motif/ORF seed unions
would be additional, independently audited policies.

Fit the single shared threshold on reserved chromosomes of train species
only, with A frozen and before comparing tree arms. Proposed gate targets
before applying a compute cap are at least 99% CDS-base support, 99% support
of each donor/acceptor site's clipped 12-base flanks on each side, and 95%
whole-CDS-chain support, macro-averaged over train development species.
These percentages are design acceptance criteria, not observed recall.
Display each species and the short-exon/noncanonical strata as well.
Among thresholds satisfying all three, choose the one with the lowest
allocated token count; if none does, the gate has failed its acceptance
criterion. Do not relax the targets using held-out annotations. A remains
the usable control, and a predeclared full-tile GPU arm can isolate whether
the comparative encoder works when gate recall is removed as a variable.

**Bounded buffering.** Group 32 consecutive A cores, hence 98,304 oriented
bases and 256 B cores. Retain their fine stem and 11 score channels, obtain
12 bases of right-hand score lookahead from the next A core, and carry
12 scores of left-hand history. The next A core is computed once and kept
for the next group. One group plus that entire lookahead A core reserves
`(98_304 + 3_072) * (16 + 11) * 4 = 10_948_608` bytes for float32
stem/scores, by this inventory. Rank/refine the current group's tiles,
spool its combined emissions in coordinate order, and release its buffer.
Process the opposite orientation separately. This does not require
whole-chromosome A-score storage or neural regeneration. B's alignment
halos can be fetched independently; charge their input buffers, packed
tensors, activations, ranking arrays and I/O separately. The buffer
arithmetic is not a total peak-memory measurement.

**CPU cap.** Use a predeclared `alpha <= 0.05` and K <= 8. A group with n
real bases in that orientation receives
`q = floor(alpha*n/384)` new tile slots. Start a scalar credit balance at
zero for each (sequence, orientation). Add q to the balance, select
`min(balance, eligible_tile_count)` highest-ranked tiles in the current
group, then subtract that selected count. Carry unused slots forward on that same oriented
sequence; discard the final balance. Do not transfer credits across
sequences or orientations, borrow from future groups, revisit released
groups, or refill credits after alignment filtering. The same tile
choices are held fixed across tree arms.

At the illustrative alpha = 0.05 a full group earns 12 slots for its 256
tiles, a 0.046875 incremental allocated fraction. With saved credit it
may refine up to all 256 current tiles, processed in the same bounded
B batches. No additional A group is retained. After any group prefix,
the number selected is at most the sum of q earned in that prefix.
Summing q across groups and strands gives
`384*sum(q)/(2*N) <= alpha`, including short terminal groups. Thus this
bounds allocated B work under the stated K, halo and batch assumptions.
Groups shorter than 7,680 real bases earn no new slots at that alpha,
but may spend prior credit; an entire sequence shorter than that gets
none. Report this short-sequence fallback. Padding across batches must
stay within K=8,
and unselected tiles must never acquire residual outputs from neighbors'
halos. The earlier per-group rule that forfeited all unused q is retained
as a policy ablation; section 4.4 compares it with this carried-credit
schedule using exactly the same earned slots. Credit can only help later
groups: this is not global score ranking, and early false positives can
consume slots needed by later true candidates.

This bounds work, not elapsed time. Freeze alpha from a Phase 4 pilot
that prices the full A/decoder/I/O path first; if A exhausts the CPU
budget, a positive B quota cannot make that run compliant. Alpha = 0.05
is an illustrative upper allocation, not an approved CPU runtime setting.
Keep counters for eligibility, quota truncation, unavailable alignment,
K=1 and final supported outputs, with raw and post-cap recall separately.
Also record earned/spent/forfeited credits and peak group allocation.
Do not select or retry a cap using held-out accuracy. A timer-based cutoff
would change predictions with machine/load and is not the first policy.

**Separate benchmark rows.**

| Declared regime | Comparative allocation | Interpretation |
|---|---|---|
| A-CPU | None | Common sequence/decoder control; measure its full CPU time. |
| B-capped-CPU | Above ranked quota, frozen alpha, K<=8 | Partial comparative rescue; A supplies every omitted position. |
| B-gated-GPU | All nominated tiles, fixed tau and K | Gate recall still limits comparative improvement. |
| B-full-GPU | Every tile, fixed K | Removes seed/quota omissions; still limited by actual alignment. |

GPU labels identify proposed execution regimes, not evidence of meeting
0.5 s/Mb. In particular, a CPU-hosted chromosome decoder can dominate GPU
end-to-end wall time. A dense genome need not become literally A-only
under a positive cap: some tiles still receive B. Its capped result must
not be reported as the unrestricted comparative result. Run the matched A
control on GPU too when measuring comparative overhead on that device.

### 4.3 Annotation-only support and single-path topology audit

Measured locally in the fourth pass from four cached **train** species' panel
GFFs, verifying every `md5_gff` before parsing. This audit uses the accepted
scorer's sequence and transcript filters, including the 10 kb sequence
floor and declared pseudogene/biotype exclusions. Its denominators therefore
need not equal all-primary-assembly counts in `cds_union.tsv`.
All retained CDS-bearing chains are included, including partial ones that
the scorer retains for structure metrics. No held-out annotation was
opened for this audit; no threshold, weights or development split was fitted.

For each retained CDS interval, select every 384-base core it intersects
on its **own** orientation, with the reverse grid anchored at the reverse
sequence origin. Count every tile once per strand. This is the minimum
tile set covering these annotated CDS bases on this grid, before boundary
flanks, alignment filtering, and false-positive sequence candidates.
It is an oracle diagnostic, not a runnable prediction input or the cost
of every possible sparse comparative model. Compute K=8 estimates from
section 6 using its measured f; do not interpret them as timings.

| Train species | Scored sequences / bp | Plus / minus CDS tiles | Oracle allocated f | Counted conditional B CPU-s/Mb |
|---|---:|---:|---:|---:|
| S. cerevisiae | 16 / 12,071,326 | 14,355 / 14,030 | 0.451476 | 36.528 |
| C. elegans | 6 / 100,272,607 | 64,263 / 62,507 | 0.242737 | 23.207 |
| D. melanogaster | 119 / 140,264,371 | 45,432 / 45,944 | 0.125079 | 15.698 |
| M. musculus | 58 / 2,728,195,764 | 142,625 / 142,088 | 0.020037 | 8.994 |

The first three exceed 15 conditional CPU-s/Mb from counted encoder work
alone if every retained CDS must have comparative support. Mouse's oracle
allocation leaves room only under that arithmetic; boundary flanks,
false-positive candidates, poor alignment and omitted operations may remove
it. These results are neither a guarantee of near-complete recall at a
smaller cap nor a requirement that B refine every CDS to help accuracy.

**Strand denominator correction.** The [cost baseline][cost] reports a
physical-coordinate CDS union, whereas f divides orientation-specific
allocated bases by **2N**. If U is the physical union and U+ and U- the two
strand unions, then `U <= U+ + U- <= 2U`; before tile expansion the
corresponding fraction lies between U/(2N) and U/N. Do not directly
substitute its 72.35% yeast CDS union for f. In this filtered yeast audit
the physical union is 8,726,635 bp and the summed oriented union is
8,735,738 bp; tiling yields f=0.451476. The cost note's qualitative warning
about dense genomes survives, but the conversion requires this strand
accounting rather than an equality between physical CDS coverage and f.
No unmeasured S. pombe coding fraction is needed.

**Single-path topology.** Collapse identical (sequence, strand, CDS-block)
chains, map each remaining chain to its full first-to-last-CDS span, and
greedily pack intervals in increasing end order separately per sequence
and strand. Earliest-finish interval scheduling gives the maximum number
of pairwise non-overlapping spans: choosing an earlier finishing eligible
span cannot remove a later choice available after a later finishing one.
The CRF keeps its chain state through introns, so even a gene nested wholly
inside a same-strand intron conflicts. This calculation deliberately
ignores codon legality, partial-end restrictions and minimum gaps; it is
an optimistic bound for simultaneously retaining all these reference
chains with one path per strand.

| Train species | Retained transcript IDs | Distinct CDS chains | Maximum non-overlapping distinct chains |
|---|---:|---:|---:|
| S. cerevisiae | 6,002 | 6,002 | 5,933 |
| C. elegans | 28,590 | 28,590 | 20,081 |
| D. melanogaster | 30,746 | 22,438 | 13,849 |
| M. musculus | 97,324 | 68,828 | 22,436 |

This is **not a ceiling on benchmark transcript F1**. The accepted scorer's
`_transcripts` matches emitted isoforms within matched loci and does not
count unused reference isoforms there as false negatives. Nor is the last
column a count of distinct loci. Section 4.4 applies the actual longest-CDS
choice and a proposed conflict mask, then quantifies affected loci and
boundary supervision. Opposite-strand overlaps are already
allowed by separate paths.

**Reproduction.** Run the following read-only arithmetic/annotation audit
from this checkout after placing the panel's named GFFs under
`/tmp/bench007/data/<species>/` (the path used here). It emits the counts
above plus hashes and union lengths; it stores no full text, sequence,
model, or new data files. GFF checksums come from the panel. At this audit,
`benchmark/score.py` SHA-256 was
`c738514d1a3d6b73535c543acf8baa9e12d940ed607863182e027baa4096ea24`,
and `benchmark/panel.tsv` SHA-256 was
`c9fc2b6925ac10e94164a210e9edd62a8e0f45a4e79608e020e0531a2a9b9342`;
the code prints both to expose a changed filter or panel.

```python
from collections import defaultdict
from pathlib import Path
import csv, hashlib, json, runpy

score = runpy.run_path("benchmark/score.py")
panel = {r["species"]: r for r in csv.DictReader(open("benchmark/panel.tsv"), delimiter="\t")}
species = ["Saccharomyces_cerevisiae", "Caenorhabditis_elegans",
           "Drosophila_melanogaster", "Mus_musculus"]
print("score_sha256", hashlib.sha256(Path("benchmark/score.py").read_bytes()).hexdigest())
print("panel_sha256", hashlib.sha256(Path("benchmark/panel.tsv").read_bytes()).hexdigest())
for sp in species:
    row = panel[sp]
    assert row["split"] == "train"
    path = Path("/tmp/bench007/data") / sp / (row["ftp_dir"] + "_genomic.gff.gz")
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    assert md5 == row["md5_gff"], (sp, md5)
    ann = score["load_gff"](str(path))
    seqids = score["select_seqids"](ann)
    chains = score["select_transcripts"](
        {t: c for t, c in ann.chains().items() if c[0] in seqids}, ann)
    unique = set(chains.values())
    N = sum(ann.seq_len[s] for s in seqids)
    oriented, physical, spans = defaultdict(list), defaultdict(list), defaultdict(list)
    tiles = {"+": set(), "-": set()}
    for sid, strand, blocks in unique:
        assert strand in tiles
        L = ann.seq_len[sid]
        spans[sid, strand].append((min(a for a, b in blocks), max(b for a, b in blocks)))
        for lo, hi in blocks:
            assert 1 <= lo <= hi <= L
            oriented[sid, strand].append((lo, hi))
            physical[sid].append((lo, hi))
            a, b = (lo-1, hi) if strand == "+" else (L-hi, L-lo+1)
            tiles[strand].update((sid, t) for t in range(a//384, (b-1)//384+1))
    def union_bp(groups):
        return sum(e-s+1 for iv in groups.values()
                   for s, e in score["merge_intervals"](iv))
    max_pack = 0
    for intervals in spans.values():
        last = 0
        for lo, hi in sorted(intervals, key=lambda v: (v[1], v[0])):
            if lo > last:
                max_pack += 1
                last = hi
    f = 384 * sum(map(len, tiles.values())) / (2*N)
    result = dict(species=sp, md5_gff=md5, scored_sequences=len(seqids), scored_bp=N,
                  transcripts=len(chains), unique_CDS_chains=len(unique),
                  physical_CDS_union_bp=union_bp(physical),
                  oriented_CDS_union_bp=union_bp(oriented),
                  plus_tiles=len(tiles["+"]), minus_tiles=len(tiles["-"]),
                  oracle_CDS_allocated_f=f, B8_conditional_CPU_s_per_Mb=7.71535+1914.592*f/30,
                  max_nonoverlapping_unique_chains=max_pack)
    print(json.dumps(result, sort_keys=True), flush=True)
```

### 4.4 Representative-label losses and quota sensitivity

This fifth-pass audit uses the same four **train** GFFs and scorer/panel
hashes as section 4.3. It supplements the scorer's CDS records with the
GFF's full exon lengths and first transcript appearance in file order,
then calls the existing cutter's actual `select_isoforms(..., "longest-cds")`
function. The adapter preserves the scorer's namespaced gene IDs and
sequence/strand identity. Every retained transcript has exon records in
these files, so no fabricated CDS-only tie breaker is needed. This is an
annotation audit using a proposed GFF adapter, not evidence that the
current UCSC/Ensembl window loader already implements structured training.

**Conflict-mask contract.** After one representative per locus is chosen,
form connected components of overlapping *full CDS spans*, separately
on each sequence and strand. Spans are half-open; touching endpoints do
not overlap. Mask every representative in a component containing more
than one locus and make that component's complete span unconstrained in
the chain-loss numerator. This includes same-strand genes nested inside
another selected gene's intron. Do not mask only the CDS intersections,
assign the masked region an intergenic target, or cut and rejoin fragments
from different transcripts. Keep the original ordered chains for the
other representatives. This conservative policy discards every member
of a conflict rather than choosing an arbitrary winner; it is a training
policy, not the interval-packing optimum of section 4.3.

| Train species | Retained transcript IDs | Selected loci / representatives | Nonrepresentative IDs omitted from chain target | Conflicting components / representatives masked | Representatives after topology mask |
|---|---:|---:|---:|---:|---:|
| S. cerevisiae | 6,002 | 6,002 | 0 | 69 / 138 | 5,864 |
| C. elegans | 28,590 | 19,971 | 8,619 | 65 / 132 | 19,839 |
| D. melanogaster | 30,746 | 13,965 | 16,781 | 234 / 530 | 13,435 |
| M. musculus | 97,324 | 22,183 | 75,141 | 137 / 370 | 21,813 |

The masked representatives are respectively 2.2992%, 0.6610%, 3.7952% and
1.6679% of the selected loci, derived from the preceding counts. The union
of their component spans is 263,347 / 604,025 / 6,070,755 / 25,019,190
oriented bases: 1.0908% / 0.3012% / 2.1640% / 0.4585% of the corresponding
**2N** scored-base denominators in section 4.3. These are full-span masks,
including introns; they are neither physical CDS coverage nor an accuracy
ceiling. Nonrepresentative transcript IDs include identical CDS chains
with different UTRs, so their count is not the number of lost distinct
coding structures.

**Boundary supervision stays separate.** Count distinct
(sequence, strand, kind, position) targets over all accepted transcripts,
over the selected representatives, and over representatives surviving
the topology mask. Here start and stop mean the first and last CDS base
under the cutter's boundary convention, omitting ends declared partial by
the scorer. Donor/acceptor are the first/last intron base between CDS
blocks with a gap of at least 20 bases; UTR junctions and shorter gaps are
not included. Missing declarations are not a sequence-confirmed intact
codon: this audit opens no FASTA.

| Train species / chain-label stage | Start sites | Stop sites | CDS donors | CDS acceptors |
|---|---:|---:|---:|---:|
| S. cerevisiae: all / selected | 5,960 | 6,001 | 281 | 281 |
| S. cerevisiae: after topology mask | 5,864 | 5,864 | 274 | 274 |
| C. elegans: all | 24,350 | 21,195 | 103,593 | 103,556 |
| C. elegans: selected | 19,951 | 19,965 | 101,609 | 101,611 |
| C. elegans: after topology mask | 19,824 | 19,835 | 100,950 | 100,950 |
| D. melanogaster: all | 16,316 | 16,111 | 44,113 | 44,500 |
| D. melanogaster: selected | 13,943 | 13,954 | 41,268 | 41,273 |
| D. melanogaster: after topology mask | 13,432 | 13,432 | 39,354 | 39,354 |
| M. musculus: all | 35,881 | 30,159 | 187,365 | 189,666 |
| M. musculus: selected | 22,108 | 22,084 | 178,957 | 178,919 |
| M. musculus: after topology mask | 21,771 | 21,785 | 177,047 | 177,047 |

Train the auxiliary boundary losses against the **all** sets, on each
orientation and with independent boundary channels; the topology and
isoform policies must not turn known alternative sites into auxiliary
negatives. These are soft auxiliary targets, not simultaneous hard
constraints forcing every site into the CRF path. In particular, preserve
sites inside chain-masked components. The reproduction prints those
counts too. Translate cutter positions into the decoder's oriented
coordinates: its acceptor score is attached to the first following CDS
base, one oriented base after the cutter's last-intron-base mark.
Duplicated start/stop positions shared by different genes count once;
this explains why even yeast's all-site counts need not equal its loci.

The surviving representatives are **not yet certified legal CRF targets**.
There are 0 / 15 / 3 / 45 with a declared partial CDS end after the
topology mask (0 / 17 / 3 / 51 before it). A partial end away from an
actual sequence edge cannot simply receive section 3.1's edge-entry rule.
Phase/translation consistency, recoding, short gaps and sequence ambiguity
still need the declared grammar audit before Phase 4 fitting. Preserve
these records as exceptions or unconstrained spans; never force an
incompatible path into a numerator with zero compatible probability.
The counts above isolate topology/isoform losses without claiming that
these remaining checks have passed.

**Quota comparison.** [Marx's independent audit][quota-note] reproduced
section 4.3 and identified wasted per-group slots. Recomputing from the
GFFs with perfect annotation-only tile ranking gives the following counts.
All three policies receive the *identical* sum of per-group
`floor(0.05*n/384)` slots; changing rounding is not part of this comparison.
Pooling is an offline upper bound per (sequence, orientation), independent
of candidate scores. The carried policy follows section 4.2 with zero
initial credit; an eligible tile can be served only while its group is
buffered.

| Train species | Oracle CDS tiles | Served by old per-group rule | Served by carried credit | Served by offline pooling |
|---|---:|---:|---:|---:|
| S. cerevisiae | 28,385 | 2,950 (10.3928%) | 2,950 (10.3928%) | 2,950 (10.3928%) |
| C. elegans | 126,770 | 24,336 (19.1970%) | 24,476 (19.3074%) | 24,480 (19.3106%) |
| D. melanogaster | 91,376 | 26,776 (29.3031%) | 31,333 (34.2902%) | 32,929 (36.0368%) |
| M. musculus | 284,713 | 168,501 (59.1828%) | 282,530 (99.2333%) | 284,700 (99.9954%) |

For mouse, the two online policies have exactly 666,054 earned slots.
Carrying credit changes where those slots can be spent; it neither
increases alpha nor uses a whole-sequence stem/score buffer. Pooled mouse
coverage rounds to 100% in the inbox note but leaves 13 tiles unserved
under the exact per-sequence/orientation quota. The carried schedule
leaves 2,183, reflecting its additional prefix-budget constraint. It
cannot promise complete coverage or score-optimal spending in a real
scan with false positives. Its forward direction is a declared property
of the oriented schedule; test coordinate/strand mapping and sequence-edge
sensitivity rather than assuming it ranks an entire chromosome globally.

Dilating each annotated CDS block by 12 real bases on either side
reproduces the inbox note's small oracle-density increase:

| Train species | Flanked oracle tiles | Flanked allocated f | Served by old / carried / pooled quotas |
|---|---:|---:|---:|
| S. cerevisiae | 28,670 | 0.456010 | 2,950 / 2,950 / 2,950 |
| C. elegans | 129,379 | 0.247732 | 24,340 / 24,476 / 24,480 |
| D. melanogaster | 92,821 | 0.127057 | 26,844 / 31,348 / 32,937 |
| M. musculus | 295,051 | 0.020765 | 170,880 / 292,587 / 295,036 |

Those counts increase oracle allocated support by 1.0041% / 2.0581% /
1.5814% / 3.6310% relative, computed from the two tile tables. Dilation of
every CDS block is an oracle flank sensitivity, not the actual learned
score-seed dilation in section 4.2. No alignment filtering, trained seed
scores, held-out data or runtime measurements enter either table.

**Reproduction.** The following read-only script uses the same cached
paths as section 4.3, verifies GFF MD5s before parsing, asserts train-only
membership and prints the counts above. It calls the existing selector
and boundary helper, with the full CDS block list temporarily supplied
as exons only for *CDS-boundary counting*. The longest-CDS selection uses
the real full exon list. At this revision `scripts/data/cut_windows.py`
has SHA-256
`235af56d0a605b2b6bc560aa034c8a7e13fe56c0b8a25c1d77a89e40fab5a6ce`;
the other two hashes remain those in section 4.3. All are printed.

<!-- representative-quota-audit -->
```python
from collections import Counter, defaultdict
from pathlib import Path
import bisect, csv, gzip, hashlib, json, runpy, sys
sys.dont_write_bytecode = True
score = runpy.run_path("benchmark/score.py")
cut = runpy.run_path("scripts/data/cut_windows.py")
panel_path = Path("benchmark/panel.tsv")
panel = {r["species"]: r for r in csv.DictReader(panel_path.open(), delimiter="\t")}
for source in ("benchmark/score.py", "scripts/data/cut_windows.py", "benchmark/panel.tsv"):
    print(json.dumps({"source": source, "sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest()}), flush=True)

def components(txs):
    # Zero-based, half-open FULL CDS spans, grouped by sequence AND strand.
    by = defaultdict(list)
    for t in txs:
        by[t["seqid"], t["strand"]].append(t)
    masked, masks = set(), defaultdict(list)
    for key, group in by.items():
        members, lo, hi = [], None, None
        for t in sorted(group, key=lambda t: (t["cds"][0][0], t["cds"][-1][1], t["id"])):
            a, b = t["cds"][0][0], t["cds"][-1][1]
            if members and a >= hi:
                if len(members) > 1:
                    masked.update(members); masks[key].append((lo, hi))
                members, lo, hi = [], None, None
            if not members:
                lo, hi = a, b
            members.append(t["id"]); hi = max(hi, b)
        if len(members) > 1:
            masked.update(members); masks[key].append((lo, hi))
    return masked, masks

def sites(txs):
    # CDS-chain boundary convention, including short-gap exclusions.
    out = set()
    for t in txs:
        cds_only = dict(t, exons=t["cds"])
        out.update((t["seqid"], kind, strand, pos)
                   for kind, strand, pos in cut["distinct_sites"]([cds_only])
                   if kind in ("start", "stop", "donor", "acceptor"))
    return out

def count_sites(ss):
    c = Counter(kind for sid, kind, strand, pos in ss)
    return {kind: c[kind] for kind in ("start", "stop", "donor", "acceptor")}

def load_txs(path, ann, chains):
    # Keep scorer gene/sequence/strand namespacing; supplement full exon lengths
    # and first appearance in GFF order for the cutter's real tie breakers.
    extras = {}
    with gzip.open(path, "rt") as fh:
        for line_no, line in enumerate(fh):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            raw = (score["_attr"](f[8], "Parent") or score["_attr"](f[8], "transcript_id")) if f[2] in ("exon", "CDS") else (score["_attr"](f[8], "ID") or score["_attr"](f[8], "transcript_id"))
            if not raw:
                continue
            ids = raw.split(",") if f[2] == "exon" else [raw.split(",")[0]]
            for raw_id in ids:
                tid = "\x00".join((f[0], f[6], raw_id))
                if tid not in chains:
                    continue
                meta = extras.setdefault(tid, {"order": line_no, "exons": []})
                if f[2] == "exon":
                    meta["exons"].append((int(f[3])-1, int(f[4])))
    assert set(extras) == set(chains), "missing GFF transcript order"
    txs = []
    for tid in sorted(chains, key=lambda tid: extras[tid]["order"]):
        sid, strand, blocks = chains[tid]
        cds = [(a-1, b) for a, b in blocks]
        exons = sorted(extras[tid]["exons"])
        assert exons, ("missing full exons for longest-cds tie breaker", tid)
        txs.append(dict(id=tid, gene=ann.gene_of[tid], seqid=sid, strand=strand,
                        start=exons[0][0], end=exons[-1][1], cds=cds, exons=exons,
                        _partial5=tid in ann.partial5, _partial3=tid in ann.partial3))
    return txs

def tile_audit(txs, lengths, flank):
    tiles = defaultdict(set)
    for t in txs:
        sid, strand = t["seqid"], t["strand"]
        L = lengths[sid]
        for a, b in t["cds"]:
            a, b = max(0, a-flank), min(L, b+flank)
            if strand == "-":
                a, b = L-b, L-a
            tiles[sid, strand].update(range(a//384, (b-1)//384+1))
    served_local = served_carried = served_pooled = slots = groups_exceeding = 0
    for sid, L in lengths.items():
        for strand in ("+", "-"):
            counts = Counter(t//256 for t in tiles[sid, strand])
            credit = local = carried = qtotal = 0
            for g, lo in enumerate(range(0, L, 98304)):
                n = min(98304, L-lo)
                q = n//7680  # exact floor(0.05*n/384)
                demand = counts[g]
                qtotal += q; local += min(q, demand)
                credit += q
                take = min(credit, demand)
                carried += take; credit -= take
                assert carried <= qtotal
                groups_exceeding += int(demand > q)
            served_local += local; served_carried += carried
            served_pooled += min(qtotal, len(tiles[sid, strand]))
            slots += qtotal
    total = sum(map(len, tiles.values()))
    return dict(flank=flank, tiles=total, oracle_f=384*total/(2*sum(lengths.values())),
                slots=slots, served_local=served_local, served_carried=served_carried,
                served_pooled=served_pooled, groups_over_local_quota=groups_exceeding)

for sp in ("Saccharomyces_cerevisiae", "Caenorhabditis_elegans",
           "Drosophila_melanogaster", "Mus_musculus"):
    row = panel[sp]
    assert row["split"] == "train"
    path = Path("/tmp/bench007/data") / sp / (row["ftp_dir"] + "_genomic.gff.gz")
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    assert md5 == row["md5_gff"], (sp, md5)
    ann = score["load_gff"](str(path))
    seqids = score["select_seqids"](ann)
    chains = score["select_transcripts"]({t:c for t,c in ann.chains().items() if c[0] in seqids}, ann)
    txs = load_txs(path, ann, chains)
    selected, dropped, fallback = cut["select_isoforms"](txs, "longest-cds")
    assert not fallback
    assert len(selected) == len(score["loci_of"](ann, chains, seqids))
    masked, masks = components(selected)
    retained = [t for t in selected if t["id"] not in masked]
    assert not components(retained)[0]
    all_sites, selected_sites, retained_sites = sites(txs), sites(selected), sites(retained)
    assert retained_sites <= selected_sites <= all_sites
    lengths = {sid:ann.seq_len[sid] for sid in seqids}
    N = sum(lengths.values())
    cds_bp_masked = sum(sum(b-a for a,b in t["cds"]) for t in selected if t["id"] in masked)
    selected_cds_bp = sum(sum(b-a for a,b in t["cds"]) for t in selected)
    starts = {key:[a for a,b in iv] for key,iv in masks.items()}
    def in_mask(site):
        sid, kind, strand, pos = site
        key = sid, "+" if strand == 1 else "-"
        i = bisect.bisect_right(starts.get(key, []), pos)-1
        return i >= 0 and pos < masks[key][i][1]
    result = dict(species=sp, md5_gff=md5, transcripts=len(txs), loci=len(selected),
                  nonrepresentative_transcripts=len(dropped),
                  overlapping_components=sum(map(len,masks.values())),
                  masked_representatives=len(masked), retained_representatives=len(retained),
                  mask_oriented_span_bp=sum(b-a for iv in masks.values() for a,b in iv),
                  scored_oriented_bp=2*N, selected_cds_bp_sum=selected_cds_bp,
                  masked_selected_cds_bp_sum=cds_bp_masked,
                  all_sites=count_sites(all_sites), selected_sites=count_sites(selected_sites),
                  retained_sites=count_sites(retained_sites),
                  auxiliary_sites_in_mask=count_sites({s for s in all_sites if in_mask(s)}),
                  selected_partial_transcripts=sum(t["_partial5"] or t["_partial3"] for t in selected),
                  retained_partial_transcripts=sum(t["_partial5"] or t["_partial3"] for t in retained),
                  quota=[tile_audit(txs,lengths,flank) for flank in (0,12)])
    print(json.dumps(result,sort_keys=True),flush=True)
```

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

### 6.1 Sensitivity to stem efficiency, decoding and scratch

[Marx's cost note](../../messages/20260915T025415Z-marx-0044.md)
identified three costs hidden by a common matrix-throughput divisor.
The following independently recomputed scenarios expose their effect;
the efficiency factors, cycle costs and bandwidth are **assumptions**,
not measured or asserted typical values. Let eta be the stem's throughput
relative to 30 GFLOP/s, c the average CPU cycles per counted decoder
candidate, and v the sustained scratch bandwidth in decimal MB/s:

`T_A = 218.892/30 + 12.5685/(30*eta) + (1.648/2.8)*c + 176/v`

`T_B(K=8,f) = T_A + (1914.592/30)*f`.

The 218.892 GFLOP/Mb term is the remainder of A after its 12.5685 stem;
the candidate and traffic counts come from section 3.3. The 2.8 GHz clock
is the named cost-baseline CPU, not a decoder throughput measurement.
These are additive serial CPU scenarios, with v=500 below. They still
omit nonlinearities, normalization, packing, donor buffers and other
decoder operations, reading alignments, and output overhead; no entry is
a full runtime forecast. More decoder successors at ambiguous bases are
also absent from the canonical-input candidate count.

| Stem eta | Cycles/candidate c | A conditional CPU-s/Mb | B at f=0.05 | Maximum f left under 15 s/Mb, before remaining omissions |
|---|---:|---:|---:|---:|
| 1 | 3 | 9.833064 | 13.024051 | 0.080961 |
| 1 | 10 | 13.953064 | 17.144051 | 0.016405 |
| 0.1 | 3 | 13.603614 | 16.794601 | 0.021880 |
| 0.1 | 10 | 17.723614 | 20.914601 | None; A already exceeds target |
| 1/30 | 3 | 21.982614 | 25.173601 | None; A already exceeds target |

Reproduce each row with the two expressions above and
`f_max = (15-T_A)*30/1914.592` when T_A <= 15. For example, decreasing
the stem rate tenfold while assuming three cycles per candidate removes
the 5% cap's conditional headroom. These scenarios justify measuring A's
full execution first and retaining alpha=0 as an explicit fallback;
they do not identify a feasible positive cap without that measurement.

The portable CPU target in [cost] section 5.2 is 1/11 of AUGUSTUS's
CPU-s/Mb on **S. pombe on the same machine**. It is not an 11-fold speedup
promise for each genome: the baseline's human-chromosome-21 comparison
is about 1.6-fold at the absolute ceiling. Preserve that denominator in
the eventual proposal. On GPU, the dense K=8 counted encoder rate is
0.214605 conditional GPU-s/Mb (section 6), but a serial host decoder
plus synchronization and scratch cannot be inferred to fit the remaining
wall time. GPU-resident decoding or a different decoder is a Phase 4
engineering question, not an accomplished optimization.

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
  partial-end rules and scratch/checkpoint accounting. Sections 4.3 and
  4.4 now quantify topology, longest-CDS selection, full-span conflict
  masks and boundary supervision. Specify the remaining grammar-exception
  handling, especially partial CDS ends inside sequences; retain the
  Phase 4 correctness checks in section 3.4 as prerequisites
  to any encoder comparison, without implementing the prototype now.
- A/B layer and token inventories are now explicit (sections 3.5, 4.1
  and 6). Review their fixed-grid/halo assumptions and streaming buffers;
  price non-matrix operations, decoder and I/O by measurement only after
  Phase 4 authorization. Specify C's edge scorer and an observable
  candidate/edge-density audit; its 1.20 M allocation remains a ceiling.
- Review the now-specified score scan, bounded buffering and carried-credit
  tile quota (sections 4.2 and 4.4). Threshold/alpha fitting and sequence-gate
  recall must wait for the Phase 4 trained-A pilot; the annotation-only audit establishes
  neither. Preserve the distinction between gate targets, work bounds
  and measured runtime, including the sensitivities in section 6.1.
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
proposed specification and arithmetic, not implemented-model results.
Sections 4.3 and 4.4 separately report reproducible measurements of cached
training annotations, under the named scorer/panel hashes; they are not model accuracy
measurements. Section 4.4 also pins the cutter used for representative and
boundary selection. The public axomeme page was rechecked in the fourth pass.

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
[quota-note]: ../../messages/20260915T035634Z-marx-0045.md
