# Candidate A pilot: implementation and end-to-end measurement (T-human-014)

Status: **in progress** (owner `lenin`). This document is the T-human-014
deliverable and is revised in place each tick. Measurements against the
budget are pending a gagarin fitting/inference run; this first revision
records the implementation and the measurement plan so reviewers can check
the design before any GPU time is spent.

## 1. What is implemented

The neural half of candidate A (proposal section 3), on
`work/T-human-014-lenin`:

- `model/a/features.py` — the eight-channel per-base featurizer: ACGT
  indicators, ambiguity, soft-mask, local GC over a centred 129-base window
  (ambiguous bases excluded from the denominator; fixed 0.0 when the window
  has no unambiguous base), and real-sequence availability. No reference
  annotation, no trainable parameters.
- `model/a/encoder.py` — the section 3.5 encoder: a width-16 nucleotide stem
  (kernel-9 conv + three dilated depthwise/pointwise residual blocks,
  dilations 1/2/4, GELU in the residual branch, channel layer norm), a
  context path mean-pooled by 12 and projected to width 96, four pre-norm
  local-attention/MLP blocks (four heads, relative offsets −8..+7, a learned
  per-head offset bias, expansion-4 MLP; attention forms scores and value
  products only over the 16 permitted offsets per query, O(T·W) not O(T²)),
  fine/context fusion of the 16-channel
  stem with the repeated 96-channel context, and the 11-channel emission head.
  Plus the 54 pooled-decoder scalars (`DecoderParams`) that feed the existing
  `model.grammar` reference and delayed decoders.
- `model/a/inventory.py` — the section 3.5 parameter arithmetic, torch-free.

**Parameter count is exactly 455,841**, matching the section 3.5 inventory
row for row (`tests/test_a_encoder.py`, `Inventory` cases). The eleven
emission channels are exactly the channels `model.grammar.Scores` consumes:
U, CDS by prefix length (3), intron by prefix length (3), start, stop, donor,
acceptor. Dependency radius is 491 bases, below the proposed 516-base halo.

| Component | Formula | Count |
|---|---|---:|
| Stem convolution | `8*16*9 + 16` | 1,168 |
| Three residual blocks | `3*(16*9 + 16 + 16*16 + 16 + 2*16)` | 1,392 |
| Pooled-context projection | `16*96 + 96` | 1,632 |
| Four attention/MLP blocks | `4*(12*96**2 + 9*96 + 4*96 + 4*16)` | 447,616 |
| Fine/context fusion | `112*32 + 32` | 3,616 |
| Emission projection | `32*11 + 11` | 363 |
| Pooled decoder | `3*3*2 + 2*16 + 4` | 54 |
| **Total** | | **455,841** |

## 2. What is not yet implemented (next ticks)

1. **Structured loader / admission adapter** (proposal section 3.6): rejoin
   the checksummed raw GFF3/FASTA (manifests omit exon geometry and
   auxiliary-site coordinates), build topology and metadata masks, apply the
   partial/exception masking policy, and emit the required per-species
   audit counts. Reuse the source-join fixtures verified in the prior owner's
   notes (engels-0050 … engels-0057) and the boundary contracts they check.
2. **Differentiable chain loss**: the CRF numerator (constrained to admitted
   chains, summing latent duration components) and denominator over the same
   grammar as `model.grammar`, checked against the expanded reference decoder
   as the specification oracle (section 3.2).
3. **Training and inference entry points** with a config and a run manifest
   (data, seed, commit, hardware), then a gagarin compute-request `alert`.

## 3. Measurement plan (pending gagarin)

Per the task's definition of done and proposal section 6:

- Fit on train-species development chromosomes only; select checkpoints on
  train development chromosomes; run `benchmark/leakage_check.py` before any
  held-out evaluation. S. pombe stays held out and is used only for the
  section 6.1 runtime-normalization evaluation after freezing.
- Measure preprocessing, encoder, decoder, traceback, scratch I/O and output
  separately with `/usr/bin/time -v`, peak host and device memory, and the
  hardware, following the `docs/cost-baseline/measured.tsv` conventions. CPU
  regime on one core; GPU regime with the section 6.1 multi-worker decoder
  accounting (worker count and aggregate memory included).
- Compare against the accepted budgets (15 CPU-s/Mb, 0.5 GPU-s/Mb, 8 GB) and
  against AUGUSTUS on the same machine using the S. pombe normalization; score
  with `benchmark/score.py` on the development chromosomes.
- Add rows to `docs/cost-baseline/measured.tsv` for A on at least S. pombe and
  one metazoan development chromosome, both regimes, and fill the section 6.1
  sensitivity table here with measured numbers.

## 4. Budget and caps

Model ≤ 5 M parameters (455,841 ✓) and ≤ 8 GB. Fitting uses train species
only. T-human-014 compute cap: ≤ 24 GPU-hours; one request may exceed 24 h
wall clock only with a coordinator `decision`. Sampled bases and repeats per
fitting run are declared in the run manifest before the run; actual GPU/CPU
hours are recorded in this task's log after every run.

## 5. Section 6.1 sensitivity table (measured)

_Pending the gagarin run; to be filled with measured stem-efficiency,
decoding and scratch numbers, not the proposal's arithmetic estimates._

## 6. Review responses (PR #38)

First-round encoder review findings and their resolution on
`work/T-human-014-lenin`:

- **Local attention did dense quadratic work** (engels-0059 P2). `LocalAttention.forward`
  now gathers the 16 permitted offsets per query (padded unfold) and forms
  scores/value products over that window only — O(T·W·hd), not O(T²·hd).
  `_dense_forward` retains the masked dense implementation; a torch-gated test
  asserts the two agree in output and input gradient. The fusion/emission heads
  still run over all positions; cropping the halo before those heads is the
  core-interface change deferred to the section 3.6 loader tick.
- **Padding contaminated neighbouring GC** (engels-0059 P2). `gc_track` now
  takes the availability mask and excludes unavailable positions from both GC
  counts, so an invented padding letter cannot shift a real base's window. New
  tests: `gc_track` invariance to the padded letter, and channel-6 invariance
  across `AC`/`AT`/`AN` under `available=[True, False]`.
- **`model.a` was absent from the package list** (engels-0059 P2). Added to
  `[tool.setuptools] packages` so an installed wheel imports it.
- **Duration components initialized identically** (stalin-0062 P2). `DecoderParams`
  now seeds distinct per-component hazard logits (`HAZARD_LOGIT_INIT =
  (-1, 0, 1)`, broadcast across phases; recorded here for the run manifest),
  breaking the symmetry that gave the mixture/hazard logits zero gradient.
  Uniform mixture weights are kept. A test asserts the three components start
  distinct.
