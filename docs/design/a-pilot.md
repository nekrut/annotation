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
- `model/a/loss.py` — the chain-loss **oracle** (sections 3.1–3.2), torch-free:
  `numerator_scores` builds the gold chain's hard `-inf` support mask over the
  eleven emission channels (U only intergenic, a coding channel only on a CDS
  base, an intron channel only on an intron base, `start`/`stop`/`donor`/
  `acceptor` only at their gold coordinates), and `chain_nll` returns
  `log Z − log Z_num` on the existing `model.grammar` decoders — the free
  partition minus the partition restricted to that support, which marginalizes
  exactly the latent frame (grammar-fixed) and the intron duration mixture.
  `tests/test_a_loss.py` pins the mask, `loss ≥ 0`, that the masked numerator
  Viterbi-decodes back to the exact gold chain (single- and two-exon fixtures),
  the `log 2` zero-emission case, and the gradient signs a training loss must
  have (finite-difference against the marginal difference). This is the
  standard-library reference the PyTorch loss matches on the same fixtures.
- `model/a/torch_loss.py` — the **differentiable (PyTorch) chain loss**,
  torch-gated. `chain_nll(x, emissions, cds_ranges, intron_ranges)` returns
  `log Z − log Z_num` on a `(11, n)` emission tensor (`CHANNEL_ORDER`, the same
  eleven channels the encoder head emits), computed with `torch.logsumexp` so
  autograd yields `dL/de = P_free(e) − P_num(e)`, the CRF marginal difference a
  training step applies to the emission head. It is **parity-by-construction**:
  the forward reuses the reviewed grammar state machine unchanged
  (`ReferenceDecoder.initial`/`transitions`/`terminal`) over a `TorchScores`
  view of the emission tensor, so the only substitution against the oracle is
  Python-float `+`/`logsumexp` for torch `+`/`torch.logsumexp`; `support_mask`
  builds the additive `-inf`/`0` mask that reproduces `numerator_scores`.
  `tests/test_a_torch_loss.py` (torch-gated) pins mask/oracle agreement, loss
  value parity on the section-3.4 fixtures (zero and non-zero emissions), the
  `log 2` case, `loss ≥ 0`, and that autograd's gradient at the gold start base
  equals the oracle's central difference. This has the *reference* recurrence's
  cost (explicit mandatory-intron states, a Python dict per boundary), so it is
  the differentiable **reference** the fast vectorized delayed-entry kernel is
  checked against — the same relationship `DelayedEntryDecoder` has to
  `ReferenceDecoder` — not yet the training-window kernel. Complete targets
  only: an edge-enabled decoder is rejected, as in the oracle.
- `model/a/dataset.py` — the section 3.6 structured **training-window loader**,
  torch-free at its core. The species `*.summary.json` **is** the loading
  interface: `iter_windows(summary, gff, fasta)` calls `verify_source` (a hard
  checksum gate recomputing the GFF3/FASTA MD5s and refusing any file that does
  not match the summary's digests, FASTA presence included) before yielding, and
  reads the audit `m`/`table` from that same summary, so it can never audit
  unpinned inputs or use settings that disagree with the committed manifest.
  Admission is delegated to `model.labels.admission.audit_species`, so the
  yielded set is exactly the admitted representatives — the loader never
  re-derives the metadata/sequence audit. `oriented_chain` supplies each chain's
  merged, strand-corrected, flank-padded CDS/intron coordinates; the window
  itself is rebuilt case-preserving from the raw FASTA slice (`_oriented_window`)
  so the featurizer's soft-mask channel survives orientation, with a guard that
  the two agree up to case. Each `WindowExample` carries the oriented window, the
  merged half-open CDS/intron ranges (exactly `numerator_scores`' convention),
  the genetic-code id and the source `(seqid, strand, tid)` identity, with
  `.support()` delegating to `numerator_scores` and a torch-gated `.features()`.
  `iter_windows` is scoped to the current loss's domain: it yields only **clean**
  windows (no transcript of a different gene overlaps the window; the rest are
  counted in `LoaderStats.skipped_neighbor`, since `numerator_scores` would
  otherwise force a neighbour's coding/masked bases to intergenic `U`), skips
  (and counts) edge-partial admitted chains — `chain_nll` scores only complete
  targets — and optionally skips (and counts) windows longer than a caller's
  `max_window`, never turning the encoder core into a gene-length cap. Because the
  audit loads `benchmark/score.py` by path, the loader is a checkout-time tool;
  `model.a.__init__` imports it lazily so `import model.a` needs neither
  `model.labels` nor `benchmark/`. `tests/test_a_dataset.py` builds synthetic
  species and checks the checksum gate (match, tampered FASTA, FASTA-presence
  mismatch, and the gate enforced through `iter_windows`), `m`/`table` binding,
  soft-mask preservation on both strands, the clean-window neighbour skip, the
  window coordinates and identity, the `max_window` skip counter, and the full
  round trip: the loaded example's support-masked numerator Viterbi-decodes back
  to the exact gold CDS/intron chain through `chain_nll`.

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

1. **Structured loader — boundary and crop extensions** (proposal section 3.6):
   the complete-chain loader (`model/a/dataset.py`, above) rejoins the
   checksummed source and yields admitted complete windows. What remains is the
   boundary-support path — edge-partial admitted chains (compatible entry/exit
   families and a first-row phase carried alongside the CDS intervals, scored by
   an edge-enabled decoder) which are currently skipped and counted — and
   cropping whole genes longer than the encoder core into chunks with retained
   intron/codon/duration state (the crop-integration contracts checked in the
   prior owner's notes engels-0045…engels-0057). Batched multi-window collation
   for the torch training step also belongs here.
2. **Differentiable (PyTorch) chain loss** — the reference forward is done
   (`model/a/torch_loss.py`, above): its autograd yields
   `dL/de = posterior_free − posterior_num` and it matches the `model/a/loss.py`
   oracle on the section-3.4 fixtures. What remains is the **fast vectorized
   kernel** — a delayed-entry forward on the encoder emissions with batched
   multi-window collation and chunk-seam handling, checked against this
   differentiable reference — and the **boundary conditioning at training-crop
   edges** (section 3.6, edge-partial numerators). The current forward has the
   reference recurrence's cost and is unsuitable for a full training window at
   the default `m=20`; it establishes the differentiable contract, not the
   training throughput.
3. **Boundary conditioning at training-crop edges** (edge-partial numerators via
   an edge-enabled decoder) and the batched multi-window collation the fast
   kernel needs. The training entry point (item below) already runs, one window
   per accumulated gradient step, over the complete-target windows.

The **training and measurement entry point** is now implemented:
`model/a/train.py` wires the reviewed encoder, loader and reference torch loss
into a runnable program with a declared config and a run manifest. It stays
torch-free on import (torch is imported lazily inside the fitting/measurement
functions), so its config validation, stride padding, by-sequence train/dev
split and manifest construction are unit-tested here without torch
(`tests/test_a_train.py`, 14 stdlib cases). The tensor path — encoder forward,
`chain_nll` autograd, checkpoint selection, and the decode/traceback timing —
runs on gagarin; no local host has torch. See section 3.

## 3. Measurement plan and the entry point

`model/a/train.py` has two subcommands:

- `train --config CONFIG.json` fits candidate A on the *train* sequences of the
  configured train species and selects the checkpoint on the declared
  `dev_seqids` only (the split is by sequence, so a dev chromosome never
  contributes a gradient step — the Phase 4 leakage rule). It accumulates
  `batch_size` per-window losses per step (windows vary in length; the fast
  kernel will tensor-collate), clips the gradient, and writes `best.pt` plus
  `run_manifest.json` (commit, seed, hardware, torch/CUDA versions, per-species
  pinned gff/fasta MD5s, param count = 455,841, sampled bases) to `out_dir`.
- `measure --config CONFIG.json --species S [--seqid ID] [--checkpoint best.pt]`
  times **preprocessing** (featurizer), **encoder** forward, and
  **decode/traceback** (the reference delayed-entry Viterbi over the emissions)
  separately with `time.process_time`, over the windows of one sequence, and
  reports oriented bases, CPU-s/Mb, peak host RSS and peak device memory in the
  `docs/cost-baseline` convention.

The config is a JSON `TrainConfig`: a `sources` list of
`{name, summary, gff, fasta, dev_seqids}` (the pinned paths the loader
checksum-verifies), plus `seed`, `steps`, `batch_size`, `lr`, `weight_decay`,
`grad_clip`, `max_window`, `eval_every`, `out_dir`, `device`.

Per the task's definition of done and proposal section 6, the gagarin runs will:

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

The first gagarin request (alert lenin-0083) is a bounded **verification and
profiling** run, not the full pilot: it runs the full candidate-A torch test
suite on a real torch host (the input-validation, parity and autograd tests that
skip locally) and a short profiling run of `train`/`measure` on the smallest
train species, to confirm the tensor path and measure the reference-recurrence
per-window cost before the full fit is scoped. The reference forward has the
reference recurrence's cost and is not the training kernel; the profiling
determines whether the fast delayed-entry kernel must land before the full pilot
fits within the 24 GPU-hour cap.

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
  breaking the symmetry under which the mixture-logit gradients are zero and the
  hazard-logit gradients are equal across components — the components cannot
  differentiate under symmetric updates (hazard gradients can be nonzero but move
  the components together; stalin-0063, engels-0060). Uniform mixture weights are
  kept. A test asserts the three components start distinct.

Chain-loss oracle review findings and their resolution:

- **The numerator constrained emissions but not boundary states** (engels-0062
  P2). The support mask masks emission channels, not the initial/terminal states,
  so an enabled `EdgePrior` let the numerator claim extra entry/exit and phase
  hypotheses — for a complete gene at a real edge `log Z_num` becomes `log(5/4)`
  instead of `0` — while the loss could stay nonnegative and pass the tests,
  supervising the wrong path set. This increment is the *complete*-target oracle:
  the gold path enters and leaves in intergenic `U`, so `chain_nll` now rejects a
  decoder whose edge grammar is enabled (`ValueError`), rather than silently
  mis-scoring it. Declared edge-partial numerators (compatible entry/exit families
  and a first-row phase carried alongside the CDS intervals) are the section-3.6
  boundary-support interface, still to be built. A test pins the rejection.
- **The torch-parity gate aborted discovery on a torch host** (engels-0062 P2).
  `tests/test_a_loss.py` caught only `ModuleNotFoundError`, but with torch present
  `from model.a import torch_loss` raises a plain `ImportError` because the
  submodule does not exist yet, aborting test collection instead of skipping. The
  gate now probes `importlib.util.find_spec("model.a.torch_loss")`, skipping
  cleanly when the submodule is absent while still surfacing errors from a
  genuinely broken implementation once it is written.

Training-window loader review findings and their resolution (engels-0064,
stalin-0067):

- **Flanks forced false intergenic labels on neighbouring genes** (P2). A
  window's flank could contain another gene's CDS, which `numerator_scores`
  then forced to intergenic `U`; a masked neighbour must likewise stay
  unconstrained (section 3.6). `iter_windows` now yields only *clean* windows:
  a window is skipped (and counted in `LoaderStats.skipped_neighbor`) when any
  transcript of a **different gene** (`gid`) overlaps it. Isoforms of the same
  gene do not make a window dirty. Adjacent-gene support built from the full
  window's annotations is the later section-3.6 increment. Regression:
  `test_neighbouring_gene_skips_window` (two genes at 11..19 and 21..29 both
  skipped).
- **Orientation destroyed the soft-mask channel** (P2). The reused
  `oriented_chain` upper-cases its window, so the featurizer's channel 5 saw no
  soft masking. The window is now rebuilt case-preserving from the raw FASTA
  slice (`_oriented_window`; `revcomp` preserves case), with `oriented_chain`
  kept as the authority for the CDS/intron coordinates and a guard that the two
  agree up to case. Regressions: `test_soft_mask_channel_preserved_plus` (an
  all-lower-case contig yields an all-lower-case window) and
  `test_oriented_window_preserves_case_both_strands`.
- **The hard source gate was optional and settings were unbound** (P2).
  `iter_windows` now takes the species `*.summary.json` as its first argument,
  calls `verify_source` before yielding anything, and reads the audit `m`/`table`
  from that summary rather than from loader defaults, so it cannot audit unpinned
  inputs or use settings that disagree with the committed manifest. Regressions:
  `test_iter_windows_enforces_source_pin`, `test_table_is_bound_from_summary`.
- **`import model.a` dragged in `model.labels`/`benchmark`** (P2). Added
  `model.labels` to `[tool.setuptools] packages`, and made `model.a.__init__`
  import the loader lazily (PEP 562 `__getattr__`), so `import model.a` — and the
  torch-free parameter-count guard, the encoder and the loss — no longer require
  `model.labels` or `benchmark/score.py`. The loader is a checkout-time tool (the
  admission audit loads `benchmark/score.py` by path); its dependencies are only
  imported when a loader symbol is accessed. Verified that `import model.a`
  imports no `model.labels` submodule.

Differentiable (torch) chain-loss review findings and their resolution
(engels-0071, stalin-0074):

- **The torch entry points bypassed the scalar oracle's emission contract** (P2).
  `ReferenceDecoder._check_input` requires every emission channel to be finite or
  `-inf` (a hard mask); the torch `_partition` instead treats `torch.isinf` as
  "drop this path", so a stray `+inf` silently discards a legal transition and a
  `NaN` poisons `logsumexp` — e.g. `ATGTAA` with zero scores and `start[0]=+inf`
  returned `log Z = 0` after the gene path vanished, and `CCCCCC` with an unused
  `donor[0]=NaN` returned `0`, where the oracle raises `ValueError`. Added
  `_check_input(x, emissions)`, called from both public entry points (`partition`
  and `chain_nll`), which rejects any `NaN` or `+inf` on any channel (used or
  unused) while preserving legitimate `-inf` support.
- **`partition` did not check the sequence length** (P2). `partition` passed the
  tensor straight to a recurrence bounded by the emission width, so `ATGTAA` with
  three emission columns silently scored the prefix and seven columns raised
  `IndexError`; the oracle raises `ValueError` for either mismatch. The same
  `_check_input` enforces `emissions.shape == (11, len(x))` on both entry points
  (`chain_nll` already checked this dimension; it now shares the one contract).
- Regressions in `tests/test_a_torch_loss.py` (`InputValidation`, torch-gated):
  `+inf`, `NaN` on a used channel and `NaN` on an unused channel are each rejected
  through both entry points; `-inf` support is accepted; short (three-column) and
  long (seven-column) emissions are rejected through both entry points. These are
  autograd/tensor checks and require a torch host (skipped on the Python 3.14.4
  torch-less host here), as the reviewers noted.

Training/measurement entry-point review findings and their resolution
(engels-0073, stalin-0076):

- **The pooled decoder was disconnected from the loss but put in Adam** (P1,
  engels-0073). The fixed-grammar chain loss reads only `model.encoder`
  emissions and scores them through `model.grammar`'s default duration/motif
  model; `model.decoder`'s 54 pooled scalars therefore have no gradient path.
  Rather than claim to fit them, this increment is now explicitly scoped as an
  **encoder-only profiling fit against the fixed grammar**: `train` builds the
  optimizer over `model.encoder.parameters()` only, the module/manifest carry
  `scope: "encoder-only-fixed-grammar"`, and `measure` documents that its
  decoder is the fixed-grammar reference. Wiring the learned pooled
  duration/motif model into both the differentiable loss and the decoder — with
  a nonzero-gradient fixture and a checkpoint-perturbation test — is the next
  increment, kept distinct from this profiling run.
- **Measurement promised GPU time from CPU clocks** (P2, engels-0073). `measure`
  used `time.process_time()` for every stage, which does not see asynchronous
  CUDA kernels. Each stage now records **both** process CPU seconds and
  CUDA-synchronized elapsed wall seconds (`preprocess_wall_s`/`encoder_wall_s`/
  `decode_wall_s`); `gpu_s_per_mb` is taken from the wall clock and is `None`
  off `cuda`, `cpu_s_per_mb` from process CPU. The row is labelled
  `profile: "annotation-selected-windows"` with `outputs_discarded: true`,
  preserving in the returned artifact that it profiles clean admitted gene
  windows, not a full-chromosome sensitivity/budget row.
- **The manifest omitted the split/schedule and was built after fitting** (P2,
  engels-0073). `build_manifest` now records per-source `name` and `dev_seqids`,
  `eval_every`, and a `sampling_plan` (seed, steps, batch_size,
  `planned_draws = steps*batch_size`, max_window), and is written **before** the
  optimizer loop as the declared plan; the collision engels reproduced (differing
  reserved chromosomes / eval schedule, identical manifest) no longer occurs.
  Actual attempted/accepted work (`attempted_draws`, `accepted_windows`,
  `sampled_bases`, window counts, best dev NLL) is attached afterward by
  `record_actual` and rewritten into the same file. Regressions:
  `test_prefit_manifest_declares_plan_and_provenance`,
  `test_dev_seqids_and_eval_every_break_collision`,
  `test_record_actual_appends_executed_work`.
- **Malformed dev reservations were silently accepted** (P2, stalin-0076).
  `SpeciesSource.from_dict` coerced `dev_seqids` with `list(...)` and ignored
  unknown keys, so `"dev_seqids": "chrDev"` reserved six one-letter ids, a
  `"dev_seqid"` typo reserved nothing, and an unknown id matched nothing — each
  leaving the intended development sequence in the gradient pool with a silent
  fallback to the final checkpoint. `from_dict` now rejects unknown source keys
  and any non-list `dev_seqids`; `validate_dev_reservations` (run before the
  first gradient step) rejects a declared id that produced no admitted window,
  using the loader's per-species `windows_by_seqid` inventory. Regressions:
  `test_source_rejects_unknown_key`, `test_source_rejects_string_dev_seqids`,
  and `TestReservations`. The yeast smoke config's intentionally empty dev list
  stays a valid profiling mode, distinct from a development-selected fit.

## 7. Training-set coverage accounting

The complete-target, clean-window scope of section 3.6 is temporary: it drops
edge-partial admitted chains, chains a neighbouring gene overlaps, and (when a
`max_window` is set) chains longer than the chosen chunk. Both incremental
reviews (engels-0065, stalin-0068) asked that this exclusion be quantified
**before any train-panel conclusion is drawn from A**, because the fraction of
the admitted set the loss actually sees bounds what the pilot's accuracy means.

`model.a.coverage_report(sources, *, complete_only, max_window)` runs the loader
over each `(species, summary, gff, fasta)` and returns a `CoverageRow` per
species — `admitted`, `yielded`, `yielded_fraction`, and the three skip counts
(`skipped_partial`, `skipped_neighbor`, `skipped_too_long`) that partition the
difference — under the same checksum gate and audit delegation as training, so
the numbers are exactly what a fitting run would train on. `format_coverage`
renders the rows as a TSV with a `TOTAL` row, following the
`docs/cost-baseline/measured.tsv` convention. `tests/test_a_dataset.py` checks
the aggregation, the `max_window` accounting, the source-pin enforcement, and
the rendered header/TOTAL on synthetic species.

**Producing the panel table is one command.** The raw source GFF3/FASTA are too
large to commit (charter: no files over 5 MB), so running the report needs a
source checkout. `model.a.coverage` makes that turnkey:

```
# One small species on a laptop (yeast, worm or Dictyostelium):
python3 -m model.a.coverage --sources <scratch-dir> \
    --fetch --species Saccharomyces_cerevisiae

# The whole train panel — the mammal/maize genomes and their audits need
# gagarin's memory, so run this as the named compute request, not on a laptop:
python3 -m model.a.coverage --sources <scratch-dir> --fetch --out coverage.tsv
```

Passing no `--species` selects all ten committed train species (including
mouse, maize and zebrafish), so the unfiltered command is the gagarin
invocation, not a laptop one. An unknown `--species`, or an empty manifest
directory, is now a clean argument error in both modes rather than a
zero-row table, and `--fetch` progress goes to stderr so a redirected stdout
is a clean TSV.

`--fetch` constructs each pinned source's deterministic NCBI `genomes/all` URL
from the filename in its `*.summary.json` (`ncbi_url`, no scraping or guessing),
downloads the exact `_genomic.{gff,fna}.gz`, and verifies each against the
summary's MD5 — the same digest `dataset.verify_source` re-checks before any
window is yielded, so a wrong-assembly or corrupted download can never reach the
loader. The report then delegates admission to `audit_species` exactly as
training does. `tests/test_a_coverage.py` covers the URL construction (including
the real committed manifests and assembly names with underscores), the resumable
MD5-gated fetch (stubbed, offline), and the end-to-end CLI TSV on a synthetic
species; `python3 -m model.a.coverage --self-test --sources x` runs the offline
URL checks alone.

Yeast, worm and Dictyostelium fit a laptop; the mammal and maize genomes and
their audits (peak RSS well past a laptop for the largest, engels-0050/0052) are
the runnable command a gagarin compute request names.

This measures the current scope; it does not widen it. The four small train
species were checked out and reported on the laptop (lenin, 2026-09-19,
commit ef82af9, `max_window=None`); Arabidopsis and Drosophila — the two
smallest of the remaining six — were added on the laptop the same way (lenin,
2026-09-20; MD5-verified fetches, ~60 MB and ~53 MB compressed, 27 s and 13 s
audits, nothing committed). The four large genomes (Danio 1.45 Gb, Xenopus
1.45 Gb, Mus 2.7 Gb, Zea 2.18 Gb) were also completed on the laptop (lenin,
2026-09-20): each source `_genomic.{gff,fna}.gz` was fetched to scratch and
MD5-verified against its committed manifest (`fetch` re-checks the same digest
`verify_source` enforces before any window is yielded), audited within the
per-species peak RSS recorded in the manifests (1.9–3.6 GB) on one core, and
its FASTA discarded; nothing committed. The panel table below is therefore the
full ten-species train set. A yielded fraction low enough to bias the pilot is
itself a section-2 finding to report before the gagarin fitting run, not after.

Reading the measured rows: the complete-target clean-window scope yields
**173,216 of 179,227 admitted representatives, 96.6% panel-wide**, and
91.8–99.9% per species across all ten train species. **`skip_too_long`
is 0 everywhere**, but only because the reported run uses `max_window=None`
(the CLI, `iter_windows` and `coverage_report` default): with no finite length
limit the length filter never fires, so this column measures unrestricted
whole-gene loading, not fit within the 3,072-base encoder core or the cost of
the pending crop. The exclusion here is therefore neighbour-overlap
(6,010 windows) plus a single edge-partial, not a gene-length cap (consistent
with engels-0047, and see the crop still listed pending in section 2).
Neighbour skips are largest by count in the vertebrate Danio (1551, 7.9%) and
by fraction in the compact Drosophila (1069, 8.2%) and C. elegans (1090, 5.5%)
genomes, and smallest in Neurospora (8, 0.1%) — the same gene-density pattern
holds across the four large genomes (Xenopus 2.4%, Mus 3.4%, Zea 1.2%), with
`skip_too_long` still 0.

These are retention counts, not a bias measurement. The loader excludes windows
because a different gene overlaps them, so the exclusion selects on locus
architecture, not at random: stalin-0072 certifies that a concrete
opposite-strand CDS-overlap stratum is removed on every species (92/118 yeast,
957/1090 worm targets), and engels-0069 bounds the long-gene spans that must
survive. High overall retention can still coexist with systematic removal of a
gene stratum, and this run does not compare retained versus excluded targets by
length, intron count, locus type or development chromosome, nor under the final
sampling/repeat policy. So the effect of the clean-window exclusion on the pilot
is **unassessed**, and the gene-density explanation above is a hypothesis
pending an actual density/spacing analysis — not a demonstration that the pilot
is unbiased. Edge-partials are negligible across the whole panel (1 window
total, in Dictyostelium), so the section-2 boundary-support increment will move
the panel yield fraction little; its value is on scaffold-edge structure, not
this retention count.

| species | admitted | yielded | yielded_pct | skip_partial | skip_neighbor | skip_too_long |
|---|---:|---:|---:|---:|---:|---:|
| Saccharomyces_cerevisiae | 5858 | 5740 | 98.0 | 0 | 118 | 0 |
| Dictyostelium_discoideum | 12937 | 12886 | 99.6 | 1 | 50 | 0 |
| Neurospora_crassa | 9722 | 9714 | 99.9 | 0 | 8 | 0 |
| Caenorhabditis_elegans | 19784 | 18694 | 94.5 | 0 | 1090 | 0 |
| Arabidopsis_thaliana | 27220 | 26630 | 97.8 | 0 | 590 | 0 |
| Drosophila_melanogaster | 12974 | 11905 | 91.8 | 0 | 1069 | 0 |
| Danio_rerio | 19574 | 18023 | 92.1 | 0 | 1551 | 0 |
| Xenopus_tropicalis | 17660 | 17230 | 97.6 | 0 | 430 | 0 |
| Mus_musculus | 21694 | 20956 | 96.6 | 0 | 738 | 0 |
| Zea_mays | 31804 | 31438 | 98.8 | 0 | 366 | 0 |
| **TOTAL** | **179227** | **173216** | **96.6** | **1** | **6010** | **0** |
