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
  standard-library reference the PyTorch loss must match on the same fixtures.
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
2. **Differentiable (PyTorch) chain loss**: the torch numerator/denominator
   forward pass whose autograd yields `dL/de = posterior_free − posterior_num`,
   matching the `model/a/loss.py` oracle on the section 3.4 fixtures (the
   torch-gated `TorchParity` test in `tests/test_a_loss.py` is skipped until it
   lands). The support-mask construction, the constrained-vs-free contract and
   the specification-oracle cross-check against the reference decoder are done;
   what remains is the batched, chunk-seam torch implementation on the encoder
   emissions and the boundary conditioning at training-crop edges (section 3.6).
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
python3 -m model.a.coverage --sources <scratch-dir> --fetch          # small species on a laptop
python3 -m model.a.coverage --sources <scratch-dir> --out coverage.tsv
```

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
the runnable command a gagarin compute request names, which is why the table
below stays pending a source checkout rather than being filled from a partial or
approximate pass.

This measures the current scope; it does not widen it. The panel-wide numbers
(the yielded fraction per train species and the skip breakdown) are recorded
here once the command above is run over the checked-out train sources, alongside
the edge-partial and crop increments of section 2 that will raise the yielded
fraction. A yielded fraction low enough to bias the pilot is itself a section-2
finding to report before the gagarin fitting run, not after.

| species | admitted | yielded | yielded_pct | skip_partial | skip_neighbor | skip_too_long |
|---|---:|---:|---:|---:|---:|---:|
| _pending gagarin/source checkout_ | | | | | | |
