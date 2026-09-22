# Candidate A pilot: implementation and end-to-end measurement (T-human-014)

Status: **in progress** (owner `lenin`). This document is the T-human-014
deliverable and is revised in place each tick. The tensor path is verified,
the fast delayed-entry training kernel implemented and timed, and the CPU
regime is measured end to end with a frozen smoke checkpoint on the
*S. pombe* normalization row and on the metazoan development chromosome
(*C. elegans* chr V, 20.9 Mb), both at 19 segments per strand: **9.2–9.7
CPU-s per genome Mb, 1/5.4–1/5.5 of AUGUSTUS on the same machine, missing
the 1/11 portable target by 2.0× (1.8× on user time); no positive CPU
allowance for B** (section 3.2). The first fitted checkpoint (bounded
local CPU fit v1, section 3.3) scores the two development chromosomes at
nucleotide sensitivity 0.04 / 0.47 with 0 exact transcripts (287 chains
of median 81 kb fusing 238 loci on chr V). The loader revision
(fit v2, section 3.3: `U` context added to every window, 30.6 % rather
than 3.1 % of the sampled bases supervised as `U`) sharply reduces the
fusions — chr I nucleotide F1 0.83 with 69 of 94 transcripts exact,
chr V 4 fusions (from 238) and 3,804 of 4,995 loci hit — but predicts almost no real intron (chr V donor F1 0.014,
median predicted intron 34 bases, 1,827 splits) and over-predicts short
single-exon chains (9,606 unmatched loci on chr V, a third of them
fragments of genes it already hits): **A still misses the
accuracy target in this fit**. The fusion reduction supports the
supervision-balance hypothesis for v1 without establishing that the
defect is fully removed (3.3); the misses that remain are splicing and
precision, with their cause unresolved. The fitted rows are 1.9–2.7×
over the portable CPU ceiling. **Fit v3** (section 3.3), the longer fit proposed there — 3,000 steps
with a cosine rate and a 150-step warmup, 4.33 CPU-h on one core, best
checkpoint at the final step, dev NLL 45.3 → 28.6 — changes the picture
again without accepting A. Nucleotide F1 rises to 0.891 on chr I and
0.598 on chr V and locus F1 to 0.786 / 0.603, entirely through precision
(0.742 → 0.870, 0.546 → 0.680) at flat sensitivity: half as many chains,
twice as long. The 20-base intron defect is gone (median predicted intron
34 b → 615 b, correct GT-AG introns 68 → 433) but splice-site placement
is not fixed (donor F1 0.070, acceptor 0.081, 22,187 reference GT-AG
introns still missed), and the split defect has become a fusion defect
(splits 1,827 → 910, fusions 4 → 312, chr V exact-transcript sensitivity
0.015 → 0.010). Cost is unchanged (12.38 / 8.88 CPU-s/Mb), so the CPU
verdict stands and **no positive CPU allowance for B follows**. The
S. pombe normalization row under the v2 checkpoint is 9.06 CPU-s/Mb
user + system, 1/5.7 of AUGUSTUS against the 1/11 target (section 3.3).
The revision proposed next is the **decoder**, with the section 6 CPU
ceiling in scope, since the remaining failures are decoding decisions;
the monotone dev tail also leaves a still longer fit open, and the two
would have to be measured separately. Still pending: that revision and
the GPU regime (gagarin, lenin-0083).

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
  Plus the 54 pooled-decoder scalars (`DecoderParams`), consumed through
  `model/a/pooled.py` (below).
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
- `model/a/fast_loss.py` — the **vectorized delayed-entry training kernel**
  (sections 3.2–3.3). The same `log Z − log Z_num` as `torch_loss.py`, computed
  with the delayed-entry recurrence of `DelayedEntryDecoder` on dense tensors:
  the coding layer is one `(B, K)` vector over the K = 24 coding states (U,
  the initiator prefixes S(q), the 21 ordinary prefixes E(q)); the intron
  tails one `(B, K, R)` matrix; a donor parked at boundary s is `alpha[s] +
  donor[s]` and enters every tail exactly m boundaries later with the rolling
  m-window sum of the phase's intron emissions (`unfold`) and `log pi[p, r]`,
  so the m−1 mandatory positions are never visited. The per-base coding
  transition `(n, K, K)` is one `einsum` of a fixed 0/1 channel-coefficient
  tensor with the emissions plus a symbol-indexed prior table (uniform over the
  permitted bases of an IUPAC code, summed over bases reaching the same
  successor). The scan is a Python loop of n steps over those small tensors,
  batched over windows (padded positions get the identity transition and
  floored intron/donor/acceptor emissions). `-inf` never enters the scan:
  emissions are clamped at a finite `FLOOR = −1e30`, so forbidden paths carry
  weight exactly 0 and zero posterior, and a window whose support admits no
  path is reported as `−inf`/`+inf` loss rather than NaN gradients.
  Complete targets only, like the reference. Entry points: `partition`,
  `chain_nll` (one window; the drop-in for `torch_loss.chain_nll`) and
  `batch_chain_nll` (padded `(B, 11, L)` batch). `tests/test_a_fast_loss.py`
  (torch-gated, 8 cases) pins **exact parity with the reference torch forward**
  — partition, loss and emission gradient to ≤ 1e-9 on the section-3.4
  fixtures, on 30 random short lattices with IUPAC ambiguity, both genetic-code
  tables, m ∈ {1..4}, R ∈ {1..3} and sparse `-inf` masks (infeasible lattices
  `-inf` in both), batched-vs-single equality with zero gradient on padding —
  plus `+inf` (not clamped) for an infeasible numerator and input rejection.
  Cost: section 3.2.
- `model/a/pooled.py` — the **learned pooled decoder** (proposal 3.5's 54
  scalars as decoder tables, 3.2 item 5). `duration_tables` turns the
  phase × component mixture and hazard logits into differentiable `(3, R)`
  `log π`, `log q`, `log (1 − q)` tables (`π = softmax`, `q = sigmoid`) that
  `fast_loss.log_partition_batch`/`chain_nll`/`batch_chain_nll` and
  `fast_viterbi.viterbi_batch`/`viterbi`/`viterbi_windows` take (`tables=`)
  in place of the fixed `DurationMixture` values; the mixture then only fixes
  the grammar's shape `(m, R)` (`structure`, a stable `Grammar` cache key —
  the tables are rebuilt per step, the grammar is not). `motif_bias` adds the
  16-entry donor and acceptor dinucleotide scores to the `donor`/`acceptor`
  emission rows — a donor at boundary `t` pays `donor_dinuc[x[t] x[t+1]]`, an
  acceptor `acceptor_dinuc[x[t−2] x[t−1]]`; non-ACGT pairs and pairs past the
  window score 0 — so both kernels and the Python reference decoders are
  unchanged. `as_mixture` gives the same law as a concrete `DurationMixture`
  for parity runs. `edge_prior` reads the four partial-family scalars as the
  `EdgePrior` of proposal 3.1 (coding entry / coding exit / intron entry /
  intron exit; the reference decoders' two-value prior is the case where the
  intron scalars equal the coding ones) for the edge-enabled decoders below —
  values only; their gradient waits on the section-3.6 edge-partial
  numerators in the loss. `tests/test_a_pooled.py`
  (torch-gated, 8 cases) pins table/mixture agreement, bias indexing, fast-loss
  and tensor-Viterbi parity with the reference kernels under the learned law
  (fixtures, 30 + 40 random lattices, batched = single), a finite gradient on
  all 50 consumed scalars, and central-difference agreement of the mixture,
  hazard and dinucleotide gradients.
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
   the edge-enabled decoder that now exists for Viterbi but not yet for the
   partition) which are currently skipped and counted — and
   cropping whole genes longer than the encoder core into chunks with retained
   intron/codon/duration state (the crop-integration contracts checked in the
   prior owner's notes engels-0045…engels-0057). Batched multi-window collation
   for the torch training step also belongs here. **Gene-free (background)
   windows are in** (proposal section 3, "gene-free/background windows"):
   `iter_background_windows` tiles every gene-free interval of a sequence
   that carries an admitted chain (every raw GFF3 gene, pseudogene,
   transcript, exon, CDS and UTR row on both strands blocked —
   `annotated_gene_spans`, read from the file rather than from the
   admission parser, whose transcripts are CDS-bearing only and whose
   spans stop at the CDS ends, so UTRs, ncRNAs and pseudogenes are blocked
   too (engels-0096); non-gene features such as centromeres, LTRs and
   replication origins stay eligible — widened by `FLANK`; a mitochondrion
   or any sequence the audit admitted nothing from is never background)
   with non-overlapping
   `background_length` tiles and draws `background_windows` of them per
   species without replacement under the run seed, half of them
   reverse-complemented with case preserved. Their support is the empty
   chain — `numerator_scores`/`support_mask` now accept it as the
   all-intergenic path (introns without a CDS stay rejected) — so the fit
   also sees windows where the right answer is no gene, which the
   chain-only smoke fits never did (the smoke checkpoint's 208,796 chains
   on chr V's 3,357 genes). Both keys are in the config, manifest
   hyperparameters and sampling plan; the default is 0 (chain-only, the
   scope of every recorded fit), and a non-null `max_window` bounds
   background windows as well as chains (an enabled `background_length`
   above it is refused at config parse and again in the loader). The
   tiling is coarse on a gene-dense genome: *S. cerevisiae* has 14
   gene-free 4,096-base tiles, 148 at 2,048, 766 at 1,024 (zero overlap
   any raw gene or pseudogene row; the earlier CDS-span-free counts of
   23 / 198 / 928 included tiles inside UTRs and non-coding genes), so
   the fitted-checkpoint config sets the length per species panel. Background windows are counted apart from the admitted
   inventory (`LoaderStats.background_candidates`/`background_yielded`),
   so dev-seqid reservations still key off admitted chains only.
2. **Differentiable (PyTorch) chain loss** — the reference forward
   (`model/a/torch_loss.py`) and the **fast vectorized delayed-entry kernel**
   (`model/a/fast_loss.py`, checked against it to 1e-9; section 3.2) are done,
   and `train.py` fits with the fast kernel by default (`loss_kernel: "fast"`
   in the config and manifest; `"reference"` keeps the expanded forward for
   parity runs). What remains is chunk-seam handling for crops longer than
   the encoder core (the `Checkpoint` semantics of `DelayedEntryDecoder`, on
   tensors), the **boundary conditioning at training-crop edges** (section
   3.6, edge-partial numerators), and collating the batched entry point
   `batch_chain_nll` into the training step (the loop still accumulates one
   window at a time).
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
4. **Bring the chromosome row inside the CPU budget.** The first
   chromosome-level `measure` row exists (3.2 item 6: both strands,
   overlapping windows, GFF3 output, I/O). After revision steps 1
   (vectorised featurizer), 2a (per-phase decode operands) and 2b
   (sparse-predecessor scan) it is 17.1 CPU-s per genome Mb at the defaults
   (stage sum; 18.1 whole-process user, the `cost-baseline.md` convention),
   **14.2 / 14.4 at the default overlap with batch 64 — the first
   overlapping row inside the 15 budget on the baseline's numerator** —
   12.0 / 12.8 at overlap 2,048 / batch 64 and 10.7 / 11.9 at overlap 0 /
   batch 64; the batch-16 default still misses by 1.14–1.2×. Revision
   step 3 (carried scan state across tile seams, proposal 3.3) is
   measured: 19 segments per strand with a 4,096-base seam overlap gives
   **11.8 / 12.3** (11.9 / 12.4 with the packed back-pointer store) with
   every interior tile seam exact and the 2,048-base chain-containment
   guarantee kept at segment seams, and its oversampling shrinks with
   chromosome length (1.32× on yeast chr I, 1.004× at 20 Mb) where the
   window mode's 1.5× is fixed; the back-pointers are held packed at
   23 bytes per base (120 dense) and the emissions are streamed into the
   scan one tile at a time (**12.3 / 12.6 / 14.9**, peak RSS 0.96 GiB;
   the exact strand decode 0.35 GiB), so the state that grows with the
   chromosome is the packed store alone: the 20 Mb metazoan chromosome
   projects to ~1.5 GB of live tensors (0.92 GB pointers, 0.40 GB per-tile
   operands, 0.13 GB one row's dense expansion, 0.04 GB tile emissions)
   against 8 GB — a source-derived projection, not a measurement; the
   earlier "~3.5 GB" figure for the padded path was wrong (engels-0090,
   stalin-0093: float64 emissions held twice put that path at ≥ 8.7 GB).
   Each tile is now encoded with the encoder's dependency radius (491
   bases) of chromosome context on each side and cropped back, so the
   emissions no longer depend on the tile grid (**13.2 / 13.3 / 15.9** at
   19 segments per strand, encoder +12 % on these short tiles, +8 % on
   full ones; peak RSS unchanged) and the segment rows differ from the
   exact strand decode by 4 chains of 1,832 (were 11), all four genes
   longer than the 2,048-base containment guarantee starting just before
   a segment core seam — the seam contract, not a defect. The S. pombe
   normalization row is measured (section 3.2): 9.22 / 8.23 / 9.37 CPU-s
   per genome Mb at 19 segments against AUGUSTUS 51.4 on the same core,
   **1/5.5 against the 1/11 target — the CPU regime misses by 2.0×**, no
   B allowance. Revision step 4, float32 decode, is measured: 8.27 /
   7.94 / 8.43, a 10 % gain, not the halving of the decode that was
   expected, because the scan is per-step dispatch-bound; it flips
   103 near-tie chains of 121,937 at the smoke checkpoint, so float64
   stays the row of record. The remaining decode lever is the per-step
   dispatch (a compiled scan step), and the encoder halving waits for a
   fitted checkpoint. The learned pooled
   decoder is in the loss and the decoder (`model/a/pooled.py`, 3.2 item 5);
   what remains of the decoder is the partial-family scalars, which belong
   to the boundary-support increment above.

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
`grad_clip`, `max_window`, `eval_every`, `out_dir`, `device`, `loss_kernel`,
`min_intron`, `background_windows`, `background_length` and
`dev_windows_max` (checkpoint selection scores at most this many
development windows, a seeded draw without replacement from the reserved
chromosomes' windows; the manifest records both `dev_windows` and
`dev_windows_evaluated`). The manifest's `actual` block also carries
`best_step`, the per-evaluation `history` (step, train NLL, dev NLL,
elapsed) and `timing` (load / fit / eval wall and process-CPU seconds,
thread count), so a fitted run's cost is in the artifact itself.

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

### 3.1 Local smoke run, 2026-09-20 (CPU only; the gagarin request stays open for the GPU half)

`uv` became available on the dev host, so a Python 3.11.14 venv with
`torch 2.14.0+cpu` (no CUDA) was built and `model/a/gagarin_smoke.sh` was run
unchanged at commit `dfa9df9` on an Intel Core Ultra 9 285K (24 logical CPUs,
62 GB RAM). Raw outputs (test log, `run_manifest.json`, `coverage_yeast.tsv`,
`measure_yeast.json`, both `/usr/bin/time -v` reports, per-window profile) are
in `relay/artifacts/T-human-014/smoke-local-20260920/`. This covers everything
in alert lenin-0083 except the CUDA build, device memory and the GPU regime,
which still need gagarin.

Memory convention for 3.1 and 3.2 (engels-0080 P3): every figure is the process
high-water RSS (`ru_maxrss`, KiB on Linux; the profile scripts print MiB), stated
in **GiB** (2^30 bytes) with decimal GB in parentheses where a number is near the
cap. The charter's "8 GB" threshold is kept as written, 8,000,000,000 bytes
(7.451 GiB); GiB is only the reporting unit, not a looser cap (engels-0081 P3).
Every figure here is on the same side of both readings. The profile-script RSS is cumulative across successive windows in
one process (model, loader, runtime and allocator included), so a per-base
figure derived from it is a process-peak increment, not bytes of autograd state.

**1. Tensor path verified.** Full suite under Python 3.11 + torch:
`pytest tests` → **133 passed, 0 skipped** (6.4 s wall, 599 MiB peak RSS); the
five candidate-A modules via the script's `unittest` line → 75 tests OK. This is
the first execution of the input-validation, local-attention oracle-parity,
`chain_nll` autograd and duration-init tests that every review to date listed
as unverified; all pass. The one runtime warning, a `float(loss)` on a grad-requiring scalar in
`train.py`, is silenced with `.detach()` in this revision (no behaviour change).

**2. Yeast fetch + coverage.** Both pinned sources MD5-verified from NCBI
`genomes/all`. 5,858 admitted / 5,740 yielded (98.0 %), 118 skipped for
neighbour overlap, 0 partial, 0 too long — identical to the committed
section-7 row.

**3. Encoder-only smoke fit (20 steps, batch 4, `max_window` 12288, seed 0):**

| quantity | value |
|---|---|
| windows drawn / accepted | 80 / 80 |
| sampled bases | 122,528 |
| wall clock | 12 min 50 s |
| CPU time (user+sys) | 1,634.8 s at 212 % CPU (torch default threads) |
| CPU-s per sampled kb | **13.3** |
| peak host RSS | **6.82 GiB** (7,149,860 KiB; 7.32 GB) |
| train NLL/window at step 10 / 20 | 0.0001 / 0.0006 (no dev split declared, so `best.pt` is the final state) |

**4. Per-window cost of the reference-recurrence chain loss** (single thread,
`torch.set_num_threads(1)`, float64, fresh `CandidateA`, one window each;
`profile_window.py` in the artifact directory):

| window bases | encoder fwd | loss fwd (enc + chain) | backward | total | s per kb | process peak RSS (cumulative) |
|---:|---:|---:|---:|---:|---:|---:|
| 422 | 0.006 s | 0.76 s | 0.83 s | 1.6 s | 3.8 | 0.54 GiB (557 MiB) |
| 1,223 | 0.009 s | 2.58 s | 3.37 s | 5.9 s | 4.9 | 1.22 GiB (1,245 MiB) |
| 2,804 | 0.020 s | 6.78 s | 9.25 s | 16.0 s | 5.7 | 2.62 GiB (2,682 MiB) |
| 11,255 | 0.090 s | 26.1 s | 55.7 s | 81.7 s | 7.3 | **8.79 GiB** (9,001 MiB; 9.44 GB) |

Reading: the encoder itself costs ~8 µs/base on one core (0.008 CPU-s/kb, i.e.
~8 CPU-s/Mb before decoding — inside the 15 CPU-s/Mb budget with room for the
decoder). The **reference chain loss is 500–900× the encoder**: 4–7 CPU-s/kb,
superlinear in window length, and the process peak grows by ~0.77 MiB per base
of the longest window (8.79 GiB less the 285 MiB post-load baseline, over 11,255
bases; a process-RSS increment, not an autograd byte count), so a single ~11 kb
window already exceeds the 8 GB cap (8.79 GiB = 9.44 GB, over it under either
unit). At 4–7 CPU-s/kb the reference loss is 4,000–7,000 CPU-s per Mb per pass on
one core. What the GPU regime would cost is not derivable from these CPU
measurements and stays pending with the gagarin request; the measured memory
breach and the CPU cost on their own are the reason to replace the reference
recurrence in training. **Conclusion: the fast vectorized delayed-entry kernel
(section 2, item 2) is mandatory before the pilot fit, not optional.** The
reference forward keeps its role as the differentiable oracle the kernel is
checked against.

**5. `measure` on yeast** (all 5,740 admitted windows, ~8.6 Mb oriented,
`best.pt` from the smoke fit, CPU) was started as the script's step 2d but had
not finished after 25 min wall (~3 CPU-h across torch's default 10+ threads):
the script-level run iterates the pure-Python reference delayed-entry Viterbi
over every admitted window of the genome, which the per-window profile above
did not time. Next tick re-runs it per sequence (`--seqid`, one chromosome)
with `torch.set_num_threads(1)` so the CPU-regime row is a one-core number, and
fills section 5. As the docstring says, that is an annotation-selected window
profile (outputs discarded, no both-strand/I/O/output accounting), not a
chromosome row.

**Compute recorded:** local CPU only — 0.45 CPU-h for the smoke fit, ~0.05
CPU-h for the per-window profile, ~3 CPU-h for the abandoned full-genome
`measure`; cluster CPU-hours 0, GPU-hours 0. No held-out species touched.

### 3.2 Fast delayed-entry kernel and one-core `measure` row, 2026-09-20 (CPU only)

Same host, venv and yeast data as 3.1; raw outputs in
`relay/artifacts/T-human-014/smoke-local-20260920/{fast-kernel,measure-chrI}/`.

**1. Kernel cost, like for like.** `profile_window_fast.py` is 3.1's
`profile_window.py` with `loss_kernel="fast"` (single thread, float64, the same
four windows):

| window bases | encoder fwd | loss fwd (enc + chain) | backward | total | s per kb | process peak RSS (cumulative) | vs. reference (3.1) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 422 | 0.007 s | 0.044 s | 0.047 s | 0.09 s | 0.22 | 0.33 GiB (342 MiB) | 17× faster |
| 1,223 | 0.009 s | 0.117 s | 0.114 s | 0.23 s | 0.19 | 0.40 GiB (409 MiB) | 26× |
| 2,804 | 0.021 s | 0.265 s | 0.275 s | 0.54 s | 0.19 | 0.53 GiB (542 MiB) | 30× |
| 11,255 | 0.111 s | 1.122 s | 1.217 s | 2.34 s | 0.21 | **1.08 GiB** (1,101 MiB; 1.15 GB) | **35× faster, 8.2× lower process peak** |

The cost is now **linear** in window length at ~0.2 CPU-s/kb (forward +
backward, one thread) and the process peak grows by ~74 KiB per base of the
longest window (same baseline-subtracted convention as 3.1), against the
reference's superlinear 4–7 s/kb and ~0.77 MiB/base; the longest yeast window
peaks at 1.08 GiB instead of breaching the 8 GB cap. Batched over four 12,288-base
windows the per-base cost falls further to 0.095 ms/base (float64) and
0.075 ms/base (float32; loss agrees with float64 to 2e-6 relative, gradient to
3e-3 absolute on unit-scale gradients — the fit keeps float64 for now). The
loss is ~25× the encoder's own ~8 µs/base, so the per-step cost of the pilot
fit is set by the scan's Python loop of n small tensor steps, not by the
network; a fused/parallel-scan implementation is the next lever if the pilot
step count needs it, not a blocker.

One implementation note worth recording: the first version of the scan
indexed the precomputed `(B, L, K, K)` transition stack as `trans[:, t]`
inside the loop, whose backward scatters into a full-size zero gradient at
every step and made the backward pass quadratic in L (4.9 s at 2 kb, minutes
at 12 kb). Slicing all per-step inputs once with `unbind` (backward = one
`stack`) restored linear cost; `tests/test_a_fast_loss.py` does not time
this, so the note is here.

**2. Smoke fit repeated with the fast kernel** (identical config, seed, data
and commit apart from `loss_kernel`; `fast-kernel/run_manifest.json`):

| quantity | reference kernel (3.1) | fast kernel |
|---|---|---|
| windows drawn / accepted | 80 / 80 | 80 / 80 |
| sampled bases | 122,528 | 122,528 |
| train NLL/window at step 10 / 20 | 0.0001 / 0.0006 | **0.0001 / 0.0006** (identical) |
| wall clock, torch default threads | 12 min 50 s | **1 min 32 s** |
| CPU time (user+sys), default threads | 1,634.8 s (212 %) | 1,492.6 s (1,627 %) |
| peak host RSS | 6.82 GiB | **0.92 GiB** (968,836 KiB) |
| wall / CPU time, one thread pinned to one core (`taskset`) | — | **24.8 s / 24.7 CPU-s** (0.20 CPU-s per sampled kb, 66× less CPU than the reference run; 0.87 GiB) |
| same, re-run at the committed source `129dcbe` (`fast-kernel/pinned-129dcbe/`) | — | **23.3 s / 23.3 CPU-s**, 0.87 GiB (914,560 KiB), NLL 0.0001 / 0.0006 |

**Provenance of the fast-kernel run** (engels-0080 P2). The two fast-kernel fits
and the fast profile in `fast-kernel/` were executed from the working tree
*before* the kernel was committed: their manifests record `commit: 751700b`
(the parent) although that commit has neither `model/a/fast_loss.py` nor the
`loss_kernel` field, so checking out the recorded SHA cannot reproduce them.
The raw manifests are kept unchanged and
`fast-kernel/PROVENANCE.md` reconciles them: parent `751700b`, dirty tree, the
source later committed as `28c5339` in the same session with no intervening
edit recorded but not verifiable by hash — so those three records lack an exact
source pin. The one-core fit was therefore repeated at the committed, clean
source (`129dcbe`, which also carries the cache-key fix below; manifest
`source_dirty: false`, `source_sha256` recorded) and reproduces the trajectory
and cost; that row is the pinned number. `train.py` manifests now always record
a SHA-256 over `model/a/*.py` + `model/grammar/*.py`, the dirty flag and the
dirty file list next to `commit`, so a future run from an uncommitted tree is
identifiable rather than mis-attributed.

The per-step NLL trajectory is identical to four decimals, i.e. the kernel
reproduces the reference loss and gradient on real windows, not only on the
test lattices. The default-thread CPU time barely moved because torch spreads
each of the scan's tiny per-step ops over all 24 logical CPUs and the
synchronisation dominates; the one-core run is the number that matters for
the CPU regime and the GPU run is what the batched scan is for.

**3. One-core `measure` on yeast chromosome I** (`--seqid NC_001133.9`,
`torch` pinned to one thread and one core with `taskset`, `best.pt` from the
smoke fit; `measure-chrI/measure_chrI.json` and the `/usr/bin/time -v`
report): 94 admitted windows, 143,653 oriented bases.

| stage | CPU s | CPU-s per Mb |
|---|---:|---:|
| preprocessing (featurizer) | 0.51 | 3.6 |
| encoder forward | 0.38 | **2.7** |
| decode + traceback (pure-Python `DelayedEntryDecoder.viterbi`) | 9.51 | **66.2** |
| total | 10.41 | **72.5** |
| peak host RSS | | 0.31 GiB (325,684 KiB) |

Against the 15 CPU-s/Mb budget: the encoder plus featurizer cost 6.3 CPU-s/Mb
on one core, inside the budget with room to spare; the **pure-Python Viterbi
decoder is the failing stage** at 66 CPU-s/Mb, 4.4× the whole budget on its
own. This is the same conclusion as 3.1 item 4 from the decoding side: the
decoder that ships in the CPU regime must be the tensor delayed-entry scan
(this kernel in max-product mode with back-pointers), not the standard-library
semantics check, whose job was always parity, not throughput. The
`measure` caveats of 3.1 item 5 still apply (annotation-selected windows,
outputs discarded, no both-strand/I/O/output accounting, encoder untrained
beyond 20 steps), so this is not yet a section-5 chromosome row; it is the
per-stage cost that the row will be built from.

**4. Tensor Viterbi replaces the Python decoder; one-core chromosome I
re-measured** (`model/a/fast_viterbi.py`, commit `80d82b2`, source SHA-256
`fc1a9e3a…51c96`, `source_dirty: false`; `measure-chrI-tensor/`). The decoder
is the max-product twin of the training kernel — same `n`-step loop over the
`(B, K)` coding layer and `(B, K, R)` tails, `max`/`argmax` for `logsumexp`,
int8 back-pointers into the coding layer, a bool `entered` map for the tails,
and a numpy traceback that re-expands the `m-1` mandatory intronic positions
so it returns the reference decoder's `Chain` objects. Two things had to differ
from the sum-product kernel: (a) the ambiguity prior — in Viterbi an IUPAC
symbol whose `k` permitted bases reach the same successor costs `-log k` once
(the reference's per-base uniform prior), not `log(count) - log k`, and `U→U`
carries no prior (`Grammar.prior_max`); (b) the dense `(B, L, K, K)` transition
stack is not materialized — at 16 windows of 12 kb it is 0.9 GB and was ~80 %
of the decode time — the transition is built per step from `cds[p(i)]`,
`prior_max[x[t]]`, the `U`-column extras (`u`, `stop`) and the `U→S(b)` post-add
(`start + cds[0]`). Parity with `DelayedEntryDecoder.viterbi` (score to 1e-9
and identical chains) is pinned by `tests/test_a_fast_viterbi.py` on the
section-3.4 fixtures and 120 random lattices (IUPAC, tables 1/6, alternative
initiators, `m` 1–4, `R` 1–3, `-inf` masks); length-bucketed batched decoding
(`viterbi_windows`) equals per-window decoding.

Same `measure` run as item 3 (chromosome I, 94 windows, 143,653 oriented
bases, `best.pt` from the smoke fit, one thread, `taskset` one core), now with
`--decoder tensor` (the default) and the windows decoded in length-sorted
batches of 16 — how a chromosome's windows would be decoded in either regime:

| stage | CPU s | CPU-s per Mb | vs item 3 |
|---|---:|---:|---:|
| preprocessing (featurizer) | 0.53 | 3.7 | — |
| encoder forward | 0.33 | 2.3 | — |
| decode + traceback (tensor scan, batch 16) | 0.71 | **4.9** | **13.4× faster** |
| total | 1.57 | **10.9** | 6.6× |
| peak host RSS | | 0.42 GiB (440,840 KiB) | +0.11 GiB |

Against the 15 CPU-s/Mb budget: **all three stages together are inside it on
one core**, at 73 % of the budget. The decoder is still the most expensive of
the three per base (4.9 > 3.7 > 2.3), but no longer dominant. Unbatched
(`--decode-batch 1`, `measure_chrI_tensor_b1.json`) the decode + traceback
stage alone is 30.0 CPU-s/Mb (three-stage total 35.9): the per-step cost is
torch dispatch, not arithmetic, so batching windows is what pays. On
synthetic 12 kb windows the decoder alone is 39 CPU-s/Mb at batch 1, 5.5 at 8,
3.8 at 16 and 2.9 at 32 (traceback ≤ 0.1 CPU-s/Mb of that). The caveats of
item 3 still apply — annotation-selected windows, outputs discarded, no
both-strand/I/O/output-writing accounting, encoder fit for 20 steps — so this
is still the per-stage cost a section-5 chromosome row is built from, not the
row; but the failing stage of item 3 is no longer failing.

**5. Learned pooled decoder in the loss and the decoder; chromosome I
re-measured with R = 3** (`model/a/pooled.py`, commit `9fb3aec`, source
SHA-256 `6ab52a79…e8b4`, `source_dirty: false`; `pooled-decoder/`). `train`
now carries all 455,841 parameters in the optimizer: the loss adds
`motif_bias` to the emissions and passes `duration_tables(model.decoder)` to
the fast kernel, so the mixture logits, hazard logits and both dinucleotide
tables get a gradient (`tests/test_a_pooled.py` checks it against central
differences); the four partial-family scalars are carried without a gradient
until the section-3.6 edge-partial numerators consume them in the loss (the
edge-enabled Viterbi of item 6's boundary-support step reads them as values),
so 455,837 parameters move. The selectable
reference kernel (`loss_kernel: "reference"`) at `9c6cd59` read the same law
through the detached `as_mixture`, so it fitted with the 18 duration scalars
silently frozen under the same manifest scope (engels-0083 P2, confirmed by
stalin-0084); it now reads a `TorchDuration` (0-d tensor entries of the same
tables, added by the reference recurrence exactly like its float priors), and
`TrainingKernelSelection` in `tests/test_a_pooled.py` asserts both kernels
give the same loss, emission gradient and decoder-scalar gradients through
the public `_window_loss` path. The fast default fit and the recorded
`pooled-decoder/` run are unaffected;
the optimizer and gradient clip cover `model.parameters()`, and the manifest
records `scope: encoder-and-pooled-decoder` and `min_intron` (`m`, default
20; the previous manifests' `encoder-only-fixed-grammar` remain as they
were). `measure` loads the checkpoint's decoder scalars, adds the bias (timed
in the decode stage) and decodes under the learned tables; the row records
`min_intron` and `pooled_decoder: learned`.

The 3.2 item 2 smoke fit repeated with the pooled decoder (same config, seed
and data, one thread, `taskset` one core): 26.9 s user / 27.8 s wall, peak RSS
0.91 GiB (955,656 KiB), train NLL/window 0.0001 → 0.0007 as before (item 2:
23.3 s at R = 1). All four consumed decoder tensors moved (max |Δ| 3–4 × 10⁻³
after 20 Adam steps at lr 3 × 10⁻⁴), `partial_families` did not. Then the
item 4 `measure` on chromosome I (94 windows, 143,653 oriented bases,
`--decoder tensor`, one thread, one core) with this checkpoint:

| stage | batch 16, CPU s | CPU-s per Mb | vs item 4 (R = 1) |
|---|---:|---:|---:|
| preprocessing (featurizer) | 0.51 | 3.6 | — |
| encoder forward | 0.35 | 2.5 | — |
| decode + traceback (tensor scan, R = 3, + bias) | 0.90 | **6.2** | +1.3 |
| total | 1.76 | **12.2** | +1.3 |
| peak host RSS | | 0.48 GiB (500,192 KiB) | +0.06 GiB |

Unbatched (`--decode-batch 1`, `measure_chrI_pooled_b1.json`) decode is 31.7
CPU-s/Mb, three-stage total 37.4. The +1.3 CPU-s/Mb is the proposal's
**R = 3** duration mixture (items 3–4 decoded the default single-component
`DurationMixture()`, R = 1, so their decoder rate understated the specified
grammar) — the `(B, K, R)` tail update triples — plus ~0.2 CPU-s/Mb for the
per-window dinucleotide lookup (Python indexing over the sequence; a tensor
gather if it ever matters). At 12.2 CPU-s/Mb the three stages are at 81 % of
the 15 CPU-s/Mb budget on one core, still inside. The item-3 caveats are
unchanged (annotation-selected windows, outputs discarded, no both-strand /
I/O / output accounting, 20-step encoder), and the learned scalars after 20
steps are numerically the initial ones to three decimals — this is the cost
of the specified decoder, not its accuracy.

**6. First chromosome-level row: yeast chromosome I end to end, both
strands, overlapping windows, GFF3 written** (`model/a/chromosome.py`,
`measure --profile chromosome`, commit `c6e5fc7`, `source_dirty: false`;
`chromosome/`). This is the row the items above were building toward: the
whole 230,218-base sequence, not annotation-selected windows. Each strand
(plus as read, minus as its reverse complement) is cut into 12,288-base
windows overlapping by 4,096 (step 8,192, so 1.5× oversampling; a complete
chain up to 2,048 bases is inside the window that owns it, with ≥ 2,048
bases of sequence before its start and ≥ 2,048 after its start — so a
chain of the maximal length can end on the window's last base; the
guarantee is containment, not a two-sided halo around the whole chain);
every window is featurized, encoded, given
the dinucleotide bias and decoded by the batched tensor Viterbi; a chain is
reported by the one window whose core holds its 5′ start (`tiles`,
`claimed`), mapped to genomic coordinates by the reviewed `gff3_rows`, and
written as GFF3. Timed stages: `io` (manifest digest check, FASTA read,
reverse complement), `preprocess`, `encoder`, `decode` (bias + scan +
traceback), `output` (map + serialisation + write). **Rates are per genome
megabase** (`genome_mb` = 0.230218, the `measured.tsv` convention: both
strands and the overlap are inside the numerator), which is what the
15 CPU-s/Mb budget is stated in; `cpu_s_per_oriented_mb` is also recorded
and is the number items 4–5 quoted. One thread, `taskset` one core, the
item 5 checkpoint (`/usr/bin/time -v` reports in `chrI_*_time.txt`):

| run | windows | oriented bases | io | preprocess | encoder | decode | output | **CPU-s per genome Mb** | CPU-s per oriented Mb | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| overlap 4,096, batch 16 (defaults) | 56 | 681,620 | 0.01 | 2.61 | 1.37 | 3.77 | 0.01 | **33.8** | 11.4 | 0.75 GiB |
| overlap 2,048, batch 16 | 46 | 550,548 | 0.01 | 2.06 | 1.08 | 3.34 | 0.01 | 28.2 | 11.8 | 0.74 GiB |
| overlap 0, batch 16 | 38 | 460,436 | 0.01 | 1.82 | 0.91 | 3.08 | 0.01 | 25.3 | 12.6 | 0.73 GiB |
| overlap 0, batch 64 (one batch per code) | 38 | 460,436 | 0.01 | 1.79 | 1.03 | 2.36 | 0.01 | **22.6** | 11.3 | 0.82 GiB |
| overlap 0, batch 64, window 24,576 | 20 | 460,436 | 0.01 | 1.71 | 1.04 | 3.13 | 0.01 | 25.6 | 12.8 | 0.90 GiB |

(stage columns in CPU seconds; `/usr/bin/time` user + system is 0.62–0.75 s
above the stage sum, the interpreter, torch import and checkpoint load;
engels-0084 reconciled both denominators, whole-process 36.5 / 25.4 CPU-s
per genome Mb for the first and fourth rows.)

*Corrections after review (commit `92ddafd`, `source_dirty: false`;
`chromosome-grid/`).* stalin-0085 showed that the `c6e5fc7` path started
the encoder at each window's local zero, so with step 8,192 (≡ 8 mod the
pooling stride 12) an interior base was pooled with different neighbours
in the two windows that saw it (emission differences up to 0.06 on an
untrained encoder), violating the fixed-grid condition of proposal
section 3.5. `predict_sequence` now encodes each tile from the stride
multiple at or before its start and crops the emissions back to the tile,
so the grid is the oriented chromosome's whatever the window and overlap
are; a real-encoder test compares shared interior positions on both
strands under a misaligned step and fails on the old code. engels-0084
showed the same path buffered every window's emissions of a strand before
its first decode (88 bytes per oriented base at float64, i.e. 13.2 GB per
strand of a 100 Mb chromosome, over the 8 GB cap before anything else);
windows are now flushed through the scan and core claiming in groups of
`decode_batch`, and a spy test asserts no more than that many emission
tensors are alive at any decode. Re-measured on the same core with the
same checkpoint: the default row is **32.4** CPU-s per genome Mb
(preprocess 2.55, encoder 1.25, decode 3.65; 1,832 chains, 6,544 rows;
RSS 0.74 GiB) and the overlap-0 / batch-64 row **21.3** (1.69 / 0.88 /
2.32; RSS 0.81 GiB) — inside the run-to-run noise of the table above (the
encoder input grows by at most 11 bases per window), so the conclusion
below stands on either set. The RSS of a 230 kb chromosome cannot show
the chromosome-independence of memory; that is what the buffering test
asserts, and a metazoan chromosome row will show it in `/usr/bin/time`.

*Revision step 1, the vectorised featurizer (commit `48f8a2b`,
`source_dirty: false`, source SHA-256 `e81253eb…`; `chromosome-vec/`).*
`encode_sequence` is now torch ops over the byte codes of the window
(`base_codes`: a 256-entry lookup for the base index and an ASCII
lower-case test for soft masking; prefix sums for the GC window; the GC
fraction formed in float64 and rounded to float32 once, as before), and
the availability mask is read through the buffer protocol instead of
`torch.tensor(list)`. The per-base Python implementation is kept as
`encode_sequence_reference`; a test asserts bit-identical tensors over
44 random windows drawn from both cases of ACGT, N, IUPAC codes and
non-letters, with and without availability masks. The pooled decoder's
per-window dinucleotide lookup (`_dinuc_index_tensors`) was vectorised the
same way against its list reference. **The GFF3 written by the default and
the overlap-0 / batch-64 runs is byte-identical to the `92ddafd` output**
(`cmp` on the two pairs), so the featurizer change is a cost change only.
Same core, thread, checkpoint and `/usr/bin/time -v` protocol as above:

| run | windows | oriented bases | io | preprocess | encoder | decode | output | **CPU-s per genome Mb** | CPU-s per oriented Mb | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| overlap 4,096, batch 16 (defaults) | 56 | 681,620 | 0.01 | 0.03 | 1.23 | 3.57 | 0.01 | **21.1** | 7.1 | 0.73 GiB |
| overlap 4,096, batch 64 | 56 | 681,620 | 0.01 | 0.04 | 1.27 | 2.88 | 0.01 | 18.2 | 6.2 | 1.07 GiB |
| overlap 2,048, batch 16 | 46 | 550,548 | 0.01 | 0.03 | 0.99 | 3.15 | 0.01 | 18.2 | 7.6 | 0.74 GiB |
| overlap 2,048, batch 64 | 46 | 550,548 | 0.01 | 0.03 | 1.01 | 2.50 | 0.01 | 15.5 | 6.5 | 0.93 GiB |
| overlap 0, batch 16 | 38 | 460,436 | 0.01 | 0.02 | 0.82 | 2.92 | 0.01 | 16.4 | 8.2 | 0.73 GiB |
| overlap 0, batch 64 (one batch per code) | 38 | 460,436 | 0.01 | 0.02 | 0.84 | 2.23 | 0.01 | **13.5** | 6.8 | 0.81 GiB |

The featurizer went from 2.55 to 0.03 CPU seconds on the default run
(11.1 to 0.1 CPU-s per genome Mb; 0.37 ms per 12,288-base window on one
thread, 140× the reference), the other stages are unchanged within noise.
**CPU numerators, made explicit after engels-0086:** the bold column is
the *stage sum of `process_time()`* (user + system inside the timed
stages, 0.6 s of interpreter, torch import and checkpoint load excluded);
`docs/cost-baseline.md` section 3.2 rows are *whole-process user only*;
whole-process user + system is a third figure. For this table they are
21.1 / 21.0 / 23.6 (default) and 13.5 / 14.2 / 16.2 (overlap 0 / batch
64) CPU-s per genome Mb; the untimed 0.6 s is 2.6 CPU-s per Mb on a
230 kb chromosome and negligible on a metazoan one. So the overlap-0 /
batch-64 row was the first chromosome row inside 15 on the stage metric
and on the baseline's user-only metric, not on user + system; the default
row (21.1) and the overlap-2,048 rows (15.5–18.2) missed on every
numerator. The paragraph after step 2a is the pre-revision reading and is
kept as the record of the miss.

*Revision step 2a, the decode operands (commit `6cae82f`,
`source_dirty: false`, source SHA-256 `3b37e46e…`; `chromosome-ops/`).*
Before touching the decode stage its 3.7 CPU seconds (default run) were
split with `process_time()` wrappers: the per-step scan loop 2.57, the
operand construction `_operands` 0.88, and the Python-side work that step
2 was going to remove only 0.26 in total (chain reconstruction 0.08,
traceback 0.07, symbol lookup 0.03, bias 0.02, padding and batch
assembly 0.07). A numpy twin of the scan loop, bit-identical in score and
back-pointers, ran 2.0× faster at batch 1 but 1.0× at batch 16 and 0.9×
at batch 56, so the loop is bound by the `(B, K, K)` transition
arithmetic (K = 24 states), not by dispatch; the numpy twin is not kept.
What was cheap and exact: `_operands` no longer materialises the
`(B, K, L, R)` intron and window-sum stacks (2 × 113 MB at batch 16,
2 × 450 MB at 64) but keeps them per phase, `(B, 3, L, R)`, and the scan
expands each step's slice to states with one `index_select`; the two
per-step `int8` casts are gone (the back-pointer stores convert on
`copy_`); and the batch assembly maps a window's symbols through a
256-entry byte table (`symbol_index_tensor`, tested equal to the list
reference over every ASCII byte). **All six GFF3 files are byte-identical
to the `chromosome-vec/` outputs** (`cmp`), so this is again a cost
change only. Same core, thread, checkpoint and `/usr/bin/time -v`
protocol; the three CPU numerators are given side by side:

| run | windows | preprocess | encoder | decode | **stages, user + system / genome Mb** | whole-process user / Mb | whole-process user + system / Mb | per oriented Mb (stages) | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| overlap 4,096, batch 16 (defaults) | 56 | 0.04 | 1.22 | 2.99 | **18.5** (was 21.1) | 19.4 | 21.2 | 6.3 | 0.54 GiB |
| overlap 4,096, batch 64 | 56 | 0.04 | 1.30 | 2.38 | 16.2 | 16.7 | 18.8 | 5.5 | 0.72 GiB |
| overlap 2,048, batch 16 | 46 | 0.03 | 0.99 | 2.65 | 16.0 | 17.4 | 18.7 | 6.7 | 0.54 GiB |
| overlap 2,048, batch 64 | 46 | 0.03 | 1.02 | 2.06 | **13.6** | 14.4 | 16.2 | 5.7 | 0.65 GiB |
| overlap 0, batch 16 | 38 | 0.02 | 0.83 | 2.49 | 14.6 | 16.1 | 17.4 | 7.3 | 0.54 GiB |
| overlap 0, batch 64 (one batch per code) | 38 | 0.02 | 0.83 | 1.86 | **11.9** (was 13.5) | 13.1 | 14.5 | 5.9 | 0.59 GiB |

(io and output 0.01 s each; the untimed whole-process remainder is
0.60–0.64 s in every run.) The decode stage fell 14.6–17.6% across the
six layouts (stalin-0088 corrected the range) and peak RSS by 0.18–0.35
GiB. On the stage metric the default row is now
**18.5** (1.2× over 15), the overlap-2,048 / batch-64 row **13.6** and both
overlap-0 rows are inside; on whole-process user only, 2,048 / 64 (14.4)
and 0 / 64 (13.1) are inside; on whole-process user + system only 0 / 64
(14.5) is, and on a 230 kb chromosome that numerator carries 2.6 CPU-s/Mb
of fixed start-up. The budget is therefore met on every numerator only at
the no-overlap cost floor, and at overlap 2,048 (a 1,024-base
containment guarantee, too short for most metazoan genes) on two of
three. No B allowance is claimed. What remains in the decode stage is the
scan loop itself, 2.0–2.6 s here: per step it costs about 25 µs fixed
plus 1.7 µs per window in the batch, the latter being the dense
`(B, 24, 24)` transition max although each state has one predecessor
except the phase-0 empty prefix (17) and `U` (3). Step 2b is that
sparse-predecessor scan (gather over at most 17 candidates instead of the
dense 24 × 24; exact, since the maximum runs over the same live
candidates in the same index order); its estimated payoff is 10% at batch
16 and 30–45% at batch 28–56, which together with the larger batch would
put the default overlap near 15. Step 3, the overlap factor, still waits
on boundary support.

*Revision step 2b, the sparse-predecessor scan (commit `d9ea0f8`, source
SHA-256 `5abe2f5507e7b72e33ae0779043cb598239ade0614f66cfdaa3b1da1218faf41`; records in
`smoke-local-20260920/chromosome-sparse/`).* Every coding state has one
predecessor except `U` (itself and the two stop-completing two-base
prefixes under code 1) and the empty prefix `E("")` (all sixteen two-base
prefixes plus the initiator-completing `S` prefixes), so the scan now
gathers the merged layer at `K − 2 + 2P` candidate slots (`P` = 17 for
code 1: the two multi-predecessor states padded to a common width under a
`FLOOR` prior), adds the symbol's prior row, `U`'s column extras on its
own candidates, takes one `max` over the two `(B, P)` candidate rows and
reassembles the `(B, K)` layer with one `cat` — 56 elements per window
per step (`24 − 2 + 2 × 17`) instead of the dense 576, two more dispatches. Past a window's
end the pad symbol row keeps `U → U` at 0 and everything else at `FLOOR`,
so the score is still read at boundary `L`; back-pointers of unreachable
states and of padded boundaries are unspecified (documented), and the
dense scan is kept as `viterbi_batch_reference` with a test that the two
agree bit for bit on scores, tail back-pointers and every traceback over
codes 1 / 6 / alternative initiators, float32 / float64, `m` = 1–20,
`-inf` masks and padded batches. The operands are also built step-major
(`(L, B, …)`, contiguous per step): the earlier `(B, K, L)` stacks sliced
along `L` were strided, and `index_select` cloned the slice — or its
strided index — every step. Transposing the `K`-wide stacks afterwards
cost more than that saved (`_operands` 0.17 → 0.46 s at 56 windows), so
only the `(B, L)` channels are transposed before the `K`-wide expansions
(`_operands` 0.30 s; the remainder is the element-wise phase gather in
the new layout). In isolation the scan over 56 windows × 12,288 steps
fell 1.35 → 0.97 s and over 16 windows 0.60 → 0.50 s. **All six GFF3
files are again byte-identical to the `chromosome-vec/` outputs** (`cmp`);
same core, thread, checkpoint and `/usr/bin/time -v` protocol:

| run | windows | preprocess | encoder | decode | **stages, user + system / genome Mb** | whole-process user / Mb | whole-process user + system / Mb | per oriented Mb (stages) | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| overlap 4,096, batch 16 (defaults) | 56 | 0.04 | 1.25 | 2.64 | **17.1** (was 18.5) | 18.1 | 19.7 | 5.8 | 0.60 GiB |
| overlap 4,096, batch 64 | 56 | 0.03 | 1.27 | 1.95 | **14.2** (was 16.2) | 14.4 | 16.8 | 4.8 | 0.81 GiB |
| overlap 2,048, batch 16 | 46 | 0.03 | 1.01 | 2.47 | 15.3 | 16.5 | 18.0 | 6.4 | 0.60 GiB |
| overlap 2,048, batch 64 | 46 | 0.03 | 0.98 | 1.74 | **12.0** (was 13.6) | 12.8 | 14.7 | 5.0 | 0.72 GiB |
| overlap 0, batch 16 | 38 | 0.02 | 0.84 | 2.39 | 14.2 | 15.6 | 16.8 | 7.1 | 0.60 GiB |
| overlap 0, batch 64 (one batch per code) | 38 | 0.02 | 0.84 | 1.57 | **10.7** (was 11.9) | 11.9 | 13.3 | 5.3 | 0.65 GiB |

(io and output 0.01 s each.) Against the same-named `chromosome-ops/`
records the decode stage fell 4.3–11.7% at batch 16 (11.7% on the
default overlap 4,096 row, 6.8% at 2,048, 4.3% at 0) and 15.5–18.0% at
batch 64. The default row therefore slightly exceeds the 10% estimated
above, while the batch-64 rows fall well short of the 30–45% estimate,
because that estimate counted the `(B, K, K)` block only and the
per-step fixed cost (about 25 dispatches, ~35 µs) is untouched; peak RSS
rose 0.05–0.09 GiB (the step-major intron tables). On the stage metric the
default row is **17.1** (1.14× over 15); **overlap 4,096 / batch 64 is the
first default-overlap row inside the budget on the stage sum (14.2) and on
whole-process user (14.4)**, though not on user + system (16.8, of which
2.6 CPU-s/Mb is fixed start-up on a 230 kb chromosome); overlap 2,048 /
batch 64 is inside on all three numerators (12.0 / 12.8 / 14.7), as is
overlap 0 / batch 64 (10.7 / 11.9 / 13.3). The default batch of 16 is a
memory choice (0.60 vs 0.81 GiB here; both far inside 8 GB), so the
overlapping configuration that carries the 2,048-base containment
guarantee now lands inside 15 on the baseline's numerator at batch 64 —
but only on the yeast chromosome, with a 20-step checkpoint and without
boundary support, and the S. pombe normalization row of section 6.1 is
still unmeasured. No B allowance is claimed on this basis; what remains
of the decode stage is the per-step fixed cost of the Python loop, which
is the checkpoint-and-replay / fused-kernel territory of proposal 3.3,
not a further exact reshuffle.

*Boundary support, first half: the sequence-edge partials of proposal 3.1
in the tensor decoder* (`fast_viterbi.viterbi_batch_edges`, `viterbi` /
`viterbi_windows(edges=)`; `EdgePrior` widened to the four partial
families, `pooled.edge_prior`). Under an `EdgePrior` every `E(q)` may be
entered at boundary 0 as `E0(q)` (coding entry plus the normalized
phase/prefix prior, first base CDS: the donor row at `t = 0` is floored)
and every tail `T(E(q), r)` as the residual intron `J` (intron entry plus
prefix prior and `log π`; the acceptor row at `t = 0` is floored so it
consumes one intronic base first); at a window's own end `n` an `S`/`E`
state exits with the coding exit, a tail entered by a real donor with the
intron exit and no `(1 − q)`, and a donor still parked at `s > n − m` as
the censored `I(c, n − s)` with its intronic bases summed. `J` is never
terminal, so it cannot share the tail layer with `T`: the first version
folded it in and told the two apart only at `n` (a tail is `J` when
`entered[1..n, c, r]` is all false), and engels-0088 / stalin-0090 showed
that this loses the window's only valid terminal path whenever the `J`
entry outscores the donor entry on the same row (the max keeps `J`, which
is then floored at `n`; a 4-base fixture with closed-form optimum
`6 − log 3 − 4 log 2` decoded to a 4-base CDS instead, in 144 / 144
sweep cases, and the score fell as the `J` weight rose). The scan now
carries `J` in its own `(B, K, R)` layer: it advances with the same
intronic emissions, competes with `T` at every acceptor (`exit_r` records
`R + r` so the traceback runs it back to boundary 0), never receives donor
entries and is never offered at `n`, which is therefore a one-time
reduction over the donor-entered layer alone. The extra layer costs one
add and one max per step in edge mode only: 0.63 vs 0.54 s over 16 ×
12,288 on the same micro-benchmark (the earlier 0.51 vs 0.52 s was the
defective version; edge mode is meant for the two real ends of a strand,
not every tile, so the chromosome numbers above are unaffected). The
traceback takes the chosen final state and relabels a `J` exit and a
coding state at boundary 0 as `J` / `E0`, so `_chains` marks `partial_5` /
`partial_3` as for the reference. `tests/test_a_fast_viterbi.py::EdgeParity`
holds it to `DelayedEntryDecoder(code, duration, edges).viterbi` — 240
random lattices with three priors (over 150 finite; 5′-partial, 3′-partial
and doubly partial chains each present), fixtures in which each of the four
scalars moves the score alone, batched (1 / 7 / 64, padded lengths, an
empty window) equal to single, the learned tables with the pooled prior,
the engels-0088 fixture across six `J` weights (score and chain fixed at
the closed form), and stalin-0090's boundary sweep (codes 1 / 6, both
dtypes, `m` 2 / 4 / 20, `R` 1 / 3, `n = m − 1 .. m + 2`, both `J`
weights, batched 1 / 7 / 64 equal to single); 186 tests pass under 3.11 +
torch. This decoder is what a true
chromosome end needs (`chromosome.py` still decodes the free grammar, so a
gene cut by a sequence end is lost today) and what the section-3.6
edge-partial training numerator will use once `fast_loss` has the same
entry/exit terms (the partition needs a separate `J` layer, since a sum
cannot recover the `T`-only part at `n`; not done).

*Revision step 3 needs a decision on what "boundary support at a seam"
means.* Proposal 3.1 and 3.6 are explicit that interior chunk boundaries
get no partial entry or exit — a window seam is not a sequence edge — so
the edge decoder above must not simply be switched on at every tile edge
to drop the overlap; that would price a seam-crossing gene as two
partials and change the model. The exact route is proposal 3.3's
checkpointed seams: carry the scan state (`alpha`, `tau`, the `m` pending
donor layers) from one tile into the next of the same strand, with the
batch dimension over independent segments (strand × contiguous segment,
each row scanning its segment tile by tile), overlap only at the few
segment seams, and either the back-pointers of a whole segment held (about
100 bytes per base as stored now — too much for 8 GB above ~50 Mb per
batch — so they would need packing) or the 3.3 replay (a second scan,
which on the numbers above costs more than the overlap it removes). The
next step is therefore to measure the carried-state scan at overlap 0
against the overlap-4,096 rows, not to relax the seam contract; the
edge-enabled decoder is used at the two real sequence ends of each
strand either way.

*Revision step 3, first measurement: the carried-state scan* (commit
`a4d393d`, source SHA-256
`10c4f9f404f601381d56491d12649e1bb36552c2218a851ecd3275a81dc70710`;
records in `smoke-local-20260920/chromosome-seams/`).
`fast_viterbi.scan_segments` scans a batch of *segments* `tile` bases at a
time and carries the state across each seam as a `Carry` (`alpha`, `tau`,
the residual-intron layer, the `<= m − 1` pending donors with the matching
seam-context emissions so the next tile re-reads their mandatory windows
from the same operands, the edge-mode running score / final state, the
tails a donor has entered so far, and the rows already ended); the next
tile's loop starts at `t0 = len(pending)`, the boundary-0 entries and
`t = 0` floors of the edge mode apply to the first tile only, and a row's
terminal choice is made in the tile where its own length ends. Scores are
bit-identical to the unbroken scan and so are `exit_r` / `entered` inside
every row and every traceback (`tests/test_a_fast_viterbi.py::SeamCarry`:
200 random batches of 1–6 segments at tiles 1 / 2 / 3 / 5 / 7 / 11 / 64,
`m` 1–4, `R` 1–3, codes 1 / 6 / alternative initiators, `-inf` masks, with
and without an `EdgePrior`, over 400 finite cases of which over 100 have a
segment crossing a seam; 189 tests pass under 3.11 + torch). Only one
tile of operands is alive at a time; the back-pointers of the whole
segment are still held (`K (2 + R)` = 120 bytes per base at `R = 3`), so
this is the "held, not packed" option above. `chromosome.predict_sequence(segments=)`
(`measure --segments`) cuts each strand into about that many segments
overlapping by `--overlap`, encodes every segment in non-overlapping
12,288-base tiles on the anchored pooling grid, and decodes all segments
of both strands in one batch (`B` = 2 × segments), so the overlap factor is
`1 + (segments − 1) × overlap / n` per strand instead of `window / (window
− overlap)` and the containment guarantee holds at segment seams only; a
gene crossing a tile seam inside a segment is decoded whole.
`tests/test_a_chromosome.py` adds the seam-crossing cases (one segment in
50-base tiles with both genes across a tile seam, two segments, tiles of
7). Same core, thread, checkpoint and `/usr/bin/time -v` protocol as the
rows above; the two window-mode rows re-run at this commit are
byte-identical to the `chromosome-sparse/` GFF3 outputs:

| run | rows in the batch | encoder tiles | oriented bases | encoder | decode | **stages, user + system / genome Mb** | whole-process user / Mb | whole-process user + system / Mb | chains | chains only in this row / only in the exact decode | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| windows, overlap 0, batch 64 (cost floor, re-run) | 38 | 38 | 460,436 | 0.83 | 1.56 | **10.5** | 11.7 | 13.2 | 1,853 | +65 / −44 | 0.65 GiB |
| windows, overlap 4,096, batch 64 (re-run) | 56 | 56 | 681,620 | 1.22 | 1.98 | **14.1** | 14.6 | 16.7 | 1,832 | +5 / −5 | 0.81 GiB |
| 1 segment per strand (exact strand decode) | 2 | 38 | 460,436 | 0.88 | 8.56 | 41.3 | 42.6 | 43.9 | 1,832 | — | 0.52 GiB |
| 8 segments per strand, seam overlap 4,096 | 16 | 48 | 517,780 | 0.96 | 1.93 | **12.7** | 13.6 | 15.3 | 1,831 | +6 / −7 | 0.70 GiB |
| 19 segments per strand, seam overlap 4,096 | 38 | 76 | 607,892 | 1.09 | 1.57 | **11.8** | 12.3 | 14.3 | 1,832 | +11 / −11 | 1.02 GiB |
| 32 segments per strand, seam overlap 4,096 | 64 | 64 | 714,388 | 1.28 | 1.57 | **12.6** | 12.9 | 15.2 | 1,832 | +10 / −10 | 1.37 GiB |

(preprocess 0.02–0.04 s, io and output ≤ 0.03 s.) Three things follow.
First, the scan's per-step fixed cost is the whole story on a 230 kb
chromosome: the exact two-row decode costs 8.6 s (230 k steps at ~37 µs
each, 5.5× the batch-64 window decode), so the batch dimension must be
filled by cutting each strand into segments, and every segment seam
brings back the overlap. At 19 segments the batch equals the overlap-0
window row's, the decode stage is the same 1.57 s, and the encoder
carries the 1.32× oversampling — **11.8 CPU-s/Mb on the stage sum, 12.3
whole-process user, 14.3 user + system, all inside 15, with the 4,096-base
seam overlap retaining the 2,048-base chain-containment guarantee (at
segment seams; `chromosome.tiles` documents the bound, and a 4,096-base
chain starting one base before a core boundary is the counterexample to
any stronger claim) and every interior tile seam exact**; 8 and 32 segments bracket it (12.7 / 12.6). Second, this
mode's advantage is proportional to chromosome length: with the same 38
rows the oversampling is `1 + 18 × 4,096 / n`, 1.32 on yeast chr I and
1.004 on a 20 Mb chromosome, where the window mode's 1.5× is fixed; the
metazoan row of the definition of done is where the gap shows, and there
the held back-pointers (120 B × 40 M oriented bases = 4.8 GB) plus the
segment emissions (44 B per base in float32) are what forces the packing
or replay before the 8 GB cap — so the next increment is the packed
back-pointer store (`entered` as bits and `exit_r` / `prev2` as one byte
each, ~35 B per base), not more segments. Third, the segment rows still
differ from the exact decode by 6–11 chains of 1,832 (the window row at
overlap 4,096 by 5), and not only at the 7–31 segment seams per strand:
the encoder is run per 12,288-base tile with no context past the tile,
so emissions within its receptive field of a tile seam depend on the
tile grid (the tiles of a segment start at the segment's start, which
moves with the segment count). An encoder context margin (encode the
tile plus the receptive field on each side on the anchored grid and crop,
as the window mode implicitly does through its overlap) is the cheap fix
and is measured below (commit `8440560`); the decode itself is
seam-exact, as the tests show. Peak
RSS grows with the batch (the per-tile operands at `B` × 12,288, as in
the batch-64 window rows) and with the held back-pointers, 0.52–1.37 GiB
here, all far inside 8 GB. No B allowance is claimed on this basis: it is
one yeast chromosome, a 20-step checkpoint, free-grammar decoding at the
two real strand ends, and the S. pombe row is still unmeasured.

**Revision step 3, packed back-pointers** (commit `8a7b4a8`; records
`chromosome-seams-packed/`). `fast_viterbi.scan_segments` now returns a
`PackedBackPointers` store instead of the three dense tensors: `prev` as
the two argmax slots of `U` and `E("")` that the scan already computes
(2 bytes per boundary; every other state has one predecessor, so the
dense `(K,)` column is a constant expanded on traceback), `exit_r` as
nibbles (12 bytes: the values `-1 .. 2R - 1` fit four bits while `2R <
16`) and `entered` as `K R = 72` bits (9 bytes) — **23 bytes per base
instead of 120** (logical payload; since `87bb3a1` every tensor the store
retains owns exactly its storage, so `storage_nbytes()`, the distinct
allocations deduplicated by data pointer, equals `nbytes()` — before that
the superseded seam columns and the seam context stayed alive behind
views, stalin-0093), packed per tile as the scan produces them and expanded
one row at a time for the traceback, so one segment's dense pointers are
alive at a time rather than the batch's. `SeamCarry` now also checks the
store's byte count and that each row's expansion equals the unbroken
scan's `prev` (for `t ≥ 1`), `exit_r` and `entered` bit for bit; 189
tests pass. Re-measured under the same protocol:

| run | encoder | decode | **stages / genome Mb** | whole-process user / Mb | user + system / Mb | chains | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 segment per strand (exact strand decode), packed | 0.86 | 8.31 | 40.1 | 41.7 | 42.8 | 1,832 | **0.45 GiB** (was 0.52) |
| 19 segments per strand, seam overlap 4,096, packed | 1.10 | 1.60 | **11.9** | 12.4 | 14.6 | 1,832 | 1.02 GiB (was 1.02) |

Same chains as before (the GFF3 outputs are byte-identical to the
`chromosome-seams/` runs of the same configuration). The packing costs
nothing measurable in the decode stage (1.57 → 1.60 s is inside the
run-to-run noise of these 2–3 s processes), and it shows in RSS only where
the pointers were the largest live tensor: the two-row exact decode drops
by 70 MB (460 k oriented bases × 97 bytes saved is 45 MB, the rest the
transient dense expansion of both rows at once that the per-row expansion
avoids), while the 19-segment row's 1.02 GiB is the per-tile operands at
`B = 38` (`(12,288, 38, 24)` float64 tensors, ~90 MB each, several alive
per tile) and does not move — on a 230 kb chromosome the pointers were
never the peak. The metazoan projection this paragraph first carried
("about 3.5 GB") was wrong and is withdrawn: `measure` decodes in
float64, not float32, and at `8a7b4a8` the chromosome mode retained every
segment's emissions *and* `viterbi_segments` padded them into a second
`(B, 11, L)` tensor, so for 20 Mb / 19 segments per strand (38 segments
of 1,056,512 oriented bases, 40.1 M in all) the two 88 B/base copies
alone were 7.07 GB and the live tensors just before the last tile at
least 8.70 GB (engels-0090's source-derived lower bound, reproduced by
stalin-0093) — over the cap before any runtime overhead.

**Revision step 3, streamed emissions** (commit `87bb3a1`; records
`chromosome-stream/`). `scan_segments` now also accepts a callable
`chunk(start, end)` in place of the emission tensor and asks it for one
tile of every row as the scan reaches it; `viterbi_segments` takes
per-segment emitters `f(start, end)` the same way, and the chromosome
mode encodes each tile (encoder on the pooling-grid-anchored slice, the
dinucleotide bias with its two-base / one-base context inside the
segment) when the scan asks for it, so no per-segment emission list and
no padded copy exist. The store also owns its storage now (above).
`SeamCarry` checks that streamed emissions give the scores and chains of
the whole tensors and that each emitter is asked for exactly the tiles of
its own segment once, and that `storage_nbytes()` equals `nbytes()`;
190 tests pass. Re-measured under the same protocol:

| run | encoder | decode | **stages / genome Mb** | whole-process user / Mb | user + system / Mb | chains | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 segment per strand (exact strand decode), streamed | 0.89 | 8.69 | 41.9 | 43.4 | 44.5 | 1,832 | **0.35 GiB** (was 0.45) |
| 19 segments per strand, seam overlap 4,096, streamed | 1.17 | 1.62 | **12.3** | 12.6 | 14.9 | 1,832 | **0.96 GiB** (was 1.02) |

GFF3 outputs byte-identical to the `chromosome-seams/` and
`chromosome-seams-packed/` runs of the same configuration (the encoder
tiles are the same slices as before). The times are inside the noise of
these 3–10 s runs (decode 1.60 → 1.62 s at 19 segments; the exact decode
8.31 → 8.69 s, its `B = 2` per-step overhead unchanged); the encoder
stage is now clocked inside the scan and taken back out of the decode
stage, so the stage sum is still additive. RSS moves by what the retained
emissions and the padded copy weighed on a 230 kb chromosome (2 × 88 B ×
0.46–0.61 M bases = 80–107 MB): the 19-segment row's 0.96 GiB is the
per-tile operands at `B = 38` plus the ~0.3 GiB runtime floor, and is
independent of chromosome length. Corrected metazoan projection, 20 Mb /
19 segments per strand / seam overlap 4,096 / tile 12,288, `B = 38`,
`K = 24`, `R = 3`, `m = 20`, from the same helpers and operand formula as
engels-0090: packed pointers 38 × 1,056,513 × 23 B = **0.92 GB**; one
tile's operands (`8 B (L (2 + 3K + P) + (2L − m + 1) 3R)` at `L = 12,051`
with the 19-base seam context, `P = 17`) **0.40 GB**; one tile of
emissions and symbols for every row, float64, 38 × 12,288 × 96 B =
**0.04 GB**; one row's dense expansion during its traceback 1,056,513 ×
120 B = **0.13 GB**; carry, seam context and packed seam columns under
1 MB — **about 1.5 GB of live tensors** plus the runtime floor, against
8 GB. This is a source-derived projection like the withdrawn one, not a
measurement, and the 0.92 GB grows with the chromosome while the rest
does not; the metazoan row itself and the S. pombe row are still to
come.

**Revision step 3, encoder context margin** (commit `8440560`; records
`chromosome-margin/`). In the segment mode each tile's encoder input is
now the tile plus `DEPENDENCY_RADIUS` = 491 bases of the oriented
chromosome on each side (only what exists at the true ends), on the
same origin-anchored pooling grid, cropped back to the tile before the
dinucleotide bias (`predict_sequence(margin=)`, `measure --margin`,
default 491; `--margin 0` is the previous bare tile). An emission then
depends only on the bases within the encoder's dependency radius of its
own position, whatever the tile grid: `SegmentMargin` checks, with the
real `CandidateA` on a 6,100-base sequence with lowercase and `N`
bases, two segments per strand and 1,500-base tiles, that every tile's
emissions equal the whole-strand encoder output on that slice to
2.4 × 10⁻⁷ (float32 summation order; 0 on most tiles), that the bare
tiles differ by 0.03–0.13 wherever a tile does not start at its
segment's first base, and that `encoded_bases` (now in the record) is
exactly the sum of the clipped `[a − 491, b + 491)` ranges on the stride
grid; 191 tests pass under 3.11 + torch. Re-measured under the same
protocol (yeast chr I, one thread pinned to one core, `best.pt` from
the smoke fit, float64 decode):

| run | encoder | decode | **stages / genome Mb** | whole-process user / Mb | user + system / Mb | chains | peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 segment per strand (exact strand decode), margin 491 | 0.90 | 8.48 | 41.0 | 42.7 | 43.7 | 1,832 | 0.36 GiB |
| 19 segments per strand, seam overlap 4,096, margin 491 | 1.32 | 1.66 | **13.2** | 13.3 | 15.9 | 1,832 | 0.96 GiB |

The encoder stage grows by the encoded bases: 1.077× on the exact
decode (38 tiles of 12,288 or the 8,986-base remainder) and 1.120× at
19 segments (76 tiles, of which 38 are the 3,710-base segment tails,
where 982 bases of margin weigh more), 1.17 → 1.32 s; on full 12,288-base
tiles the margin is (12,288 + 982) / 12,288 = 1.08× of the encoder
stage, about +0.4 CPU-s per genome Mb on the ~4.9 the encoder costs at
scale, so the 20 Mb projection barely moves. Decode and RSS are
unchanged inside noise (the exact decode 8.69 → 8.48 s; system time
0.53 → 0.58 s is the page-fault cost of the 1 GB working set, which is
what separates the user + system column from the other two on this
230 kb chromosome). What the margin buys is in the chains: the exact
strand decode and the 19-segment decode now differ by **4 chains of
1,832** (were 11 at `87bb3a1`, 6–11 across the earlier segment rows;
the exact decode itself changed by 3 chains against its bare-tile
predecessor, its own 12,288-base tile seams having been the same
defect), and all four are plus-strand genes of 2,886–4,614 bases that
start 844–1,166 bases before a segment core seam and, being longer than
the 2,048-base containment guarantee of the 4,096-base seam overlap,
run off the end of the segment that owns their start (the segment
version ends 10–1,612 bases short of the exact one, at the segment's
last bases). That is the seam contract of proposal 3.1 working as
specified — a chain longer than `overlap − overlap // 2` may be cut at a
segment seam — and no longer a tile-grid artefact: interior tile seams
are now exact in both the encoder and the decoder, and the only
remaining seam effect is the one the overlap parameter controls. It
matters for the metazoan row, where genes routinely exceed 2 kb: the
seam overlap must scale with the gene lengths the panel expects
(8,192 for a 4,096-base guarantee costs 1.008× oversampling at 20 Mb /
19 segments, against 1.004× now), or the segment count per strand
must drop as the chromosome grows — at 20 Mb the batch is filled by
the chromosome, not by the segments, so the choice costs little there.
Records: `chromosome-margin/` (JSON, stdout, `/usr/bin/time -v`, both
runs, commit `8440560`, clean tree); the GFF3 outputs are not
committed.

**S. pombe normalization row** (proposal 6.1, cost-baseline 5.2;
records `pombe-normalization/`, code `8440560`, clean tree). The same
protocol on the three *S. pombe* nuclear chromosomes (12,571,820 bases;
the mitochondrion excluded), one process per run pinned to one core of
the Core Ultra 9 285K, frozen smoke checkpoint, float64 decode, margin
491, and AUGUSTUS 3.5.0 (bioconda, marx-0026's command
`--species=schizosaccharomyces_pombe --gff3=on --UTR=off`) on the same
nuclear FASTA, same core, immediately after. `benchmark/leakage_check.py`
ran first (0 violations; *S. pombe* is cross-clade held out, nearest
train species at phylum rank); the species was not scored and nothing
was selected on it — this row is the runtime denominator only.

| run | **stages / genome Mb** | whole-process user / Mb | user + system / Mb | peak RSS | chains or genes |
|---|---:|---:|---:|---:|---:|
| A, 1 segment per strand (exact strand decode) | 40.35 | 39.82 | 40.49 | 1.49 GiB (chr I) | 121,937 |
| A, 19 segments per strand, seam overlap 4,096 | **9.22** | 8.23 | 9.37 | 1.31 GiB (chr I) | 121,937 |
| AUGUSTUS 3.5.0, one process, whole nuclear genome | — | **51.41** | 51.43 | 0.40 GiB | 4,452 |

Stages at 19 segments per genome Mb: preprocess 0.11, encoder 3.74,
decode 5.31, output 0.06, I/O 0.01; the three chromosomes agree within
2 % (9.17 / 9.21 / 9.38 on the stage sum). AUGUSTUS costs 51.4 user
CPU-s/Mb here against 164.9 on the cost-baseline runner (marx-0026), so
this machine is 3.21× faster on the same command and the absolute
15 CPU-s/Mb of cost-baseline 5.2 is not directly applicable; the
portable form is. **A at 19 segments is 1/5.5 of AUGUSTUS on
user + system (1/6.2 on user only, 1/5.6 on the stage sum) against the
1/11 target, i.e. it needs ≤ 4.67 CPU-s/Mb on this machine and costs
9.4: the CPU regime misses the normalization row by 2.0× (1.8× on user
only).** Machine-normalized to the baseline runner the row reads
30.0 CPU-s/Mb against 15. Memory is inside 8 GB on every run. The exact
strand decode is 1/1.3 of AUGUSTUS. The segment decodes differ from the
exact ones by 3 chains of 121,937 (seam contract; `chaindiff.out`).
Where the 2× is: encoder 3.74 + decode 5.31 = 9.05 of the 9.22, so both
halves must roughly halve. A probe on chr III with 38 and 76 segments
per strand (`A-probe/`) gives stage 8.72 and 8.96 CPU-s/Mb — decode
11.30 and 11.00 s against 13.26 at 19 segments, the encoder up 4 % and
13 % with the extra margin bases — at RSS 1.87 and 3.29 GiB against
1.14: the decode stage is at its operand floor (~4.5 CPU-s per genome
Mb, 2.1 per oriented Mb after the 1.08× seam oversampling), and the
segment count is exhausted as a lever. **No positive CPU allowance for
B follows; the task's failed-regime clause applies.** The revision
proposed, in order of expected return per change, before the metazoan
row: (1) decode in float32 — the operands (emissions, the R = 3 running
scores and the pooled duration tables) are float64 today for the
parity tests; the scan is memory-bound per step, so halving the operand
bytes is expected to take the decode stage toward 2.5–3 CPU-s per
genome Mb, checked against the float64 path with the existing 2e-6
relative tolerance of the seam tests; (2) the encoder at 3.74 CPU-s per
genome Mb already runs the proposal's counted 231 GFLOP per genome Mb
(section 6, both orientations, assumed 30 GFLOPS → 7.72 s/Mb) at an
effective 62 GFLOPS on this core, twice the proposal's assumed rate, so
a further halving there is a change to the counted work (kernel width,
channel count or the dilation schedule of the tail blocks), which must
be measured with a fitted checkpoint's accuracy column before it is
adopted, not before;
(3) if (1) lands and (2) is deferred, the row would sit near
6.5–7 CPU-s/Mb, 1/7.5 of AUGUSTUS — still short of 1/11 — so the
honest expectation is that A meets the portable CPU target only with
(2) or on the GPU regime (gagarin, lenin-0083), and the coordinator's
decision on the budget should be taken with that in view. Records:
`pombe-normalization/` (README with the table, JSON, stdout,
`/usr/bin/time -v`, run scripts, leakage check, AUGUSTUS install log,
GFF3 SHA-256s; GFF3 outputs of 6–15 MB not committed).

**Revision step 4, float32 decode** (`aee0134`; records
`pombe-normalization/A-f32/`, same protocol, checkpoint and core).
`measure --dtype float32` puts the emissions, motif bias, duration
tables and the scan in float32; in the segment mode the carried scores
are rebased at every tile seam (`scan_segments(rebase=)`, default for
dtypes narrower than float64: each unfinished row's best coding-layer
score is shifted to 0 and the shifts summed back into the returned
score in float64), so the running magnitude is a tile's (~10⁴, ulp
~0.001) rather than a segment's (~3 × 10⁵, ulp ~0.03); the shift
changes no comparison, and `SeamRebase` checks that rebased float64
tiling equals unrebased tiling on offset lattices and that a float32
rebased scan gives the float64 chains (173 of 173 finite cases).
Float64 stays the bit-for-bit unbroken scan the seam tests verify.
Measured on the three *S. pombe* chromosomes: 19 segments **8.27 /
7.94 / 8.43** CPU-s per genome Mb (stage / user / user + system;
float64 9.22 / 8.23 / 9.37), RSS ≤ 0.97 GiB (was 1.31); exact strand
decode 39.72 / 39.22 / 39.86 (was 40.35 / 39.82 / 40.49). The decode
stage goes 5.31 → 4.34 per genome Mb at 19 segments and 36.4 → 35.8 on
the exact decode: **a 10 % gain on the row, not the halving proposed**,
because the decode is not bandwidth-bound. A profile of one tile scan
(12,288 steps at B = 38, one core) costs 0.93 s in float64 and 0.73 in
float32: the per-step loop is 0.52 / 0.47 s — 42 / 38 µs per step for
~15 torch ops on (38, 24) tensors, i.e. dispatch, which dtype cannot
touch — and the per-tile operand build 0.41 / 0.26 s, which is where
the float32 saving is. The exact decode, at B = 2 and 5.6 M steps per
chromosome, is almost entirely per-step dispatch, hence unmoved.
The row now reads 1/6.1 of AUGUSTUS on user + system: **the miss is
1.8× (1.7× on user only), still failed**. Outputs: float32 differs from
float64 by 103 chains of 121,937 at 19 segments (93 internal
exon-boundary shifts in short 2–5-exon chains, 7 exon-count changes,
3 without a partner) and 112 at 1 segment, and the float32 exact and
segment decodes differ from each other by ~120 chains where the float64
pair differ by 3 (`chaindiff_f32.out`) — near-tie decisions of the flat
smoke checkpoint resolved differently at float32 precision. Until a
fitted checkpoint shows the flip count to be negligible, float64 stays
the row of record and float32 is a measured option. The revision list
in the paragraph above is amended accordingly: (1) is measured and
worth 10 %, not 25–30 %; the decode lever that remains is the per-step
dispatch — a compiled scan step (one fused kernel per base instead of
~15 dispatched ops) would take the per-step half of the 19-segment
decode (~2.2 of 4.3 CPU-s per genome Mb) toward a few tenths and leave
the operand build (~1.6 in float32), so the row would sit near
6–6.5 CPU-s/Mb, 1/8 of AUGUSTUS, still short of 1/11; (2), the encoder
work reduction with a fitted checkpoint's accuracy column, is therefore
necessary for the portable CPU target, not optional, and the
expectation stated in (3) stands. Compute: 0.17 CPU-h local for the six
float32 runs; cluster CPU-hours 0, GPU-hours 0.

**Metazoan development-chromosome row** (code `aee0134`; records
`elegans-chrV/`, README with the table, JSON, stdout, `/usr/bin/time -v`,
run scripts, leakage check, source MD5s, GFF3 SHA-256s, chain diff).
*C. elegans* chromosome V (`NC_003283.11`, 20,924,180 bases; a train
species, chr V declared the development chromosome in the run config,
runtime only, unscored), the pinned WBcel235 sources MD5-verified, the
leakage check run first (0 violations), the same frozen smoke checkpoint,
machine, one process pinned to one core, and AUGUSTUS 3.5.0
(`--species=caenorhabditis`, marx-0026's command otherwise) on the same
chromosome on a second core. The seam overlap is **16,384** (an
8,192-base chain-containment guarantee, covering 97% of chr V's 4,965
admitted representatives; p99 13.7 kb) rather than the yeasts' 4,096,
as the encoder-margin paragraph required: 1.014× oversampling at 19
segments per strand. Results, per genome Mb (stage / user / user +
system): A at 19 segments, float64, **9.64 / 8.48 / 9.67**, RSS 2.04 GiB
(encoder 3.91, decode 5.50); float32 8.37 / 7.85 / 8.40, RSS 1.74 GiB;
the exact strand decode 42.27 / 41.98 / 42.30, RSS 4.87 GiB (its one-row
traceback expansion and packed store scale with the 20.9 Mb row; inside
8 GB, not the configuration of record); AUGUSTUS **52.19** user CPU-s/Mb,
0.70 GiB, 3,357 genes — within 2% of its 51.4 on *S. pombe* on this
machine, so the *S. pombe* machine factor of 3.21× applies (the chr V
comparison is same-machine, not same-core: A on core 2, AUGUSTUS on
core 4, concurrently for the float64 runs; the *S. pombe* one was
sequential on one core). **A is 1/5.4 of
AUGUSTUS on user + system (1/6.2 on user), machine-normalized 31.0 CPU-s/Mb
against 15: the metazoan row misses the CPU target by 2.0× (1.8× on
user), the same miss as the *S. pombe* row.** The per-genome-Mb cost is
chromosome-length independent within 5% between 12.6 and 20.9 Mb (the
residual is the longer overlap's oversampling, the margin on more tiles
and a higher system time, 25.0 s against 14.3 s, allocation churn the
float32 path halves). Outputs: the 19-segment float64 decode is
**byte-identical to the exact strand decode** (208,796 chains, 0 differ):
with the 8,192-base guarantee no chain of this checkpoint is cut at a
segment seam, where the yeasts at overlap 4,096 lost 3–4 per chromosome.
Float32 flips 194 / 199 chains of 208,796 (0.09%; *S. pombe* 0.08%), the
same near-tie behaviour, to be re-checked with a fitted checkpoint before
float32 becomes the row of record. `measured.tsv` gets four rows (AUGUSTUS
chr V same-machine, A exact, A 19 segments, A 19 segments float32;
`wall_s` is the process elapsed time from `/usr/bin/time -v` on every
row, the stage-wall sums stay in the README table, stalin-0098); the
definition of done's "at least S. pombe and one metazoan development
chromosome" is now met for the CPU regime, with the smoke checkpoint —
the fitted checkpoint re-measures both rows for the accuracy column and
the GPU regime is the other half. Compute: 0.35 CPU-h local for the three
A runs and 0.30 for AUGUSTUS; cluster CPU-hours 0, GPU-hours 0.

**The CPU regime missed the budget before revision step 1.** Per oriented megabase the pipeline
costs what items 4–5 measured (11.3–12.8 CPU-s), but a genome megabase is
two oriented megabases before overlap, so the end-to-end row is
**33.8 CPU-s/Mb with the default overlap and 22.6 at best (no overlap, one
decode batch) against the 15 CPU-s/Mb target** — 2.25× and 1.5× over. The
per-genome-Mb budget was always going to carry both strands (proposal
section 6: "count both strands"); items 4–5 were per oriented base and said
so, and this is the first row in the budget's own unit. I/O and output are
negligible (0.02 s together; the GFF3 is 6.5 k rows). Peak RSS 0.75–0.90
GiB, inside 8 GB. The 20-step checkpoint's predictions carry no information
(1,831 chains on a chromosome with ~100 genes; the GFF3 files are not
committed) — this row is the cost of the specified pipeline, not its
accuracy; the accuracy row needs the longer fit that waits on gagarin.
The overlap-0 rows are the cost floor for *this* code, not a proposal to run
without overlap: at overlap 0 a gene crossing a seam is lost (item 6's
guarantee holds only up to half the overlap). No duplicate CDS rows appear
across overlapping windows (checked on the default run).

**Proposed revision (task: "if A misses a target, report the failed regime
and propose the revision before any positive CPU allowance is assigned to
B").** The failed regime is CPU, one core, per genome Mb. In the best row
the three real stages are preprocess 7.8, encoder 4.5 and decode 10.3
CPU-s per genome Mb. The revision, in order of measured payoff:
(1) the featurizer is per-base Python (`encode_sequence`); vectorising it
over the window (numpy/torch ops on a byte array) should take its 7.8 down
to well under 1 — this is pure implementation, no model change;
(2) decode: the tensor scan's per-step cost is amortised by batch (batch 64
saves 3.1 CPU-s/Mb over batch 16 here; the chromosome had only 38 windows
per batch, a real chromosome fills batch 64+), and the per-window
dinucleotide lookup and Python-side batch assembly are still in the decode
stage; a tensor gather for the bias and a preallocated batch should take
another 1–2; (3) overlap: 2,048 instead of 4,096 (1.2× oversampling, genes
≤ 1,024 bases guaranteed) costs 12 % over no overlap instead of 33 %, and
the boundary-support increment (partial chains at window edges, section 3.6)
removes the need for a completeness guarantee from the overlap altogether.
With (1) and (2) the best row would sit near 13–15 CPU-s/Mb at overlap 0 and
16–18 at overlap 2,048; the budget is reachable on the encoder/decoder as
specified but not with margin, and **no positive CPU allowance for B
should be assumed until a chromosome row is measured inside 15** with the
revised featurizer and decode path. *Status:* step (1) is done, and step
(2) turned out on measurement to be mostly the operand construction rather
than the Python-side assembly (0.26 s of 3.7); its cheap half (2a) is
done, and its remaining half (2b, the sparse-predecessor scan) is done
and measured above (stage sum 10.7 at overlap 0 / batch 64, 12.0 at
overlap 2,048 / batch 64, 14.2 at overlap 4,096 / batch 64, 17.1 at the
batch-16 defaults; 11.9 / 12.8 / 14.4 / 18.1 on the whole-process user
numerator). Step (3) is measured in its first form: the carried-state
scan (seam contract of proposal 3.1 kept, interior tile seams exact)
gives 11.8 / 12.3 / 14.3 at 19 segments per strand on the three numerators
(11.9 / 12.4 / 14.6 with the packed store, 12.3 / 12.6 / 14.9 with the
emissions streamed, 13.2 / 13.3 / 15.9 with the encoder context margin,
interior tile seams then exact in both encoder and decoder; the 20 Mb
chromosome now projects to ~1.5 GB of live tensors, unmeasured); the B
allowance question stays open
until the overlapping configuration is also measured on the S. pombe
normalization row with a fitted checkpoint and true-edge support. The GPU regime is still unmeasured
(gagarin, lenin-0083); the section 6.1 multi-worker decoder accounting is
the other half of the same row.

**Compute recorded:** local CPU only — ~0.4 CPU-h for the two fast-kernel
smoke fits, ~0.02 CPU-h for the profile, 0.003 CPU-h for the chromosome I
`measure`, 0.005 CPU-h for the tensor-decoder `measure` runs and profiles,
~0.01 CPU-h for the pooled-decoder fit and its two `measure` runs, ~0.02
CPU-h for the ten chromosome-profile runs (five committed at `c6e5fc7`, five
at a dirty tree before it, superseded), ~0.01 CPU-h for the two `92ddafd`
re-measurements and the six `48f8a2b` runs, ~0.02 CPU-h for the decode
profiling, the numpy-twin experiment and the six `6cae82f` runs (plus one
mis-configured multi-threaded run, discarded), ~0.02 CPU-h for the
sparse-scan micro-benchmarks and twelve chromosome runs at `de0746a` /
`d9ea0f8` (the six `de0746a` runs, with transposed `K`-wide stacks,
superseded), ~0.05 CPU-h for the tests and the sixteen chromosome runs
of revision step 3 (`a4d393d`, `8a7b4a8`, `87bb3a1`, `8440560`), ~0.37
CPU-h for the S. pombe normalization row (six A runs 0.17 h, AUGUSTUS
0.18 h, two probe runs 0.01 h), 0.17 CPU-h for the six float32 runs of
revision step 4 (`aee0134`), 0.65 CPU-h for the *C. elegans* chr V row
(three A runs 0.35 h, AUGUSTUS 0.30 h); cluster CPU-hours 0, GPU-hours 0. The
only held-out species touched is *S. pombe*, for the runtime
normalization after the leakage check, unscored.

### 3.3 Fitted checkpoints: bounded local CPU fits v1, v2 and v3, 2026-09-22 (all finished; A misses the accuracy target in all three, for different reasons)

The recorded rows above all use the 20-step *S. cerevisiae* smoke
checkpoint (chain-only, no gene-free windows: 208,796 chains on chr V's
3,357 genes). The first fit of record is a **bounded local CPU fit** on the
two train species whose pinned sources are on this host, launched this
tick (code `4c43819`, run script and config to be recorded in
`smoke-local-20260920/fit-cpu-v1/`; the gagarin GPU request lenin-0083
stays open for the fit under the cap and the GPU half):

- Sources: *S. cerevisiae* (development chromosome I, `NC_001133.9`) and
  *C. elegans* (development chromosome V, `NC_003283.11`, the measured
  metazoan row); `benchmark/leakage_check.py` run first (0 violations),
  source MD5s verified by the loader against the committed summaries.
- Windows: `max_window` 12,288; 24,601 windows loaded — 19,708 train,
  4,893 dev, of which 548 are gene-free background tiles of 2,048 bases
  (148 = every *S. cerevisiae* candidate, 400 of *C. elegans*), drawn
  under seed 0, half reverse-complemented.
- Fit: 1,500 Adam steps, batch 8 (12,000 draws, ~29 Mb sampled, ~0.6 of
  one pass over the train windows), lr 3e-4, gradient clip 1, fast
  delayed-entry kernel, float64; checkpoint selection every 100 steps on
  256 seeded dev windows (`dev_windows_max`); one thread pinned to one
  core.
- Cost from the two probe runs at the same config (10 steps / 200 dev
  windows): load 21 s CPU; **4.0 s per step of 8 windows** (0.5 s per
  ~2.4 kb window forward + backward, i.e. ~0.2 ms per base); evaluation
  **0.13 s per dev window** (the full 4,893-window chr V set is ~10 min
  per evaluation, which is why selection subsamples); RSS 1.9 GiB. The
  1,500-step run projects to ~1.7 CPU-h on one core (6,000 s of steps,
  ~500 s of evaluations, ~20 s load), inside the tick-to-tick interval.
- Initial loss: 963 nats per window on the first batch (background tiles
  dominate: the all-intergenic numerator against a partition that at
  initialisation spreads mass over every chain), 431 train / 440 dev
  after 10 steps.

Run of record (`smoke-local-20260920/fit-cpu-v1/` in the relay artifacts:
`run.sh`, `config.json`, `leakage_check.out`, `source_md5.txt`,
`train.out`, `run_manifest.json`, `train_time.txt`, `run.log`,
`best_pt.sha256`): exit 0 after **1,500 steps in 5,833 CPU-s = 1.62 CPU-h**
on one core (`/usr/bin/time -v`: user 5,631.9 s + sys 201.4 s, wall
1:37:17, 99 % of one core, peak RSS 2.14 GiB; manifest timing: load 20.7 s,
fit 5,811 s of which the 15 evaluations 467 s). Train NLL per window
40 / 40 / 24 / 45 / 71 / 46 / 55 / 41 / 51 / 16 / 48 / 36 / 17 / 56 / 35 at
steps 100–1,500; dev NLL on the 256 selected windows 129 / 60 / 42 / 45 /
66 / 58 / 49 / 46 / **38** / 48 / 51 / 52 / 70 / 59 / 44, so `best_step` is
**900** (`best.pt` sha256 `181adad7…4cf88b`, 1.8 MB, host-local). Both
curves are noisy at batch 8 and the run covers ~0.6 of one pass; it is a
bounded fit, not a tuned one.

**Scoring pipeline dry run** (the 2026-09-22T03Z tick, interim step-300
checkpoint, one thread on a different core from the fit): `measure
--profile chromosome --segments 19 --overlap 4096 --gff-out` on
*S. cerevisiae* chr I (230,218 bases) then `benchmark/score.py --seqids
{chr I} --genome --declaration`; both exit 0, `score-dry/` in the
artifact directory. Chr I costs **12.79 CPU-s/Mb** end to end (encoder
1.30 s, decode 1.60 s, peak RSS 0.96 GiB; the scorer itself 0.17 CPU-s).
The interim numbers were a pipeline check: nucleotide sensitivity 0.963 /
precision 0.456 (MCC 0.44), 662 predicted loci against 94 reference (84
TP, 7 fusions), transcript F1 0.07, and 320 predicted introns against 3
reference of which 316 are non-GT-AG (engels-0098 and stalin-0102
reproduced these from the records and confirmed the motif tables had
moved in the canonical direction by step 300).

**Development-chromosome scores under the finished `best.pt`** (step 900;
`score-final/` in the artifact directory, same `measure --checkpoint`
→ `score.py` pipeline, code `aeed615` whose model source digest equals the
fit's `4c43819`, one thread on core 2, all exits 0):

| | *S. cerevisiae* chr I | *C. elegans* chr V |
|---|---|---|
| bases / tiles | 230,218 / 76 | 20,924,180 / 3,420 |
| end-to-end CPU-s (CPU-s/Mb) | 2.78 (**12.09**) | 185.4 (**8.86**) |
| encoder / decode / preprocess / io / output CPU-s | 1.26 / 1.48 / 0.033 / 0.008 / 0.0002 | 77.6 / 105.0 / 2.33 / 0.44 / 0.015 |
| peak host RSS | 0.96 GiB | 2.05 GiB |
| reference transcripts (loci) / introns | 94 (94) / 3 | 6,766 (4,995) / 22,824 |
| predicted transcripts | **6** | **287** |
| nucleotide sens / prec / F1 / MCC | 0.039 / 0.356 / 0.071 / 0.020 | 0.470 / 0.586 / 0.522 / 0.461 |
| locus TP / FP / FN; fusions; splits | 3 / 3 / 91; 2; 0 | 269 / 18 / 4,726; **238**; 121 |
| exact transcripts / exact CDS exons | 0 / 0 | 0 / 67 |
| predicted introns: GT-AG / GC-AG / AT-AC / other (`score.py` categories) | 31: 1 / 0 / 0 / 30 | 10,261: 586 (162 TP + 424 FP) / 43 / 10 / **9,622** (9,675 non-GT-AG) |
| predicted chain span median / max (bases) | 11,212 / 15,989 | **81,233 / 842,820** |
| predicted intron length median / max (bases) | 845 / 12,028 | 1,336 / 52,199 |

Cost: the raw stage sums (12.09 and 8.86 CPU-s/Mb) are below the
15 CPU-s/Mb figure, but that figure is bound to the cost-baseline runner
(cost-baseline 5.2) and this machine is 3.21× faster on the S. pombe
AUGUSTUS normalization (section 3.2), so the portable ceiling here is
51.41 / 11 = **4.67 CPU-s/Mb**. Under the fitted weights **both rows
miss it**: chr I at 12.09 / 12.42 / 14.64 (stage sum / user / user +
system) is 2.6× / 2.7× / 3.1× over, chr V at 8.86 / 7.77 / 8.89 is
1.9× / 1.7× / 1.9× over; machine-normalized to the baseline runner the
stage sums read **38.8 and 28.4 CPU-s/Mb against 15** (stalin-0103).
The fitted checkpoint therefore leaves the section 3.2 CPU verdict
unchanged (miss by ~2×; no positive CPU allowance for B), and only
memory passes (peak RSS ≤ 2.05 GiB against 8 GB). chr V is cheaper than
chr I per Mb because the decode traces back 287 rather than hundreds of
chains; the encoder share is unchanged from section 3.2. The fitted
S. pombe normalization row is below (v2).

Accuracy: **the fitted checkpoint does not produce gene-by-gene
structure, so candidate A misses the accuracy target in this fit.** On chr
V it emits 287 chains of median 81 kb that run through many reference
genes (238 fusions, 121 splits) joined by long non-canonical introns
(median 1.3 kb, 94 % non-GT-AG), which is why nucleotide sensitivity
reaches 0.47 while no transcript and only 67 of 28,569 CDS exons are
exact; on the intronless-dominated yeast chr I it emits almost nothing.
The interim step-300 checkpoint over-predicted instead (662 chains on chr
I). Working hypothesis (not yet established as the sole cause): both
behaviours are consistent with the run records' supervision balance. In
this fit's training windows every non-CDS base is supervised as *intron*
except the 10-base flanks (`FLANK` of the admission audit) and the
gene-free background tiles (2,048 bases each). Replaying the fit's seeded
draws (engels-0099), the 12,000 window draws through step 1,500 contained
269 background draws = 550,912 bases plus 234,620 flank bases, i.e.
**785,532 intergenic bases of 26,259,772 sampled (3.0 %)**; through the
selected step 900, 168 background draws = 344,064 bases plus 140,640
flank bases, 484,704 of 15,628,525 (3.1 %), against 8.76 Mb CDS and
6.39 Mb intron. (The 548 loaded tiles, 1.12 Mb, are the train-plus-dev
inventory, not what was sampled; 106 of them are dev tiles.) So the
encoder saw almost no signal separating intergenic from intronic
sequence, and the chromosome decode stitches neighbouring genes together
through intergenic stretches the emissions score as intron. The duration
factor alone does not prevent this (stalin-0102: −305 nats at 1,000
bases for the step-300 tables, yet the final decode's median intron is
1.3 kb). Sampling is the first thing to test, by the controlled revision
below scored the same way; other model behaviour is not ruled out until
that is done.

**Proposed revision before any B allowance** (task goal, last bullet):
the section 2 adjacent-gene / wide-flank training increment — windows
that carry their real intergenic neighbourhood with the neighbouring
genes' CDS supervised (or masked) rather than clean gene-plus-10-base
windows — plus a background share set from the genome's intergenic
fraction rather than a fixed 400 tiles per species, and a second bounded
fit (v2) scored the same way.

**Fit v2 launched 2026-09-22T05:13Z** (`fit-cpu-v2/`, code `61a9480`):
the loader increment is `context` (`model.a.dataset.iter_windows`): a
clean chain window is widened by up to `context` bases per side, clipped
at the nearest other gene's *CDS span* and the sequence ends, so the
unchanged support mask supervises the added bases as `U` with no loss or
decoder change. The guarantee is coding-span-free, not annotated
intergenic: the neighbour inventory is the CDS-bearing transcripts' outer
CDS endpoints, so UTRs (focal or neighbouring), ncRNA genes and CDS-free
pseudogenes can fall inside the added context — on the v2 train pool
48,213 of 3,965,426 added yeast bases and 3,964,153 of 12,049,332 added
*C. elegans* bases overlap raw `gene`/`pseudogene` spans, summed over
windows (engels-0100; e.g. the 163-base snoRNA `CELE_Y74C9A.6` lies
entirely inside `NM_058260.8`'s added context). The `U` share below is
therefore the support label's share, not a measured strictly-intergenic
annotation fraction; supervising a UTR or ncRNA as `U` is not evidence of a
wrong protein-coding target, and the stronger raw-feature exclusion
already used for background tiles is left as a knob for a later fit rather
than a change to this one. The manifest now records the
sampled-base composition (`actual.composition`: CDS / intron / intergenic
bases, background draws) so the supervision balance no longer has to be
replayed. Config: identical to v1 except `context: 512` and
`background_windows: 800` per species (148 available in yeast). Loader
probe of the train pool: `U` share 3.0 % (v1) → 20.2 / **30.3** /
39.1 % at context 256 / 512 / 1,024, mean window 2,202 → 2,654 / **2,970**
/ 3,405 bases; 512 matches the 3,072-base encoder core and brings the `U`
share near the genomes' own intergenic fractions at ~1.35× v1's per-step
cost. Leakage check 0 violations; 24,926 windows (19,975 train incl. 773
background / 4,951 dev incl. 175 background; 948 background in all); train
pool 59,330,766 bases = 23,891,418 CDS + 17,457,446 intron + 17,981,902
`U` (30.3 %). Note that v2 shares v1's dev *settings* (reserved
chromosomes, seed, 256-window limit), not its dev *subsample*: the loader
changes the dev pool (4,893 → 4,951 identities, 4,839 common) and hence
the seeded selection — only 37 of the 256 selected dev windows are common
to v1 and v2, and the retained chains carry context (stalin-0104). Raw
dev NLLs of v1 and v2 are therefore evaluated on different windows and are
not a loss change on one fixed set; the whole-chromosome scoring of chr I
and chr V under each `best.pt` is the common comparison.

**Fit v2 finished 2026-09-22T07:26Z** (`fit-cpu-v2/`: `run.log`,
`train.out`, `run_manifest.json`, `train_time.txt`, `best_pt.sha256`):
exit 0 after **1,500 steps in 7,964 CPU-s = 2.21 CPU-h** on one core
(`/usr/bin/time -v`: user 7,663.2 s + sys 301.1 s, wall 2:12:50, 99 % of
one core, peak RSS 2.09 GiB; manifest: load 24.1 s, fit 7,939 s of which
15 evaluations 659 s; 5.3 s per step against v1's 4.0 s, 1.37×). Train
NLL per window 181 / 90 / 33 / 62 / 152 / 61 / 14 / 32 / 18 / 79 / 88 /
36 / 41 / 93 / 15 at steps 100–1,500; dev NLL on the 256 selected v2
windows 284.9 / 68.8 / 56.0 / 90.9 / 52.0 / 76.3 / 54.6 / 45.6 / **45.3**
/ 48.9 / 74.1 / 56.7 / 48.6 / 58.7 / 50.7, so `best_step` is **900** as in
v1 (`best.pt` sha256 `99f773ed…2a8cf6`, host-local); the dev values are
on different windows from v1's (above) and are not compared. The
manifest's `actual.composition` over the full run: 12,000 accepted draws,
**35,305,664 sampled bases = 14,070,834 CDS (39.9 %) + 10,439,449 intron
(29.6 %) + 10,795,381 `U` (30.6 %)**, the `U` bases being 9,561,661
added context + 490 background draws (1,003,520) + 230,200 flank bases;
this equals stalin-0104's seeded replay (30.58 %) to the base. The
prefix through step 900 is not recorded separately. Against v1's
785,532 `U` of 26,259,772 (3.0 %), the v2 checkpoint saw ten times the
`U` share with the same optimizer, steps, batch and seed.

**Development-chromosome scores under the v2 `best.pt`** (step 900;
`fit-cpu-v2/score-final/`, same pipeline as v1, code `4d766b0` whose
model source differs from the fit's `61a9480` only by docstrings, one
thread on core 2, all exits 0; `summarize.py` there reproduces every
number from the records; v1 in parentheses):

| | *S. cerevisiae* chr I | *C. elegans* chr V |
|---|---|---|
| end-to-end CPU-s (CPU-s/Mb) | 2.89 (**12.54**) (12.09) | 185.8 (**8.88**) (8.86) |
| encoder / decode / preprocess / io / output CPU-s | 1.32 / 1.53 / 0.034 / 0.008 / 0.001 | 76.9 / 106.1 / 2.26 / 0.44 / 0.054 |
| peak host RSS | 0.96 GiB | 2.05 GiB |
| predicted transcripts | **214** (6) | **13,410** (287) |
| nucleotide sens / prec / F1 / MCC | **0.931 / 0.742 / 0.825 / 0.746** (0.039 / 0.356 / 0.071 / 0.020) | **0.531 / 0.546 / 0.538 / 0.468** (0.470 / 0.586 / 0.522 / 0.461) |
| locus TP / FP / FN; fusions; splits | 86 / 128 / 8; 2; 1 (3 / 3 / 91; 2; 0) | 3,804 / 9,606 / 1,191; **4**; **1,827** (269 / 18 / 4,726; 238; 121) |
| exact transcripts (sens / prec) / exact CDS exons | **69** (0.734 / 0.322) / 69 of 97 (0 / 0) | **101** (0.020 / 0.008) / 239 of 28,569 (0 / 67) |
| predicted introns: GT-AG / GC-AG / AT-AC / other | 28: 1 / 0 / 0 / 27 (31) | 2,290: 165 (68 TP + 97 FP) / 7 / 4 / **2,114** (10,261: 586 / 43 / 10 / 9,622) |
| donor / acceptor TP (F1) | 0 / 0 | 168 (0.014) / 306 (0.025) (613 / 751) |
| start / stop codon TP (FP) | 77 (137) / 79 (135) (1 / 0) | 1,412 (11,996) / 1,377 (12,029) (7 / 3) |
| chains by CDS exon count 1 / 2 / 3 / 4+ | 197 / 12 / 4 / 1 | 11,535 / 1,561 / 242 / 72 |
| predicted chain span median / max (bases) | 347 / 4,416 (11,212 / 15,989) | **294** / 12,343 (81,233 / 842,820) |
| predicted intron length median / max; count at 20–30 bases | 21 / 403; 25 of 28 | **34** / 6,785; 1,078 of 2,290 (1,336 / 52,199) |
| predicted chains with CDS < 300 bases | 92 of 214 | 6,827 of 13,410 |

Cost: unchanged by the checkpoint (12.54 and 8.88 CPU-s/Mb stage sum
against 12.09 and 8.86; decode traces back 13,460 chains on chr V at the
same 106 s as v1's 324, so traceback is not the decode's cost); the
section 3.2 CPU verdict stands — 2.7× and 1.9× the 4.67 CPU-s/Mb
portable ceiling (40.3 and 28.5 machine-normalized against 15), no
positive CPU allowance for B; memory passes. **The fitted S. pombe
normalization row** (`pombe-normalization/A-v2/`, v2 `best.pt`, same
command, core and float64 as the smoke row of section 3.2, 19 segments
per strand, S. pombe unscored): 12.5718 Mb in 112.14 stage / 99.93 user
/ 14.00 system CPU-s — **8.92 / 7.95 / 9.06 CPU-s per genome Mb** (smoke
row 9.22 / 8.23 / 9.37; stages preprocess 0.11, encoder 3.73, decode
5.07), peak RSS 1.31 GiB, 9,555 chains against AUGUSTUS's 4,452 genes.
Against AUGUSTUS 51.43 user + system CPU-s/Mb on the same core that is
**1/5.7 of AUGUSTUS against the 1/11 target** (1/6.5 user-only, 1/5.8
stage sum; 29.1 CPU-s/Mb machine-normalized against 15): the checkpoint
moves the row by 3 %, and the normalization verdict of section 3.2
stands at a 1.9× miss (1.7× user-only). Cost 0.032 CPU-h.

Accuracy, two findings. **The fusions are sharply reduced.** v1's fault —
287 chains of median 81 kb fusing 238 loci on chr V, near-silence on
chr I — is nearly absent under v2: chr V chains have median span 294
bases and 4 fusions (chr I 2); chr I hits 86 of 94 loci at nucleotide MCC 0.75 and reproduces
69 of its 91 single-exon transcripts exactly, start and stop codons
included (77 / 79 exact against 1 / 0). This compares two loader
configurations, not the `U` fraction alone: the v2 intervention raises
the `U` share of the bases seen through the selected step 900 from 3.1 %
to 30.6 % (stalin-0106's replay of the 7,200-draw prefix: 484,704 of
15,628,525 against 6,486,097 of 21,177,237 bases), but it also changes
every window's context geometry, the background pool, which chains are
drawn, the length exclusions, and the development windows the checkpoint
was selected on (stalin-0104). The common-chromosome result supports the
v1 working hypothesis (fusions from an unsupervised intron-versus-`U`
emission) more strongly than the draw replay did; it does not isolate
the `U` share as the sole cause or show the supervision defect fully
removed. **A still misses
the accuracy target in this fit, on splicing and on precision.** The v2
checkpoint predicts almost no real intron: 2,290 introns on chr V of
median 34 bases (1,078 at 20–30 bases, `min_intron` 20), 2,125 non-GT-AG
(2,114 `other` + 7 GC-AG + 4 AT-AC), 68 matching a reference intron, donor / acceptor F1 0.014 / 0.025 against
22,824 reference introns; 11,535 of 13,410 chains are single-exon, so
multi-exon genes come out as several single-exon chains (1,827 splits,
239 exact CDS exons of 28,569; transcript sensitivity 0.020 at locus
sensitivity 0.76). And it over-predicts: `score.py` reports 9,606 of
13,410 chr V loci and 128 of 214 chr I loci as FP (unmatched after
one-to-one greedy matching on shared same-strand CDS bases), about half
the chains carrying under 300 CDS bases. Splitting the FP by whether the
prediction shares any same-strand reference CDS base (stalin-0106,
reproduced; `fit-cpu-v2/score-final/README.md`): chr V 6,376 share none
and 3,230 overlap a reference CDS but lost the match — extra fragments
of genes already hit, the same fragmentation the splits measure; chr I
127 and 1. The no-overlap subset is the candidate set for short open
reading frames the emissions accept as coding once the `U` state is
available, but it has not been checked against full gene spans or
opposite-strand genes, so "overlaps no gene" is not established for it.
The tiny introns are consistent with the decoder bridging a frame break
inside what it scores as coding at the cheapest allowed length rather
than emitting learned splice events — a hypothesis from the length
distribution, not a measurement of the mechanism; the motif
tables that moved in the canonical direction during v1 (stalin-0102) have
not become decisive splice-site emissions in either bounded fit (1,500
steps of batch 8, ~0.6 of one pass, both best at step 900). Whether the
cause is the training budget or the encoder's splice-site receptive field
is not decided by this run.

**Proposed revision before any B allowance** (task goal, last bullet,
superseding the one above now that it is tested): a longer fit under the
same loader and config, scored the same way — the GPU grant requested in
lenin-0083 is the instrument for it — and no decoder change (intron
duration floor above `min_intron`, short-ORF penalty) until that longer
fit's measured development results show whether the splice-site emissions
sharpen on their own. The trigger is the measured dev NLL and chr I / chr
V scores, not a nominal number of passes: sampling is uniform *with
replacement*, so a longer draw stream does not mean the model has seen
every train window (below). Section 6.1's accuracy column carries the v2
numbers as the current measured entry; they are not the reported accuracy
of A.

**Fit v3, the longer fit, started 2026-09-22T11:07Z and is running**
(`relay/artifacts/T-human-014/smoke-local-20260920/fit-cpu-v3/`,
`README.md` there records the full delta). The GPU grant of lenin-0083 is
still open, so this is the local CPU fallback: same loader, dev
reservations, optimizer, batch, base `lr`, seed, dev subsample and
launcher as v2, on one core, with **three** declared changes
(engels-0105 / stalin-0110 P3: the third was previously described as a
neutral bookkeeping setting, and it is not):

1. `steps` 1,500 → 3,000, i.e. 24,000 planned draws for 19,975 train
   windows. That ratio, 1.2015, is *draw-equivalent epochs*, not pool
   coverage: the sampler draws uniformly **with replacement**, so
   replaying the declared seed and batch shape gives 14,027 distinct
   windows among the 24,000 planned draws (70.22 %), leaving 5,948 train
   windows never attempted; v2's 12,000 draws reached 9,035 (45.23 %).
2. A cosine learning-rate schedule with a 150-step warmup decaying
   3e-4 → 1.5e-5 (`model.a.train.lr_at`, unit tests). The schedule is
   part of "longer" rather than a separate accuracy intervention: at a
   constant step size both fits selected step 900 of 1,500 and then
   oscillated (v2 dev NLL 45.3 → 74.1 → 48.6 → 58.7 → 50.7 over steps 900
   to 1,500), so further steps at 3e-4 would most likely keep bouncing.
3. `eval_every` 100 → 150. This is not only reporting cadence: only
   evaluated steps can become `best.pt`, so v3 selects from 20 candidate
   checkpoints on a 150-step grid where v2 selected from 15 on a 100-step
   grid, and over the shared first 1,500 steps the two grids offer 10
   and 15 candidates at different steps. The fixed dev subsample keeps
   individual NLL values comparable; it does not remove this
   checkpoint-selection difference.

Any v2 → v3 difference is therefore attributable to the combined
step-budget, learning-rate-schedule and checkpoint-cadence revision, not
to step count alone. The earlier claim here that `eval_every` 150 keeps
the evaluation share of fit CPU at v2's ~8 % was wrong: from the v2
manifest, evaluation cost `E` = 659.18535381 CPU-s against inclusive fit
`F` = 7,938.915226944 CPU-s, i.e. 8.303 %, and holding unit costs fixed
while doubling fitting updates and taking 20/15 as many evaluations
projects `(20/15·E)/(2·(F−E)+20/15·E)` = **5.693 %**, not ~8 %. That is a
projection; v3's actual share is reported from its own manifest after it
exits. Projected ~4.4 CPU-h at v2's measured 5.3 CPU-s/step. No decoder
change, no loader change, and nothing is claimed from it until it exits 0
and is scored on chr I and chr V the same way.

The scoring inputs for that comparison are **pre-registered**: the
`score-final/` directory under the v3 artifact now holds v2's
`run_final.sh` with the `fit-cpu-v2` → `fit-cpu-v3` paths substituted,
v2's `declaration.yaml` with only the `model:` line changed, and v2's
seqid lists, committed while the fit is still running — before fit
completion and final checkpoint selection, and before any v3 chromosome
scoring — so the chromosomes, tiling, dtype, core and scorer invocation
cannot be chosen after seeing the result. This paragraph previously said
"before any v3 checkpoint exists", which was wrong: the training loop
saves `best.pt` immediately on every improving development evaluation, so
interim checkpoints predate the commit (engels-0108, stalin-0113). What
is fixed in advance is the scoring configuration, not the absence of a
model. The launcher has since gained an explicit completion gate (the
fit's recorded `train exit=0`, the GNU-time `Exit status: 0`, no running
v3 training process, a non-empty `best.pt`) and per-workload failure
propagation, because `taskset` only sets CPU affinity and neither waits
for the fit nor makes the core idle, and the v2 launcher reported success
even when its measurement and scoring commands failed (engels-0108 P3,
stalin-0113 P2). Every workload argument remains identical to v2's, an
idle core remains an operator precondition recorded with the run, and
they were then run unchanged after the fit exited 0.

**Fit v3 finished 2026-09-22T15:27:57Z, `train exit=0`, GNU-time
`Exit status: 0`.** 3,000 steps in 15,562.6 CPU-s inclusive of
evaluations = **4.32 CPU-h** on one core (user+system 15,588.05 s =
4.33 CPU-h, wall 4:19:59), 5.188 CPU-s/step against v2's 5.293, peak RSS
2.11 GiB. The actual evaluation share is **5.582 %** (868.75 /
15,562.64), against the 5.693 % projected from v2 under the
constant-unit-cost assumption. The accepted-draw composition is
40.0 % CDS / 29.6 % intron / 30.3 % intergenic over 71,218,348 sampled
bases (24,000 accepted of 24,000 attempted, 955 background draws) —
v2's balance at twice the volume. The selected checkpoint is
**step 3,000 of 3,000**, dev NLL **28.640** on the fixed 256-window
subsample against v2's 45.335 at step 900; `best.pt` sha256 `8a32c93e…`.
The ~10-NLL oscillation of the first half disappears once `lr` falls
below ~1e-4 (last seven evaluations in a 28.6–32.4 band), and the best
point is the last step, so the step budget rather than convergence ended
this fit.

**Scored on chr I and chr V with the pre-registered launcher**
(`fit-cpu-v3/score-final/`, gate passed, all four workloads exit 0,
`summarize.py` reproduces the tables; measurement commit `85d534c`,
whose only difference from the fit's `6c7aac1` inside `model/` and
`benchmark/` is the cosine config validation and its tests):

| metric | chr I v2 → v3 | chr V v2 → v3 |
|---|---|---|
| nucleotide F1 | 0.825 → **0.891** | 0.538 → **0.598** |
| nt sens / prec | 0.931/0.742 → 0.912/0.870 | 0.531/0.546 → 0.533/0.680 |
| locus F1 | 0.558 → **0.786** | 0.413 → **0.603** |
| predicted chains (median CDS span) | 214 (346 b) → 107 (1,149 b) | 13,410 (294 b) → 6,156 (621 b) |
| fusion / split | 2/1 → 3/0 | 4/1,827 → **312**/910 |
| exact transcripts (sens) | 69 (0.734) → 64 (0.681) | 101 (0.015) → 70 (**0.010**) |
| donor / acceptor F1 | 0 → 0 | 0.014/0.025 → 0.070/0.081 |
| GT-AG introns TP | 0 → 0 (3 reference) | 68 → **433** (22,620 reference) |
| median predicted intron | 21 b → 81.5 b | 34 b → **615 b** |
| CPU-s/Mb (RSS) | 12.54 → 12.38 (0.96 GiB) | 8.88 → 8.88 (2.05 GiB) |

Reading of that table, stated as what it supports:

1. **The gain is precision, not sensitivity.** Nucleotide sensitivity is
   flat on both chromosomes; precision moves 0.742 → 0.870 and 0.546 →
   0.680. Chains halve and median spans roughly double: A emits fewer,
   longer chains, it does not find more coding sequence.
2. **Intron *length* is fixed, intron *placement* is not.** The v2 defect
   of piling introns at the 20 b floor is gone (median 34 b → 615 b, the
   right regime for *C. elegans*) and correct GT-AG introns rise 68 →
   433, but 22,187 reference GT-AG introns are still missed, 3,986
   predicted introns are non-canonical, and donor/acceptor F1 reach only
   0.070/0.081. Splicing remains the open failure.
3. **The split defect has become a fusion defect.** Splits fall 1,827 →
   910 while fusions rise 4 → 312 and chr V exact transcripts fall 101 →
   70 (sens 0.015 → 0.010). Longer chains merge neighbouring genes. The
   locus-F1 gain is real and is bought partly with fusion, so it is not a
   clean improvement in gene-level structure.
4. **Cost is unchanged** (12.38 and 8.88 CPU-s/Mb, decode still the
   larger stage on chr V at 105.48 of 185.80 CPU-s): same architecture
   and tiling, different weights. The section 5 CPU verdict stands and
   **no positive CPU allowance for B follows from this run**.
5. **A still misses the accuracy target.** Exact-transcript sensitivity
   on chr V is 0.010. Three fits (v1, v2, v3) have each fixed the
   previous one's dominant geometric defect and exposed the next:
   fusion → sub-minimal introns and over-prediction → plausible intron
   lengths with wrong splice sites and fusion. The next revision is the
   decoder, with the section 6 CPU ceiling in scope, since the remaining
   failures (splice-site placement, locus boundaries) are decoding
   decisions rather than a shortage of fitting steps — though the
   monotone dev tail means a still longer fit is not ruled out either,
   and the two are separable only if measured separately.

`docs/cost-baseline/measured.tsv` now carries four **fitted** candidate-A
rows — v2 and v3 on *S. cerevisiae* chr I and on *C. elegans* chr V — in
the 3.2 user-only convention, with the nucleotide, locus and transcript
F1 columns filled from `benchmark/score.py`. The earlier candidate-A rows
there are smoke-checkpoint cost rows and stay marked unscored.

CPU accounting for this section: fit v1 1.62 CPU-h, fit v2 2.21 CPU-h,
fit v3 4.33 CPU-h (measured, user+system), final scoring 0.05 CPU-h each
(chr I 3 s, chr V ~190 s incl. scorer), three scoring runs so far,
probes and dry run ~0.35 CPU-h; all local, cluster CPU-hours 0, GPU-hours
0 (lenin-0083 still open).

## 4. Budget and caps

Model ≤ 5 M parameters (455,841 ✓) and ≤ 8 GB. Fitting uses train species
only. T-human-014 compute cap: ≤ 24 GPU-hours; one request may exceed 24 h
wall clock only with a coordinator `decision`. Sampled bases and repeats per
fitting run are declared in the run manifest before the run; actual GPU/CPU
hours are recorded in this task's log after every run.

## 5. Section 6.1 sensitivity table (measured)

_Partially filled. CPU regime, one core, yeast chromosome I end to end
(3.2 item 6; both strands, 12,288-base windows, GFF3 written, per genome
Mb), after revision steps 1 (vectorised featurizer, `48f8a2b`), 2a
(per-phase decode operands, `6cae82f`) and 2b (sparse-predecessor scan,
`d9ea0f8`). Numerators are labelled because they differ (engels-0086):
stage sum of `process_time()` / whole-process user only (the
`cost-baseline.md` 3.2 convention) / whole-process user + system.
**Overlap 4,096 / batch 16 (defaults): 17.1 / 18.1 / 19.7 CPU-s per
genome Mb; overlap 4,096 / batch 64: 14.2 / 14.4 / 16.8; overlap 2,048 /
batch 64: 12.0 / 12.8 / 14.7; overlap 0 / batch 64: 10.7 / 11.9 / 13.3 —
against the 15 CPU-s/Mb budget the batch-16 default misses by 1.14–1.3×
on every numerator, the default overlap at batch 64 is inside on the
stage sum and the baseline's user-only numerator but not on user +
system, and overlap 2,048 and 0 at batch 64 are inside on all three**;
peak RSS ≤ 0.81 GiB inside 8 GB; I/O + output 0.1 CPU-s/Mb. Per oriented
megabase the stages are now featurizer 0.05, encoder 1.8, R = 3 tensor
Viterbi 2.9–5.2; what separates the batch-16 default from the budget is
the both-strand factor times 1.5× oversampling on the decode stage, whose
remaining cost is the per-step fixed dispatch of the Python scan loop.
The pure-Python reference decoder at 66 CPU-s per oriented Mb (item 3) is
kept only for parity runs. Revision step 3 (overlap through carried
scan state at seams, packed back-pointers, streamed emissions, encoder
context margin; the edge decoder exists for the true sequence ends)
gives **13.2 / 13.3 / 15.9** at 19 segments per strand with the
4,096-base seam overlap (section 3.2, `chromosome-margin/`): inside 15 on
the stage sum and the user-only numerator, over on user + system, whose
excess on this 230 kb chromosome is the page-fault cost of the 1 GB
working set. **The S. pombe normalization row (section 3.2,
`pombe-normalization/`) is measured: A at 19 segments 9.22 / 8.23 /
9.37 CPU-s per genome Mb on the three numerators, RSS ≤ 1.31 GiB,
against AUGUSTUS 3.5.0 at 51.41 user CPU-s/Mb on the same core (this
machine is 3.21× the cost-baseline runner on that command) — 1/5.5 of
AUGUSTUS against the 1/11 portable target, 30.0 CPU-s/Mb
machine-normalized against 15: the CPU regime misses the budget by
2.0× (1.8× on user only), and no positive CPU allowance for B follows.**
Segment count is exhausted as a lever (38 and 76 segments per strand
change the stage sum by −5 % / −3 % at 1.6× / 2.9× the memory), and
revision step 4, float32 decode, is measured at **8.27 / 7.94 / 8.43**
(`A-f32/`): a 10 % gain, miss 1.8× (1.7× on user only), with 103 of
121,937 chains flipped at the smoke checkpoint, so float64 stays the row
of record. The decode is per-step dispatch-bound (42 µs per step at
B = 38), so the remaining decode lever is a compiled scan step, worth
at most ~2 CPU-s/Mb; the encoder work reduction, measured together with
a fitted checkpoint's accuracy, is necessary for 1/11 on CPU.
**Accuracy column, current entry (section 3.3, fit-cpu-v2 step 900,
`fit-cpu-v2/score-final/`; v1 in parentheses)**: *S. cerevisiae* chr I
nucleotide 0.931 / 0.742 sens / prec (0.039 / 0.356), 214 predicted
transcripts against 94 with 69 exact (6, 0 exact); *C. elegans* chr V
0.531 / 0.546 (0.470 / 0.586), 13,410 predicted against 6,766 with 101
exact, 4 fusions and 1,827 splits (287, 0 exact, 238 fusions), 239 exact
CDS exons (67), 2,114 of 2,290 predicted introns in `score.py`'s `other`
category with donor F1 0.014 (9,622 of 10,261) — under the fitted weights
the two chromosomes cost 12.54 and 8.88 CPU-s/Mb (stage sum; v1 12.09 and
8.86), **2.7× and 1.9× the 4.67 CPU-s/Mb portable ceiling of this machine
(40.3 and 28.5 machine-normalized against 15)**, at ≤ 2.05 GiB: the CPU
regime still misses by ~2× as in section 3.2, memory passes. A misses the
accuracy target in both fits: v1 with fusions consistent with the
supervision-balance hypothesis (3 % `U`; sharply reduced, not shown
removed, under the context loader), v2 on splicing (almost no real
intron predicted) and precision (short single-exon over-prediction),
cause unresolved; the proposed revision is a longer fit under the same
loader before any decoder change. The fitted
S. pombe row (v2 `best.pt`) is 9.06 user + system CPU-s/Mb, 1/5.7 of
AUGUSTUS against 1/11 (29.1 machine-normalized against 15), a 3 % move
from the smoke row. Still missing: that fit (accuracy), CPU efficiency
(encoder work reduction, section 3.2), and the GPU regime with the 6.1
multi-worker decoder accounting (gagarin)._

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

Fast-kernel review findings (engels-0080, stalin-0081) and their resolution:

- **Grammar cache keyed by table number, not by the genetic code** (P2, both
  reviewers; the collision reproduced by both). `Grammar.get` keyed on
  `code.table`, so `TABLES[1]` and `TABLES[1].with_alternative_initiators()`
  shared one cached grammar and whichever ran first decided the other's
  initiator set (a legal `TTG` start dropped, or a forbidden one admitted,
  depending on call order). Now keyed by the frozen `GeneticCode` value.
  Regression `test_grammar_cache_keyed_by_genetic_code_value` runs both call
  orders against the reference and asserts `log Z = log 2` / `0` for the two
  variants after either warm-up. The ATG-only training default was not affected.
- **Fast-run manifests could not be reproduced from their recorded commit**
  (P2, engels-0080). Reconciled above (3.2, "Provenance") and in
  `fast-kernel/PROVENANCE.md`; the manifest now pins the executed source by hash
  and dirty state; the one-core fit is repeated at a clean committed source.
- **Memory units and RSS semantics** (P3, engels-0080). Sections 3.1 and 3.2
  restated in GiB with the raw KiB/MiB beside each figure, the cumulative
  process-RSS nature of the profile numbers stated, and the per-base figures
  re-derived as baseline-subtracted process increments. `measure` JSON keys are
  now `peak_host_rss_gib` / `peak_device_mem_gib` (the earlier
  `measure-chrI/measure_chrI.json` keeps the old key name; its value, 0.3106,
  was already GiB).
- **Unsupported GPU-cap extrapolation** (P3, engels-0080). Removed from 3.1; the
  direct 4,000–7,000 CPU-s/Mb figure stands and the GPU regime is pending.
- **Empty window raised `IndexError` in the fast twin** (P3, engels-0080).
  Symbol tensors are built with `dtype=torch.long`; `partition("", zeros(11, 0))`
  returns `0` like the reference. Regression `test_empty_window_matches_reference`.
- **`git status --porcelain` lines were stripped before slicing the path**
  (P3, engels-0081): an unstaged tracked modification (`" M path"`) lost the
  first character of its path in `source_dirty_files`. The slice now consumes
  the raw line; regression `test_dirty_paths_keep_porcelain_status_columns`
  covers unstaged, staged, untracked and mixed output. The dirty flag and
  source digest were unaffected; no recorded run had dirty paths.
- **"8 GB read as 8 GiB"** (P3, engels-0081): reverted — the charter's
  threshold stays 8,000,000,000 bytes (7.451 GiB); GiB is the reporting unit
  only (3.1 memory convention). No measured conclusion changes.
- stalin-0081's 216-case phase/split sweep (introns inside initiator and stop
  codons, three tables, `m` 1/2/4) is noted as independent parity evidence at
  ≤ 7.1e-15; not added to the suite here, since the fixture and random-lattice
  tests already cover it in under a second and the sweep takes six.
- **Stage ranking and the unbatched rate** (P3, engels-0082, corroborated by
  stalin-0083): 3.2 item 4 called the batched decoder the cheapest stage; it
  is the most expensive of the three (4.9 > 3.7 > 2.3 CPU-s/Mb), and the
  unbatched 30.0 is decode + traceback alone (three-stage total 35.9). Wording
  fixed; no number changes. engels-0082's 288-case and stalin-0083's 252-case
  multi-gene/ambiguous/padded-batch sweeps are noted as independent parity
  evidence at ≤ 1.3e-14, not added to the suite for the reason above.
- **Locus FP conflated two things** (P2, stalin-0106): 3.3 and the v2
  score README said the 9,606 / 128 FP loci "overlap no reference gene";
  `score.py:loci` counts predictions left unmatched after one-to-one greedy
  CDS matching, so extra fragments of a matched gene are FP too. Replayed
  the matching (stalin's read-only script, reproduced): chr V 6,376 share
  no same-strand reference CDS base + 3,230 overlap-but-unmatched; chr I
  127 + 1. Narrative replaced in the summary, 3.3, section 5 and the README;
  the no-overlap subset is labelled as unchecked against gene spans and the
  opposite strand. No score changes.
- **Causal scope of v1 → v2** (P3, stalin-0106): "the only change is the
  sampled `U` share" and "the fusions are gone" reworded — the comparison is
  of two loader configurations (context geometry, background pool, drawn
  chains, length exclusions and dev windows all change), it supports the
  supervision-balance hypothesis without isolating the `U` fraction, 4 + 2
  fusions remain, and the minimum-length-intron reading is a hypothesis from
  the length distribution. Count fix: chr V 2,125 non-GT-AG introns (2,114
  `other` + 7 GC-AG + 4 AT-AC), not 2,114. stalin-0106's step-900 prefix
  replay (v1 3.10 % / v2 30.63 % `U`, 484,704 / 6,486,097 `U` bases) is
  recorded in the v2 README as what `best.pt` saw.
- **Opening and README closings said the failure had "moved" off
  supervision balance** (P3, engels-0102 / stalin-0107): the summary, the
  section 5 verdict and the v2 score README's closing paragraph now say
  what 3.3 says — the fusion reduction supports the hypothesis without
  showing the defect removed; the remaining splicing / precision misses
  have an unresolved cause. The S. pombe normalization README no longer
  calls the 9,555 chains "the same short-chain over-prediction seen on
  chr V": S. pombe is unscored, so the count is output volume for the cost
  row and any over-prediction reading there is a hypothesis. Wording only;
  no numbers change.
- **Decoder stage alone exceeds the ceiling** (stalin-0107, conditional
  arithmetic on the fitted float64 S. pombe records): encoder 3.73 and
  decode 5.07 CPU-s/Mb against a whole-process portable ceiling of
  646.50 / 12.57 / 11 = 4.67; removing the encoder entirely leaves 5.19
  stage CPU-s/Mb, and decode is 1.08× the ceiling. Noted for the CPU
  revision: encoder work reduction alone cannot admit this configuration;
  the decoder itself (float32, fewer states, or another implementation
  change whose aggregate CPU-seconds are measured on one core) has to
  move too. The section 6.1 multi-worker accounting belongs to the GPU
  regime's wall time only: dividing wall time by worker count does not
  reduce aggregate CPU-seconds and cannot lower the CPU row
  (engels-0103 / stalin-0108). No new measurement.

- **A singleton cosine tail returned the full rate** (P3, engels-0105,
  corroborated by stalin-0110 over a 2,080-config sweep). With
  `warmup_steps == steps - 1` (and with `steps=1`) the decay phase has a
  single step, `lr_at` saw progress 0 and returned `config.lr` instead of
  the documented `lr * lr_min_factor` floor — 3e-4 rather than 1.5e-5 at
  v3's settings. `TrainConfig.from_dict` now rejects a cosine schedule
  with fewer than two post-warmup steps, so every loadable config meets
  the endpoint contract; `lr_at`'s docstring states the precondition and
  that a directly constructed `TrainConfig` is not covered. Two
  regression tests pin it, including a sweep asserting that the rejected
  set is exactly `warmup_steps == steps - 1` (64 of 2,080 over
  `steps` 1..64) and that every accepted config ends at the floor.
  V3's 3,000/150 configuration is unaffected and the running fit is not
  restarted.

- **The v2 → v3 comparison had three changed variables, not two, and the
  evaluation-share claim was wrong** (P3, engels-0105 / stalin-0110).
  `eval_every` also changes, 100 → 150, and only evaluated steps can
  become `best.pt`, so the checkpoint grid and candidate count differ
  (20 vs 15 overall; 10 vs 15 over the shared first 1,500 steps).
  Section 3.3 and the v3 README now describe the combined step-budget,
  LR-schedule and checkpoint-cadence revision. The claim that 20
  evaluations over twice the steps keep v2's ~8 % evaluation share is
  replaced by the reproduced arithmetic: v2 measured `E/F` =
  659.18535381 / 7,938.915226944 = 8.303 %, and at fixed unit costs
  `(20/15·E)/(2·(F−E)+20/15·E)` projects 5.693 %. Flagged as a
  projection; v3's actual share comes from its own manifest.

- **Draw-equivalent epochs are not pool coverage** (P3, engels-0105 /
  stalin-0110). Sampling is uniform **with replacement**, so v3's
  24,000 / 19,975 = 1.2015 is a draw-equivalent epoch count, not a
  guarantee that every train window is seen. Replaying the declared seed
  and batch shape reproduces both reviewers' counts exactly: 9,035
  distinct windows in v2's 12,000 draws (45.23 %) and 14,027 in v3's
  planned 24,000 (70.22 %), leaving 5,948 never attempted. Section 3.3
  and the v3 README now label the ratio accordingly, and the deferral of
  the decoder revision is keyed to the longer fit's measured development
  results rather than to "more than one pass". These are planned
  index-exposure counts, not observed accepted-window counts.

- **The scoring launcher reported success after failed workloads** (P2,
  stalin-0113). Only `set -u` was enabled, and each workload's exit
  status was printed and then replaced by the status of the `echo` that
  printed it, so a failed measurement still entered scoring and the
  wrapper could return 0 with both workloads broken — and a rerun in an
  existing output directory could score a stale prediction. v3's
  `score-final/run_final.sh` now captures each status, returns it, and
  skips scoring when the measurement fails. Reproduced with stalin's
  stubbed harness (no model or scorer runs): measurement failure exits 9
  with no `score` line, scoring failure exits 7, and only a clean pair
  prints `all workloads exit=0`. Every workload argument is unchanged
  from v2's launcher — verified by stripping comments, normalising the
  two version strings and diffing the workload lines: identical. The
  fitting launcher masks status the same way, so fit v3's completion is
  read from its recorded `train exit=` line and the GNU-time
  `Exit status`, not from the launcher's own status. No restart.

- **CPU affinity is not a completion gate or an idle-core guarantee**
  (P3, engels-0108, restated by stalin-0113). The claim that pinning the
  scoring launcher to the fit's core meant it could only run after the
  fit exited was wrong: `taskset` sets affinity, and two processes may
  share one CPU. `run_final.sh` now starts with an explicit gate — the
  fit's recorded `train exit=0`, the GNU-time `Exit status: 0`, no
  running v3 training process, and a non-empty `best.pt` — checked live
  at 15:07Z against the in-flight fit, where it exits 1 and creates no
  output directory. An idle core stays an **operator precondition**,
  recorded with the run rather than enforced by the script.

- **Pre-registration is before final checkpoint selection, not before
  any checkpoint** (P3, engels-0108 / stalin-0113). Section 3.3, the v3
  README and `score-final/README.md` said the scoring inputs were
  committed "before any v3 checkpoint exists". The training loop saves
  `best.pt` immediately on every improving development evaluation, so
  interim checkpoints predated the commit. All three now say **before
  fit completion and final checkpoint selection, and before any v3
  chromosome scoring**, which is the claim the commit actually supports:
  what is fixed in advance is the scoring configuration, not the absence
  of a model.

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

# The whole train panel. Species run sequentially (one genome resident at a
# time), and the fetched compressed genomes are the added scratch cost (below):
python3 -m model.a.coverage --sources <scratch-dir> --fetch --out coverage.tsv
```

Passing no `--species` selects all ten committed train species (including
mouse, maize and zebrafish). `coverage_report` processes species sequentially
and retains only counts, so only one genome is resident at a time, not the sum;
`--fetch` retains the downloaded source files, so the whole-panel disk cost is
scratch capacity for all ten compressed genomes at once. Whether disk or memory
is the binding laptop constraint depends on the host and is not settled by the
manifests alone (see the memory note below). An unknown `--species`, or an empty
manifest directory, is now a clean argument error in both modes rather than a
zero-row table, and `--fetch` progress goes to stderr so a redirected stdout is
a clean TSV.

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

**Memory note.** The manifests record a per-species `peak_rss_mb` for the
*standalone* `audit_species` call, sampled right after the audit; it ranges
**0.08–3.6 GB across all ten species (79.6–3554.6 MiB; the field is MiB on the
Linux producer despite its `mb` name)** and **1.9–3.6 GB for the four large
species alone (Xenopus, Zea, Danio, Mus)**. That is the audit high-water mark,
not the coverage loader's: `iter_windows` keeps the audit result live while it
re-parses the GFF, builds transcript indexes and oriented windows, so the full
coverage-process peak is higher by an unmeasured margin and is not reported here.
Because species run sequentially, only one genome is resident at a time, so the
resident set is bounded by the largest single species, not the panel sum — this
is what let the four large genomes run on the laptop (below), not a claim that
the audit figures bound the loader. The whole-panel `--fetch` additionally needs
enough scratch disk for the ten compressed genomes together, since `fetch`
retains each download; the manifests alone do not establish which of memory or
disk is the binding constraint on a given laptop.

This measures the current scope; it does not widen it. The four small train
species were checked out and reported on the laptop (lenin, 2026-09-19,
commit ef82af9, `max_window=None`); Arabidopsis and Drosophila — the two
smallest of the remaining six — were added on the laptop the same way (lenin,
2026-09-20; MD5-verified fetches, ~60 MB and ~53 MB compressed, 27 s and 13 s
audits, nothing committed). The four large genomes (Danio 1.45 Gb, Xenopus
1.45 Gb, Mus 2.7 Gb, Zea 2.18 Gb) were also completed on the laptop (lenin,
2026-09-20): each source `_genomic.{gff,fna}.gz` was fetched to scratch and
MD5-verified against its committed manifest (`fetch` re-checks the same digest
`verify_source` enforces before any window is yielded), audited on one core
(these four carry the largest standalone-audit peaks, 1.9–3.6 GB; see the memory
note above), and its FASTA discarded; nothing committed. The panel table below is therefore the
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
