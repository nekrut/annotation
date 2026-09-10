# Cost baseline: what existing gene predictors cost, and the budget ours must fit

Task [T-human-009](../relay/tasks/T-human-009.md). Owner: `marx`.
Status: **in review, revised after lenin's review** (2026-09-10; the six
fixes of `relay/messages/20260910T120853Z-lenin-0033.md` are in §5 and §6,
and the isoform-summed coding fraction is replaced by a measured genomic
union, `docs/cost-baseline/cds_union.tsv`). Every figure below carries
its source, and the measured rows were produced by the commands in §3.4 on
the runner described there. Two rows, EGAPx and Tiberius end to end, need
the cluster; the `alert` requesting it is
[marx-0033](../relay/messages/20260910T115307Z-marx-0033.md) (§5.4).

Two kinds of number appear here and they are kept apart. **Documented**
figures are quoted from a paper, a README or a project's own report; they
were produced on hardware we do not control, with settings we often cannot
see, and are normalized per megabase only where the source states the
genome size. **Measured** figures were produced for this document on one
machine with one stopwatch (`/usr/bin/time -v`) and are recorded together
with the accuracy of the same run, scored by the panel scorer
(`benchmark/score.py`), so cost and accuracy are never quoted from two
different runs.

## 1. Summary

| tool | class | documented cost | measured here | fits a laptop? |
|---|---|---|---|---|
| EGAPx | evidence pipeline (miniprot/ProSplign, STAR, minimap2, Gnomon) | 71 CPU-h / 3 h wall for the 144 Mb fly; 425 CPU-h / 5.5 h wall for the 1.1 Gb chicken; wants a 32-CPU, 256 GB machine (§2.1) | not attempted; below the stated minimum machine | no |
| Tiberius | CNN + LSTM + differentiable HMM, ab initio | GPU with ≥ 8 GB; batch sizes given for A100 80 GB, RTX 3090, RTX 2070; runtime "reduced by 30%" in 2.0.0 (§2.2) | not measurable on this runner (no GPU, Python 3.11 < 3.12) | GPU laptop only |
| BRAKER3 | GeneMark-ETP + AUGUSTUS training pipeline | 20 min test runs on a Xeon E5530; full genomes hours to days (§2.3) | not attempted this tick | CPU, hours to days |
| AUGUSTUS 3.5.0 | GHMM, ab initio | 6 min per 1.6 Mb on a 2.4 GHz PC (2003); 2 h 25 min per mammal on 48 threads (2024) (§2.4) | **57 to 169 CPU-s/Mb, 0.23 to 0.97 GB** on five panel genomes, one of them the 100 Mb *C. elegans* (§3) | yes |
| Helixer 0.3.7 | CNN + biLSTM, ab initio | GPU with 8 to 11 GB; 3 to 5 min demos on a GPU (§2.5) | Docker image only; no daemon on this runner | GPU laptop only |
| GeneMark-ETP | self-training GHMM + evidence | `--cores 32` in its own usage line; licence CC BY-NC-SA 4.0 (§2.6) | not attempted; non-commercial licence recorded | CPU, hours |
| KA/KS window test (our floor) | comparative statistic | 365 s wall for 32 genes, nearly all fetch; the test itself is milliseconds per window (§4) | measured in T-human-010 | yes |

### 1.1 Documented figures, normalized where the source allows

Sources are the primary papers and READMEs; the "reviews" column names
which Phase 1 review artifacts carry the figure (see
`docs/review/README.md` §2 for the merged table). CPU-s/Mb for multi-thread
runs is wall clock × threads ÷ genome size and is an **upper bound**,
because thread occupancy is never reported. Genome sizes are from
`benchmark/panel.tsv` where the species is on the panel.

| tool | run | wall clock | hardware | CPU-s/Mb or GPU-s/Mb | memory | source | reviews |
|---|---|---|---|---|---|---|---|
| EGAPx | fly, 144 Mb, 1 RNA-seq run | 3 h | AWS Batch, 8 to 32-CPU instances | 1,775 CPU-s/Mb (71 CPU-h) | up to 256 GB machine | [README](https://github.com/ncbi/egapx/blob/main/README.md) | 002, 003 |
| EGAPx | chicken, 1.1 Gb, 20 RNA-seq runs | 5.5 h | same | 1,391 CPU-s/Mb (425 CPU-h) | same | README | 002 |
| EGAPx on Galaxy | 1,409 jobs, 290 genomes | | | 416,000 CPU-h total; **45% failure**; 24% CPU efficiency | | charter, `nekrut/scalingPaper` (not readable by agents) | 002, 005, 012 |
| BRAKER3 | *A. thaliana*, 119 Mb | 5 h 37 min | 48 threads, Xeon E5-2650 v4; RNA-seq alignment excluded, AUGUSTUS training included | ≤ 8,160 CPU-s/Mb | | [10.1101/gr.278090.123](https://doi.org/10.1101/gr.278090.123) | 003, 004, 005 |
| BRAKER3 | mouse, 2.73 Gb | 64 h 16 min | same | ≤ 4,070 CPU-s/Mb | | same | 003, 004, 005 |
| BRAKER3 | 33-species panel, mean | 2,170 min | 72 threads | | | [10.64898/2026.04.24.720536](https://doi.org/10.64898/2026.04.24.720536) | 003, 004 |
| GeneMark-ETP | fly, 144 Mb | 3 h | 64 cores; HISAT2/StringTie2 excluded | ≤ 4,800 CPU-s/Mb | | [10.1101/gr.278373.123](https://doi.org/10.1101/gr.278373.123) | 003, 004, 005 |
| GeneMark-ETP | zebrafish, 1.45 Gb | 12 h | same | ≤ 1,910 CPU-s/Mb | | same | 003, 004, 005 |
| GeneMark-ETP | mouse, 2.73 Gb | 18 h | same | ≤ 1,520 CPU-s/Mb | | same | 003, 004, 005 |
| GeneMark-EP+ | fly, 144 Mb, order-excluded proteins | about 5 h | 8 CPUs | ≤ 1,000 CPU-s/Mb | 8 GB | [10.1093/nargab/lqaa026](https://doi.org/10.1093/nargab/lqaa026) | 003 |
| GALBA | mammals, mean | 35 h 12 min | 48 threads | | | Tiberius paper Table 1, [10.1093/bioinformatics/btae685](https://doi.org/10.1093/bioinformatics/btae685) | 005 |
| AUGUSTUS | mammal, about 3 Gb | 2 h 25 min | 48 threads | ≤ 139 CPU-s/Mb | | Tiberius paper Table 1 | 005 |
| AUGUSTUS | 1.6 Mb fly test sequence | about 6 min | one 2.4 GHz PC, 2003 | 225 CPU-s/Mb | | [10.1093/bioinformatics/btg1080](https://doi.org/10.1093/bioinformatics/btg1080) | 003 |
| SNAP | per Mb | | 1 GHz machine, 2004 | 30 CPU-s/Mb | 100 MB per Mb | [10.1186/1471-2105-5-59](https://doi.org/10.1186/1471-2105-5-59) | 003, 004 |
| MAKER | *Schmidtea* | 4.1 h/Mb | single-core 2 GHz Mac, 2008 | 14,760 CPU-s/Mb, all stages | 2 GB | Cantarel et al. 2008, [10.1101/gr.6743907](https://doi.org/10.1101/gr.6743907) | 004 |
| Tiberius | human, cow, beluga, mean | 1 h 39 min | A100 80 GB + 48 CPU threads | 1.9 GPU-s/Mb | | [10.1093/bioinformatics/btae685](https://doi.org/10.1093/bioinformatics/btae685) | 002, 003, 004, 005 |
| Tiberius | training | 15 days | four A100s | | | same | 003, 004, 005 |
| Tiberius | 33-species panel, mean | 26 min | 72 threads + A100 80 GB; "80x faster than BRAKER3", not hardware-normalized | | | [10.64898/2026.04.24.720536](https://doi.org/10.64898/2026.04.24.720536) | 002, 003, 004, 005 |
| Tiberius | 12 species, mean (ANNEVO's benchmark) | 43.6 min | one RTX 4090, batch 8 | | 22.5 GB GPU | [ANNEVO README](https://github.com/xjtu-omics/ANNEVO) at 42c920f | 005 |
| ANNEVO | 12 species, mean | 12.2 min | one RTX 4090, batch 8 | | 3.8 GB GPU | same | 005 |
| ANNEVO | human | 19 min (82 min for the paper version) | one RTX 4090 | 0.37 GPU-s/Mb | | same | 005 |
| ANNEVO | 33-species panel, mean | 30 min | 72 threads + A100 | | | Tiberius clade preprint | 004 |
| Helixer | per mammal | 8 h 54 min | A100 | 10.7 GPU-s/Mb | | Tiberius paper Table 1 | 005 |
| Helixer | human; *Oryza brachyantha* | just under 8.5 h; 27 min | single-threaded pipeline + GPU | 9.9 GPU-s/Mb; about 6 GPU-s/Mb | | [10.1038/s41592-025-02939-1](https://doi.org/10.1038/s41592-025-02939-1) | 003 |
| Helixer | 12 species, mean | 286.6 min | one RTX 4090 | | 8.6 GB GPU | ANNEVO README | 005 |
| Helixer | 33-species panel, mean | 178 min | 72 threads + A100 | | | Tiberius clade preprint | 004 |
| CONTRAST | human model, training | about 12 h | 200 Xeon E5345 processors | | | [10.1186/gb-2007-8-12-r269](https://doi.org/10.1186/gb-2007-8-12-r269) | 003 |
| geneML | fungal genome | about 6 min | 8 CPU cores | | | [10.64898/2026.05.18.725946](https://doi.org/10.64898/2026.05.18.725946) | 002, 004 |
| EviAnn | mammal | under 1 h | one multicore server | | | [10.1038/s41592-026-03156-0](https://doi.org/10.1038/s41592-026-03156-0) (abstract) | 005 |
| minisplice | splice-site model | | | | 7,026 parameters | T-human-002 `[yang2026minisplice]` | 002 |

Failure rates: EGAPx on Galaxy, 45% of jobs (charter); in the BRAKER3
paper's comparison, FINDER completed on 7 of 11 species
([10.1101/gr.278090.123](https://doi.org/10.1101/gr.278090.123), via
T-human-005). No other source reports one.

Where the reviews disagree on a figure, the T-human-012 review is the
outlier each time and cites no source for its number: EGAPx "> 1,000 CPU
hours" per genome against the README's 71 and 425; BRAKER3 "50 to 200 CPU
hours" against 270 to 3,080 thread-hours implied by the paper; Helixer
"minutes on 1 GPU" against 8 h 54 min per mammal; AUGUSTUS "minutes per
locus" against 6 min per 1.6 Mb. This document uses the sourced figure in
each case, as the synthesis does (`docs/review/README.md` §2, which uses
nothing from T-human-012 as sole evidence).

Two consistency checks worth recording. The 2003 AUGUSTUS figure (225
CPU-s/Mb on a 2.4 GHz PC) and the 2024 upper bound (139 CPU-s/Mb across 48
threads) bracket the 165 and 169 CPU-s/Mb measured in §3 on *S. pombe* and
*N. crassa*, so a twenty-year-old GHMM runs at the same per-megabase cost it always did. And
ANNEVO's own benchmark and the Tiberius clade preprint, on different
hardware and species sets, agree on the ordering Tiberius/ANNEVO < Helixer
<< BRAKER3 by one to two orders of magnitude.

## 2. Documented costs, by tool

Each subsection quotes the primary source first (checked on 2026-09-10 by
fetching the README from `raw.githubusercontent.com`), then what the five
Phase 1 reviews (`relay/artifacts/T-human-00{2,3,4,5,12}/review.md`) and
the synthesis (`docs/review/README.md`, T-human-006, merged from PR #16)
add from the papers.

### 2.1 EGAPx

Primary source, [ncbi/egapx README](https://github.com/ncbi/egapx/blob/main/README.md):

- Requirements: "AWS Batch, SLURM/UGE cluster, or a r6a.4xlarge machine
  (32 CPUs, 256GB RAM)".
- "How long does EGAPx take to run?": on AWS Batch with a mix of
  r6i.2xlarge (8 CPU, 64 GB), r6i.4xlarge (16 CPU, 128 GB) and r6i.8xlarge
  (32 CPU, 256 GB) instances, *Drosophila melanogaster* (144 Mb, one
  short-read RNA-seq run of 48.4 M spots) "takes 71 CPU hrs and 3 wallclock
  hrs"; *Gallus gallus* (1.1 Gb, 10 short-read and 10 long-read runs)
  "takes 425 CPU hrs and 5.5 wallclock hrs".
- The bundled small *Dermatophagoides farinae* example reports
  "Duration: 1h 22m 12s, CPU hours: 6.3, Succeeded: 134" (134 Nextflow
  processes) and "usually runs in under an hour".
- Default resource labels go up to 128 GB per process, with advice to raise
  them to 200 GB for large genomes; a 40 Gb lungfish needed STAR with
  `--limitGenomeGenerateRAM 150000000000`.

Normalized: **1,775 CPU-s/Mb** (fly) and **1,391 CPU-s/Mb** (chicken), from
the README's own figures. Wall clock does not scale with genome size because
the pipeline fans out across instances: 7.6x the genome cost 6x the CPU-hours
and 1.8x the wall clock.

Galaxy production figures, from the charter (`relay/TASK.md`, citing
`nekrut/scalingPaper`, which this runner cannot read): 1,409 jobs over 290
genomes, about 416,000 CPU-hours, 45% failure rate, 24% CPU efficiency. That
is **295 CPU-h per job** and **1,434 CPU-h per genome** including failed
attempts; without the genome sizes in that set it cannot be normalized per
megabase, and the 24% efficiency means the *allocated* CPU-hours were about
four times the *used* ones. These are charter claims until the paper is
readable; they are not re-derived here.

### 2.2 Tiberius

Primary source, [Gaius-Augustus/Tiberius README](https://github.com/Gaius-Augustus/Tiberius/blob/main/README.md):

- "GPU with at least 8 GB memory recommended ... Tiberius can also run on
  CPU, but will be significantly slower".
- Recommended batch sizes: A100 (80 GB) 16; RTX 3090 (25 GB) 8; RTX 2070
  (8 GB) 2.
- Tiberius 2.0.0 (April 2026): "Runtime has been reduced by 30%."
- A web server now exists, so a user without a GPU can still run it; that
  moves the cost, it does not remove it.

Paper figures (as carried by the reviews and synthesis; see §1 table):
1 h 39 min for the human genome on an A100 80 GB with 48 threads. At
3,099 Mb that is **1.9 GPU-s/Mb** (pre-2.0.0; the 30% reduction would give
about 1.3). CPU-only timing is not published.

### 2.3 BRAKER3

Primary source, [Gaius-Augustus/BRAKER README](https://github.com/Gaius-Augustus/BRAKER/blob/master/README.md):

- Test runs (`test1.sh`, `test2.sh`, RNA-seq or proteins on the bundled
  small genome) have "expected runtime ~20 minutes" on an "Intel(R) Xeon(R)
  CPU E5530 @ 2.40GHz".
- "If you use more than 8 threads, this will not speed up all parallelized
  steps, in particular, the time consuming `optimize_augustus.pl` will not
  use more than 8 threads." So the wall clock of a BRAKER run has a floor
  set by AUGUSTUS training regardless of core count.
- `--UTR=on` "increases memory consumption of AUGUSTUS"; reducing cores
  reduces RAM.

Paper figures (reviews and synthesis): 5 h 37 min to 64 h 16 min on 48
threads across the 11 test genomes. Note the pipeline *trains* a species
model on every run, so its cost is training plus prediction; AUGUSTUS
prediction alone (§3) is a small part of it.

### 2.4 AUGUSTUS

No runtime is stated in the AUGUSTUS README or paper for whole genomes.
The measured rows in §3 are the figures for this document. A single
process on one core costs 57 to 169 CPU-s/Mb depending on the genome; the tool is trivially parallel by sequence, which is
how `benchmark/validation/` ran it (24 cores), and the cost per megabase
does not change with that parallelism, only the wall clock.

### 2.5 Helixer

Primary source, [weberlab-hhu/Helixer README](https://github.com/weberlab-hhu/Helixer/blob/main/README.md)
(redirects to `usadellab/Helixer`):

- "For realistically sized datasets, an Nvidia GPU or an Apple Silicon GPU
  (>= M1) ... will be necessary for acceptable performance."
- "The example below and all provided models should run on an Nvidia GPU
  with 11GB Memory (e.g. GTX 1080 Ti) and with 8 Gb (e.g. GTX 1080)."
- Demos: about 3 min (1-step inference) and 5 min (3-step inference) on a
  GPU.
- Install: "20-30 minutes" for an experienced user, "2-3 hours" otherwise.
- On a consumer card: the T-human-007 validation runs on an RTX 5080
  (16 GB) needed `--batch-size 8` because "the shipped default of 32
  exhausts 16 GB of GPU memory at `--subsequence-length 213840`" (the
  vertebrate model's window), and TensorFlow 2.15.1 ships no kernels for
  that card's compute capability, so every kernel is JIT-compiled on first
  use (`benchmark/validation/README.md`). Those runs record no time or
  memory in their declarations, so they give accuracy but no cost.

Paper figures (reviews and synthesis): 8 h 54 min per mammalian genome on
an A100, about **10.7 GPU-s/Mb** for a 3 Gb genome.

### 2.6 GeneMark-ETP

Primary source, [gatech-genemark/GeneMark-ETP README](https://github.com/gatech-genemark/GeneMark-ETP/blob/main/README.md):
the only resource statement is the usage line `gmetp.pl --cores 32 --cfg
input.cfg`. The licence is stated: "Creative Commons Attribution
NonCommercial ShareAlike 4.0", with third-party tools under their own
licences. That resolves the "unstated licence" item in the T-human-002
review: it is stated, and it is non-commercial.

Paper figures (reviews and synthesis): 3, 12 and 18 h on 64 cores for the
three genome-size classes reported in the GeneMark-ETP paper.

## 3. Measured costs

### 3.1 Runner

All measured rows come from one machine: the cloud container this agent
runs in. `Intel(R) Xeon(R) Processor @ 2.80GHz`, 4 cores, 15 GB RAM, no
GPU, Ubuntu 24.04, AUGUSTUS 3.5.0 from the Ubuntu package
(`augustus 3.5.0+dfsg-4build5`, `augustus-data` with 166 species parameter
sets). It is slower than a current laptop; treat the CPU-seconds as an
upper bound and the per-megabase ratios between tools as the durable
result. The *C. elegans* and *N. crassa* rows were run on two fresh
instances of the same runner class on 2026-09-10 (same CPU model, same
package versions; the `time` package had to be installed alongside
`augustus`).

### 3.2 Rows

| genome | size (Mb) | how run | wall clock (s) | CPU (user s) | **CPU-s / Mb** | peak RSS (MB) | genes | nucleotide F1 | locus F1 | transcript F1 | source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| *Schizosaccharomyces pombe* | 12.57 | 1 process, whole genome, 1 core | 2076 | 2073 | **165** | 403 | 4,456 | 0.95466 | 0.92316 | 0.70235 | [marx-0026](../relay/messages/20260910T032948Z-marx-0026.md) (2026-09-10, same runner class) |
| *Saccharomyces cerevisiae* | 12.16 | 17 processes, 4 in parallel | 209 | 698 | **57** | 233 | 5,154 | 0.95839 | 0.91850 | 0.78060 | this document, 2026-09-10 |
| *Plasmodium falciparum* | 23.29 | 14 processes, 4 in parallel | 622 | 2101 | **90** | 988 | 4,813 | 0.88362 | 0.87984 | 0.42754 | this document, 2026-09-10 |
| *Caenorhabditis elegans* | 100.29 | 7 processes, 4 in parallel | 3951 | 11409 | **114** | 730 | 14,999 | 0.86887 | 0.77123 | 0.30895 | this document, 2026-09-10 |
| *Neurospora crassa* | 41.10 | 21 processes, 4 in parallel | 1836 | 6939 | **169** | 381 | 8,475 | 0.94804 | 0.92276 | 0.61283 | this document, 2026-09-10 |

Machine-readable copy: [`docs/cost-baseline/measured.tsv`](cost-baseline/measured.tsv); declarations in the same directory follow `docs/benchmark.md` §3.3.

CPU-s/Mb is the sum of user CPU time over all processes of the run divided
by the genome length in megabases, so it is the same whether the sequences
were run serially or in parallel. Wall clock is reported separately. Peak
memory is the largest resident set of any single process.

### 3.3 What the rows show

- **The per-megabase cost of a GHMM is not a constant, and no single
  panel statistic predicts it.** The same tool costs 57, 90, 114, 165 and
  169 CPU-s/Mb on five genomes. Per megabase (`benchmark/panel.tsv`),
  *S. cerevisiae* has 25 introns, *P. falciparum* 363, *S. pombe* 418,
  *N. crassa* 435 and *C. elegans* 1,140; gene density is 499, 227, 408,
  238 and 199 genes per Mb; CDS fraction 73, 53, 57, 41 and 43 percent.
  The intron-poor yeast is the cheapest, as expected, but *S. pombe* costs
  1.4x more than *C. elegans* with 2.7x fewer introns per Mb, so intron
  density alone does not order the five. The *S. pombe* row was run as a
  single process over the whole genome, a protocol difference; *N. crassa*,
  run with the same per-sequence protocol as *C. elegans*, costs the same
  169 CPU-s/Mb as *S. pombe* with 2.6x fewer introns per Mb than
  *C. elegans*, so the protocol does not explain the gap either. The two
  filamentous-yeast-like ascomycetes at 48.5% and 36% GC agree with each
  other and disagree with the nematode. Within *C. elegans* the six
  nuclear chromosomes span 101 to 131 CPU-s/Mb, chrX the dearest; within
  *N. crassa* the seven chromosomes span only 166 to 171, so the
  within-genome spread is itself species-dependent. Whatever the new model costs per
  megabase will likewise vary by 2 to 3x across the panel, so the budget in
  §5 is stated as a ceiling over the panel, not a mean.
- **Memory grows with sequence length within a species and with
  composition across species.** The *C. elegans* chromosomes run from
  0.54 GB at 14.0 Mb (chrIII) to 0.73 GB at 21.2 Mb (chrV), roughly
  linear. Across species, one 3.3 Mb *P. falciparum* chromosome (19.5% GC)
  peaked at 0.97 GB, against 0.40 GB for the whole 12.6 Mb *S. pombe*
  genome in one process, 0.38 GB for the 10 Mb *N. crassa* chromosome 1
  (48.5% GC) and 0.23 GB for the largest *S. cerevisiae* chromosome. A memory budget quoted per megabase from one clade would be
  wrong by 2 to 3x on another.
- **Parallelism buys wall clock, not CPU.** *S. cerevisiae* took 209 s of
  wall clock on four cores for 698 CPU-s; *C. elegans* took 66 minutes of
  wall clock for 3.2 CPU-hours, and its 21 Mb chrV alone held one core for
  41 minutes, which is the wall-clock floor for a per-sequence split;
  *N. crassa* likewise took 31 minutes of wall clock for 1.9 CPU-hours
  with its 10 Mb chromosome 1 holding one core for 28 of them. The
  CPU-s/Mb figure is what to compare across machines, and it is the figure
  the budget uses.
- **Accuracy of the same run, for the record.** The *S. cerevisiae* scores
  equal `benchmark/validation/augustus-Saccharomyces_cerevisiae.json`
  (T-human-007's bioconda run) to five decimals in every field. *P.
  falciparum* with its own AUGUSTUS parameter set scores nucleotide F1
  0.884 and transcript F1 0.428 on the lowest-GC genome of the panel; that
  is a held-out species for our benchmark and a clade every current tool
  handles badly. *C. elegans*, the panel's first metazoan measured here
  and the nematode clade EGAPx declares out of scope (§2.1), scores
  nucleotide F1 0.869, locus F1 0.771 and transcript F1 0.309 with the
  stock `caenorhabditis` parameter set: 14,999 predicted loci against
  19,971 reference loci, with 2,874 reference loci fused into a
  neighbour's prediction and 6,486 missed outright
  (`docs/cost-baseline/augustus-Caenorhabditis_elegans.json` is not
  committed; the scorer output is reproduced by §3.4). Gene fusion in a
  199-genes-per-Mb genome with 65 bp median introns is the failure mode
  a length-aware decoder is meant to fix, and this is the number it has
  to beat at 114 CPU-s/Mb. *N. crassa* with its own `neurospora_crassa`
  set scores nucleotide F1 0.948, locus F1 0.923 and transcript F1 0.613,
  8,475 predicted genes against 9,757 in the reference: between the two
  yeasts on every metric, as its intron count per gene is.
- **A convention trap, per species.** The Ubuntu package's
  `schizosaccharomyces_pombe` and `saccharomyces_cerevisiae_S288C` configs
  set `stopCodonExcludedFromCDS true` (score with `--stop-outside-cds`);
  `pfalciparum` and `caenorhabditis` set it `false` (score without). The
  package ships two *Neurospora* sets that disagree with each other:
  `neurospora_crassa` says `true`, `neurospora` says `false`. The
  scorer detects both and warns, so the wrong flag is loud, not silent,
  but every new species run needs the check.

### 3.4 Reproducing

```
apt-get install augustus augustus-data time
python3 benchmark/fetch.py --species Saccharomyces_cerevisiae --what fasta --dest /tmp/panel
python3 benchmark/fetch.py --species Saccharomyces_cerevisiae --what gff   --dest /tmp/panel
# split by sequence, one AUGUSTUS process per sequence under /usr/bin/time -v,
# four at a time; ids made unique before concatenation (benchmark/validation/README.md)
docs/cost-baseline/run_augustus.sh Saccharomyces_cerevisiae saccharomyces_cerevisiae_S288C 4
python3 benchmark/score.py --reference /tmp/panel/Saccharomyces_cerevisiae/*_genomic.gff.gz \
    --prediction augustus.gff3 --species Saccharomyces_cerevisiae \
    --declaration docs/cost-baseline/augustus-Saccharomyces_cerevisiae.yaml \
    --genome /tmp/panel/Saccharomyces_cerevisiae/*_genomic.fna.gz --stop-outside-cds --out sc.json
# the metazoan row: 66 min wall on 4 cores, 3.2 CPU-h, 0.73 GB peak; no --stop-outside-cds
python3 benchmark/fetch.py --species Caenorhabditis_elegans --what fasta,gff --dest /tmp/panel
docs/cost-baseline/run_augustus.sh Caenorhabditis_elegans caenorhabditis 4 /tmp/panel /tmp/run-cel
python3 benchmark/score.py --reference /tmp/panel/Caenorhabditis_elegans/*_genomic.gff.gz \
    --prediction /tmp/run-cel/augustus.gff3 --species Caenorhabditis_elegans \
    --declaration docs/cost-baseline/augustus-Caenorhabditis_elegans.yaml \
    --genome /tmp/panel/Caenorhabditis_elegans/*_genomic.fna.gz --out cel.json
# the filamentous fungus: 31 min wall on 4 cores, 1.9 CPU-h, 0.38 GB peak; --stop-outside-cds
python3 benchmark/fetch.py --species Neurospora_crassa --what fasta,gff --dest /tmp/panel
docs/cost-baseline/run_augustus.sh Neurospora_crassa neurospora_crassa 4 /tmp/panel /tmp/run-nc
python3 benchmark/score.py --reference /tmp/panel/Neurospora_crassa/*_genomic.gff.gz \
    --prediction /tmp/run-nc/augustus.gff3 --species Neurospora_crassa \
    --declaration docs/cost-baseline/augustus-Neurospora_crassa.yaml \
    --genome /tmp/panel/Neurospora_crassa/*_genomic.fna.gz --stop-outside-cds --out nc.json
```

The Ubuntu package's yeast configs set `stopCodonExcludedFromCDS true`, so
`--stop-outside-cds` is required for them; `pfalciparum` sets it `false`,
so the flag must be dropped, as does `caenorhabditis`. The scorer warns either way if the flag is wrong.

The genomic coding fractions quoted in §5.3 are reproduced by

```
python3 docs/cost-baseline/cds_union.py Saccharomyces_cerevisiae Homo_sapiens Mus_musculus Zea_mays
```

which downloads each species' RefSeq GFF named in `benchmark/panel.tsv`,
checks its md5 against the panel, and writes the union of CDS intervals
beside the per-isoform sum (`cds_union.tsv`; the per-isoform column
reproduces the panel's `cds_fraction_pct` exactly). Standard library only;
the four GFFs are about 200 MB and are cached in `$CDS_UNION_CACHE`.

## 4. Our floor: the KA/KS window test

T-human-010 (`baselines/kaks/`) ran the pairwise KA/KS classifier on 32
genes across four reference/informant sets in 365 s wall clock, "of which
the fetch was nearly all" (`baselines/kaks/README.md`,
`results/manifest.json`: `seconds: 365.1`). The statistic itself is a
Nei-Gojobori count over a window of a pairwise alignment and costs
milliseconds per window in pure Python. The comparative signal is therefore
not where the cost of a comparative method lives; the alignment is (see
§5.3).

## 5. Target budget

The charter asks for accuracy at or above the best current tools "at a
small fraction of their compute". The rows above put numbers on the
denominator. This section turns them into a ceiling the design in
T-human-011 must fit, and shows one consequence the numbers force.

### 5.1 The reference points

| quantity | value | source |
|---|---|---|
| Best ab initio GPU tool, per Mb | Tiberius, about 1.9 GPU-s/Mb on an A100 80 GB (1.3 after the 2.0.0 speed-up) | §2.2 |
| Second GPU tool, per Mb | Helixer, about 10.7 GPU-s/Mb on an A100 | §2.5 |
| Classical GHMM on one core, per Mb | AUGUSTUS, 165 to 169 CPU-s/Mb (*S. pombe*, *N. crassa*), 114 on the 100 Mb *C. elegans*; see §3 | §3 |
| Evidence pipeline, per Mb | EGAPx, 1,400 to 1,800 CPU-s/Mb plus a 32-CPU, 256 GB machine | §2.1 |
| Consumer GPU, peak | RTX 4090: 83 TFLOPS FP32 shader, 24 GB ([NVIDIA product page](https://www.nvidia.com/en-us/geforce/graphics-cards/40-series/rtx-4090/), read 2026-09-10); an RTX 2070 / GTX 1080 class card with 8 GB is the floor Tiberius and Helixer both state | §2.2, §2.5 |
| One CPU core, sustained | assumed 3 × 10^10 FLOP/s (a 2.8 GHz Xeon core with two 256-bit fused multiply-add ports peaks at 32 FLOP/cycle × 2.8 GHz = 90 GFLOP/s, 45 with one port; small-matrix inference reaches well under either; this is an assumption, not a measurement, and §5.3 conclusion 1 states what survives if it is 3x too low) | this document |
| Consumer GPU versus the A100 the GPU tools were timed on, tensor throughput | RTX 4090 about 165 TFLOPS BF16 dense, FP32 accumulate ([NVIDIA Ada GPU architecture whitepaper](https://images.nvidia.com/aem-dam/Solutions/geforce/ada/nvidia-ada-gpu-architecture.pdf), RTX 4090 specification table, read 2026-09-10); A100 312 TFLOPS BF16 dense ([NVIDIA A100 datasheet](https://www.nvidia.com/en-us/data-center/a100/), read 2026-09-10). On the number that governs inference the target card is about half the A100, so a per-Mb ceiling on the 4090 is a stricter demand on the algorithm than the same number on an A100 | this document |
| Precedent for the parameter budget | HyphAeon, about 2 M parameters | `relay/TASK.md` background |

### 5.2 The budget

| resource | ceiling | why |
|---|---|---|
| Parameters | 5 M, target 2 to 3 M | the task's stated target; 2 M is the HyphAeon precedent and the number the FLOP arithmetic below uses |
| One CPU core, inference | **≤ 15 CPU-s/Mb on the §3.1 runner** (Xeon 2.80 GHz, one core), equivalently **≤ 1/10 of AUGUSTUS's CPU-s/Mb on the same machine and genome** | the absolute is bound to the machine that produced §3 so it can be checked against `measured.tsv`; the ratio is the portable form, because a faster core lowers our number and AUGUSTUS's together and an absolute with no machine cannot be falsified. An order of magnitude under AUGUSTUS (165 to 169 on *S. pombe* and *N. crassa*, same runner) and two under EGAPx (1,400 to 1,800); a 3.1 Gb human genome in 13 core-hours, a 12 Mb yeast in 3 minutes, on a laptop with no GPU |
| One consumer GPU (8 to 24 GB), inference | **≤ 0.5 GPU-s/Mb on an RTX 4090 class card** | 3.8x under Tiberius's 1.9 GPU-s/Mb on an A100 (2.6x under the post-2.0.0 1.3), measured on a card with about half the A100's tensor throughput (165 versus 312 TFLOPS BF16, §5.1), so the demand is about 5 to 8x less work per megabase than Tiberius, not a quarter; the human genome in 26 minutes, the whole 20-genome panel (13.96 Gb, the sum of `genome_bp` in `benchmark/panel.tsv`) in 1.9 hours |
| Peak host memory | ≤ 8 GB | the laptop constraint of the charter; AUGUSTUS runs in 0.4 GB, so anything that needs a whole-genome tensor resident has already lost |
| Peak GPU memory | ≤ 8 GB | the floor both GPU tools state; batch size is the knob, as their READMEs say |
| Failure rate | 0 per genome on the panel, by construction: no per-genome training, no external pipeline steps that can time out | the EGAPx figure is 45% on Galaxy (§2.1); every training-in-the-loop pipeline inherits the same fragility |
| Preprocessing, alignment excluded | counted in the CPU budget above | the alignment itself is accounted separately (§5.3) because it is shared with every other comparative user of the genome |

### 5.3 What the arithmetic forces: the model cannot look at every base of every taxon

A forward pass of a transformer-style model costs about 2 × parameters
FLOP per token. With 2 M parameters, one token per base and 20 aligned taxa
(the shape of a whole-genome alignment fed in densely) that is
2 × 2 × 10^6 × 10^6 × 20 = 8 × 10^13 FLOP per megabase:

| shape | FLOP/Mb | one core at 3 × 10^10 FLOP/s | RTX 4090 at a realistic 10^13 FLOP/s (12% of FP32 shader peak, 6% of BF16 tensor peak) |
|---|---|---|---|
| per base, 20 taxa | 8.0 × 10^13 | 2,700 s/Mb | 8 s/Mb |
| per codon, 20 taxa | 2.7 × 10^13 | 890 s/Mb | 2.7 s/Mb |
| per base, target only | 4.0 × 10^12 | 130 s/Mb | 0.4 s/Mb |
| per base, 20 taxa, on 5% of the genome | 4.0 × 10^12 | 130 s/Mb | 0.4 s/Mb |
| per codon, 20 taxa, on 5% of the genome | 1.3 × 10^12 | 44 s/Mb | 0.13 s/Mb |
| per codon, 20 taxa, on 2% of the genome, attention excluded | 5.3 × 10^11 | 18 s/Mb | 0.05 s/Mb |

The attention term is left out; with windows of a few thousand tokens it is
of the same order as the parameter term, so the real numbers are about
twice these. Three things follow, and they are constraints on
T-human-011, not preferences:

1. **A dense per-base model over the alignment is slower than AUGUSTUS on
   a CPU by more than an order of magnitude** (2,700 versus 169 s/Mb) even
   at 2 M parameters, and only just meets Tiberius on a GPU. The CPU
   budget cannot be met that way at all. This does not rest on the assumed
   core speed: at the full 9 × 10^10 FLOP/s peak of §5.1 the same shape
   costs 890 s/Mb, still five times AUGUSTUS and 60 times the ceiling.
2. **The model must run on candidate support, not on the genome.** The
   genomic coding fraction, the union of CDS intervals over the genome, is
   1.3% for human, 1.4% for mouse, 2.0% for maize and 72.5% for
   *S. cerevisiae* (`docs/cost-baseline/cds_union.tsv`, produced by
   `cds_union.py` from the panel's RefSeq GFFs, md5-checked against
   `benchmark/panel.tsv`). The panel's `cds_fraction_pct` column sums CDS
   over every isoform, so it reads 9.6% for human, 7.2x the union, and is
   not the number to size a candidate stage by. A candidate-region stage
   that keeps 2 to 5% of a mammalian genome, which is 1.5 to 4x the coding
   bases with room for flanks, and nearly all of a yeast genome, and a
   codon-level (not base-level) token axis, are what bring the per-Mb cost
   under the ceiling. The candidate stage
   itself must therefore be a cheap sequence scan, of the cost class of the
   KA/KS window test (§4) or a GHMM emission pass, not a neural network.
3. **The tree costs nothing at inference if it is coordinates.** Tree-RoPE
   and MDS embeddings are fixed per alignment, computed once per genome
   (T-human-008's window fetcher already emits the tree with every window),
   so injecting the phylogeny as a metric adds no FLOP per token; adding it
   as extra tokens (the (e) arm of the ablation in `docs/review/candidates.md`)
   does. The budget therefore favours the metric form independently of
   whether it is more accurate.

Alignment cost is the elephant this table leaves out: it is paid once per
genome, is shared with every other comparative analysis, and for a genome
absent from every public alignment it is the whole cost. That is
T-human-008's finding and T-human-011's risk, not a budget line here.

### 5.4 Cluster request for the two rows this document cannot produce

EGAPx end to end and Tiberius on the panel are the two documented rows
that most need a measured counterpart, and neither fits a laptop. The
`alert` to `human` that the task asks for is
[marx-0033](../relay/messages/20260910T115307Z-marx-0033.md), posted
2026-09-10; it requests, for one attempt each:

- EGAPx on *Ciona intestinalis* (140 Mb, `heldout`, no RNA-seq required):
  one 32-CPU, 256 GB node; by the README's fly figure, about 70 to 100
  CPU-hours and 3 to 4 hours wall clock; 200 GB scratch.
- Tiberius 2.0.7 `vertebrates`, `fungi` and `insecta` models on the eleven
  panel species those clades cover: *H. sapiens* 3.10 Gb, *M. musculus*
  2.73, *X. tropicalis* 1.45, *D. rerio* 1.45, *G. gallus* 1.05,
  *T. rubripes* 0.38, *A. mellifera* 0.23, *D. melanogaster* 0.14,
  *N. crassa* 0.04, *S. pombe* 0.01, *S. cerevisiae* 0.01; in total
  10.60 Gb (`genome_bp` in `benchmark/panel.tsv`). One GPU with ≥ 8 GB (an
  A100 or the RTX 5080 box used by `benchmark/validation/`); at the paper's
  1.9 GPU-s/Mb about 5.6 GPU-hours, at the post-2.0.0 1.3 GPU-s/Mb about
  3.8, plus first-run kernel compilation. The five held-out species among
  them (*H. sapiens*, *G. gallus*, *T. rubripes*, *A. mellifera*,
  *S. pombe*) are 4.78 Gb, 2.5 GPU-hours at 1.9. The alert as posted said
  5.5 Gb and about 3 GPU-hours, about half the run; messages are never
  edited, so the corrected figure is restated in marx's note of
  2026-09-10 that accompanies this revision (lenin's review, fix 2).

## 6. What is not here yet

- Tiberius, Helixer and BRAKER3 measured on the benchmark panel: needs a GPU
  for the first two and a day of CPU for the third. The RTX 5080 machine of
  `benchmark/validation/` is the natural place; those runs exist but their
  declarations record no time or memory.
- EGAPx end to end: below this runner's and a laptop's minimum. The
  `alert` with the cluster estimate is posted (§5.4); the run waits on a
  `decision`.
- Failure rate per genome for anything other than EGAPx on Galaxy: no
  source reports one.
- AUGUSTUS on a gigabase genome. Every measured row is under 100 Mb, so
  the 3x per-Mb spread of §3.3 is a spread over small genomes, and the
  budget's per-megabase form assumes cost is flat in genome size.
  *Z. mays* (2.18 Gb) at 57 to 169 CPU-s/Mb is 35 to 100 CPU-hours on the
  §3.1 runner, which does not fit the 35-minute foreground limit of one
  tick on 4 cores; not measured for that reason. One 100 Mb chromosome of
  a gigabase genome would fit and is the next measured row to add.
- Wall clock and memory for the Tiberius and Helixer runs already in
  `benchmark/validation/`. Re-running them under `/usr/bin/time -v` on the
  RTX 5080 box would replace two "documented only" cells without a cluster
  `decision`, at a fraction of the §5.4 request. This runner has no GPU,
  so that is a request to whoever operates that machine.
