# S. pombe normalization row (proposal 6.1 / cost-baseline 5.2), 2026-09-21

Candidate A (frozen smoke checkpoint `pooledrun/best.pt`, 20 steps on yeast
chr I windows; code `8440560`, PR 38 at `cb2b27d`, clean tree) and AUGUSTUS
3.5.0 on the three *S. pombe* nuclear chromosomes (12,571,820 bases; the
19.4 kb mitochondrion excluded), **same machine, one process each pinned
to core 2**: Intel Core Ultra 9 285K, 62 GB, no GPU, Python 3.11.14,
torch 2.14.0+cpu, `OMP_NUM_THREADS=1`. `benchmark/leakage_check.py` ran
first (`leakage_check.out`: 0 violations; *S. pombe* is a cross-clade
held-out species, nearest train species at phylum rank). S. pombe was not
scored and no checkpoint was selected on it; this is the runtime
normalization only. The AUGUSTUS command is marx-0026's
(`--species=schizosaccharomyces_pombe --gff3=on --UTR=off`, bioconda
3.5.0 in a micromamba env, `augustus_install.log`). The GFF3 outputs
(A: 14.5 MB per chromosome; AUGUSTUS: 6.1 MB) are not committed;
SHA-256s are in `gff3_sha256.txt`, `A-f32/gff3_sha256.txt` and
`augustus/outputs_sha256.txt`. Each A row aggregates three processes
(one per nuclear chromosome, both strands inside); the excluded
mitochondrion is NC_088682.1 (19,433 bases; engels-0093).

| run | genome Mb | stage CPU-s (sum) | process user s | system s | process wall s | **stage / Mb** | user / Mb | user + system / Mb | peak RSS (max over chromosomes) | chains / genes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A, 1 segment per strand (exact strand decode), margin 491 | 12.5718 | 507.24 | 500.61 | 8.46 | 509.30 | 40.35 | 39.82 | 40.49 | 1.49 GiB (chr I) | 121,937 |
| A, 19 segments per strand, seam overlap 4,096, margin 491 | 12.5718 | 115.96 | 103.47 | 14.32 | 117.89 | **9.22** | 8.23 | 9.37 | 1.31 GiB (chr I) | 121,937 |
| AUGUSTUS 3.5.0, one process, whole nuclear genome | 12.5718 | – | 646.33 | 0.17 | 646.75 | – | **51.41** | 51.43 | 0.40 GiB | 4,452 |

Stages at 19 segments per genome Mb: preprocess 0.11, encoder 3.74,
decode 5.31, output 0.06, I/O 0.01. Per chromosome (stage / user /
user+system per Mb): chr I 9.17 / 8.17 / 9.28, chr II 9.21 / 8.21 /
9.34, chr III 9.38 / 8.42 / 9.61; exact decode 40.15–40.56 / 39.58–40.02
/ 40.26–40.70. Encoded bases are 1.080–1.083× the oriented bases
(`measure_*.json`).

**Normalization.** AUGUSTUS costs 51.4 user CPU-s/Mb here against 164.9
on the cost-baseline runner (marx-0026, Xeon 2.80 GHz): this machine is
3.21× faster on the same command. A at 19 segments is **1/5.5 of
AUGUSTUS** on user+system (9.37 / 51.43), 1/6.2 on user only, 1/5.6 on
the stage sum; the portable target of cost-baseline 5.2 is **1/11**,
i.e. ≤ 4.67 CPU-s/Mb on this machine. Machine-normalized to the baseline
runner, 9.37 × 3.21 = 30.0 CPU-s/Mb against the 15 ceiling. **A misses
the CPU target on the normalization row by 2.0× (1.8× on user only).**
Memory is inside 8 GB. The exact strand decode is 1/1.3 of AUGUSTUS.

**Segment-count probe** (`A-probe/`, chr III, 38 and 76 segments per
strand): stage 8.72 and 8.96 CPU-s/Mb (decode 11.30 and 11.00 s against
13.26 at 19 segments; the encoder rises 9.30 → 9.66 → 10.53 s with the
extra margin bases), RSS 1.87 and 3.29 GiB (1.14 at 19). The decode
stage is at its operand floor near 4.5 CPU-s per genome Mb; more
segments no longer buy time and cost memory linearly in B. Segment
decodes differ from the exact strand decode by 3 chains of 121,937 at
19 segments (chr II 2, chr III 1), 1 at 38 and 4 at 76 on chr III
(`chaindiff.out`), consistent with the seam contract.

**Float32 decode** (`A-f32/`, revision step 4, code `aee0134`: `measure
--dtype float32`; emissions, motif bias, duration tables and the scan in
float32, with the carried scores rebased at every tile seam in the segment
mode so the running magnitude is a tile's, not a segment's; same
checkpoint, protocol and core, 2026-09-21 10:11–10:22 UTC):

| run | stage CPU-s (sum) | process user s | system s | stage wall s | process wall s | **stage / Mb** | user / Mb | user + system / Mb | peak RSS (max) | chains |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A, 1 segment per strand, float32 | 499.30 | 493.04 | 8.09 | 499.52 | 501.37 | 39.72 | 39.22 | 39.86 | 1.49 GiB (chr I) | 121,935 |
| A, 19 segments per strand, float32 | 103.99 | 99.87 | 6.06 | 104.04 | 105.98 | **8.27** | 7.94 | 8.43 | 0.97 GiB (chr I) | 121,938 |

Stages at 19 segments per genome Mb: preprocess 0.11, encoder 3.76,
decode **4.34** (float64: 5.31), output 0.05, I/O 0.01; per chromosome
stage / Mb 8.24 / 8.16 / 8.54. The exact strand decode barely moves
(decode 36.4 → 35.8 CPU-s per Mb). Float32 takes 10 % off the 19-segment
row (9.22 → 8.27 stage, 9.37 → 8.43 user + system: 1/6.1 of AUGUSTUS,
**miss 1.8× on user + system, 1.7× on user**) and 26 % off its RSS; it does
not halve the decode as proposed, because the decode is not
bandwidth-bound: a tile scan of 12,288 steps at B = 38 costs 0.93 s in
float64 and 0.73 in float32, of which the per-step loop is 0.52 / 0.47 s
(42 / 38 µs per step, ~15 torch ops on (38, 24) tensors — dispatch, not
arithmetic) and the per-tile operand build 0.41 / 0.26 s. **Outputs:**
float32 differs from float64 by 103 chains of 121,937 at 19 segments
(47 / 40 / 16 per chromosome; 93 are internal exon-boundary shifts in
short 2–5-exon chains, 7 change the exon count, 3 have no partner) and
by 112 at 1 segment, and the float32 exact and segment decodes differ
from each other by ~120 chains where the float64 pair differ by 3
(`chaindiff_f32.out`; classification by overlap-matching the
differing chains). These are near-tie decisions of the flat smoke
checkpoint at float32 resolution; the agreement must be re-measured with
a fitted checkpoint before float32 is adopted for a reported row. The
float64 rows above remain the rows of record.

Compute: ~0.37 CPU-h local for the float64 rows (nine process records on
user + system, 0.3661 h: A 0.17 h, AUGUSTUS 0.18 h, probe 0.01 h), plus
0.17 CPU-h for the six float32 runs; cluster CPU-hours 0, GPU-hours 0.

## Fitted checkpoint row (`A-v2/`, 2026-09-22T08:08Z–08:10Z)

Same machine, core, command and float64 dtype as the row of record above
(`run_A_v2.sh`; code `98c832c`, clean tree, source digest `56e472d1…`),
with `fit-cpu-v2/best.pt` (step 900, sha256 `99f773ed…2a8cf6`) in place
of the smoke checkpoint; 19 segments per strand only (the exact strand
decode was the smoke row's control and cost is checkpoint-independent, as
chr I / chr V showed). S. pombe stays unscored and unselected-on; the
leakage check on file (`../leakage_check.out`, 0 violations) covers the
same inputs. GFF3s are not committed (SHA-256 in `A-v2/gff3_sha256.txt`).

| run | genome Mb | stage CPU-s (sum) | process user s | system s | process wall s | **stage / Mb** | user / Mb | user + system / Mb | peak RSS (max over chromosomes) | chains |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A, fit-cpu-v2 step 900, 19 segments per strand, float64 | 12.5718 | 112.14 | 99.93 | 14.00 | 114.01 | **8.92** | 7.95 | 9.06 | 1.31 GiB (chr I) | 9,555 |
| A, smoke checkpoint, same settings (row of record above) | 12.5718 | 115.96 | 103.47 | 14.32 | 117.89 | 9.22 | 8.23 | 9.37 | 1.31 GiB (chr I) | 121,937 |
| AUGUSTUS 3.5.0 (above) | 12.5718 | – | 646.33 | 0.17 | 646.75 | – | 51.41 | 51.43 | 0.40 GiB | 4,452 genes |

Per chromosome (stage / user / user+system per Mb): chr I 8.91 / 7.92 /
9.02, chr II 8.88 / 7.90 / 9.01, chr III 9.02 / 8.11 / 9.26. Stages per
genome Mb: preprocess 0.11, encoder 3.73, decode 5.07, output 0.00, I/O
0.01. The fitted checkpoint is **1/5.7 of AUGUSTUS** on user+system
(9.06 / 51.43), 1/6.5 on user only, 1/5.8 on the stage sum; machine-
normalized 9.06 × 3.21 = 29.1 CPU-s/Mb against 15. The row moves by 3 %
against the smoke row (traceback of 9,555 rather than 121,937 chains),
so the normalization verdict is unchanged: **A misses the portable
1/11 target by 1.9× on the fitted checkpoint** (1.7× on user only);
memory 1.31 GiB. The chain count (9,555 against 4,452 AUGUSTUS genes and
the reference's ~5,100 protein-coding genes) is reported for the cost
row only, as output volume: S. pombe is unscored here, so chain lengths,
unmatched loci and fragmentation are not measured, and reading the count
as the short-chain over-prediction seen on chr V is a hypothesis, not an
accuracy measurement on S. pombe.

Compute for this row: 0.032 CPU-h (113.9 user + system s); cluster
CPU-hours 0, GPU-hours 0.
