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
SHA-256s are in `gff3_sha256.txt` and `augustus/outputs_sha256.txt`.

| run | genome Mb | stage CPU-s (sum) | process user s | system s | wall s | **stage / Mb** | user / Mb | user + system / Mb | peak RSS (max over chromosomes) | chains / genes |
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

Compute: ~0.3 CPU-h local (A 0.17 h, AUGUSTUS 0.18 h, probe 0.01 h);
cluster CPU-hours 0, GPU-hours 0.
