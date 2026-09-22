# Metazoan development-chromosome row: *C. elegans* chromosome V, 2026-09-21

Candidate A (frozen smoke checkpoint `pooledrun/best.pt`, 20 steps on yeast
chr I windows; code `aee0134`, PR 38 at `d78114c`, clean tree) and AUGUSTUS
3.5.0 on *C. elegans* chromosome V (`NC_003283.11`, 20,924,180 bases, the
largest chromosome of a train species, chosen to match the 20 Mb projection
of a-pilot 3.2), **same machine**: Intel Core Ultra 9 285K, 62 GB, no GPU,
Python 3.11.14, torch 2.14.0+cpu, `OMP_NUM_THREADS=1`. The A runs were
pinned to core 2 and ran one after another; AUGUSTUS was pinned to core 4
and ran concurrently with the first two A runs (one single-threaded process
per core, no shared core; see `A/run.log` and `augustus/run.log` for the
timestamps). Sources are the pinned WBcel235 FASTA/GFF (`sources_md5.txt`,
both MD5s equal to `model/labels/manifests/Caenorhabditis_elegans.summary.json`);
`benchmark/leakage_check.py` ran first (`leakage_check.out`: 0 violations).
*C. elegans* is a train species; chromosome V is declared the development
chromosome in `config.json` (`dev_seqids`) for this row, a runtime
measurement only — nothing was scored and no checkpoint was selected on it.
Seam overlap is **16,384** (an 8,192-base chain-containment guarantee:
97% of chr V's 4,965 admitted representatives span ≤ 8,192 bases, p99
13.7 kb, max 81 kb) instead of the 4,096 used on the yeasts, since a-pilot
3.2 noted the overlap must scale with metazoan gene lengths; at 19 segments
per strand that is 1.014× oversampling (oriented 42,438,184 bases against
2 × 20,924,180) and 1.096× encoded bases with the 491-base margin. The
AUGUSTUS command is marx-0026's with `--species=caenorhabditis`. GFF3
outputs (A: 61.8 MB; AUGUSTUS: 5.5 MB) are not committed; SHA-256s are in
`gff3_sha256.txt`.

| run | genome Mb | stage CPU-s (sum) | process user s | system s | stage wall s | process wall s | **stage / Mb** | user / Mb | user + system / Mb | peak RSS | chains / genes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A, 1 segment per strand (exact strand decode), float64, margin 491 | 20.9242 | 884.38 | 878.39 | 6.70 | 884.74 | 885.46 | 42.27 | 41.98 | 42.30 | 4.87 GiB | 208,796 |
| A, 19 segments per strand, seam overlap 16,384, float64, margin 491 | 20.9242 | 201.67 | 177.43 | 24.98 | 201.70 | 202.44 | **9.64** | 8.48 | 9.67 | 2.04 GiB | 208,796 |
| A, 19 segments per strand, seam overlap 16,384, float32 (seam rebase), margin 491 | 20.9242 | 175.07 | 164.17 | 11.60 | 175.16 | 175.86 | 8.37 | 7.85 | 8.40 | 1.74 GiB | 208,801 |
| AUGUSTUS 3.5.0, one process, chr V | 20.9242 | – | 1092.12 | 0.34 | – | 1092.48 | – | **52.19** | 52.21 | 0.70 GiB | 3,357 |

Stages per genome Mb at 19 segments, float64: preprocess 0.12, encoder
3.91, decode 5.50, output 0.09, I/O 0.02; float32: 0.11 / 3.74 / 4.41 /
0.08 / 0.02; exact decode: 0.13 / 3.94 / 38.07 / 0.10 / 0.02. The
19-segment float64 row is the row of record, as on *S. pombe*.

**Against the budget.** AUGUSTUS costs 52.2 user CPU-s/Mb on chr V here,
within 2% of its 51.4 on the *S. pombe* nuclear genome on this machine,
so the same-command machine factor of the *S. pombe* row (3.21× faster than
the cost-baseline runner) is used. A at 19 segments is **1/5.4 of AUGUSTUS
on user + system** (9.67 / 52.21), 1/6.2 on user only (8.48 / 52.19), 1/5.4
on the stage sum; the portable target of cost-baseline 5.2 is **1/11**,
i.e. ≤ 4.75 CPU-s/Mb on this machine. Machine-normalized, 9.67 × 3.21 =
31.0 CPU-s/Mb against the 15 ceiling. **The metazoan row misses the CPU
target by 2.0× (1.8× on user only), the same miss as the *S. pombe* row
(2.0× / 1.8×); the float32 option gives 8.40 (1/6.2, miss 1.8×).** Memory:
2.04 GiB at 19 segments and 1.74 in float32 (a-pilot's source-derived
projection for the streamed path was ~1.5 GB of live tensors at 20 Mb /
19 segments, before the runtime floor of ~0.35 GiB and the longer overlap);
the exact strand decode reaches 4.87 GiB, because its one-row traceback
expansion and packed store scale with the 20.9 Mb row — inside 8 GB but
not the configuration of record. The per-genome-Mb cost is
chromosome-length independent within 5% between 12.6 Mb of *S. pombe*
(9.22 / 8.23 / 9.37) and 20.9 Mb of *C. elegans* (9.64 / 8.48 / 9.67); the
difference is the 1.4% seam oversampling of the longer overlap, the
margin on the larger tile count and a higher system time (25.0 s against
14.3 s over 12.6 Mb, allocation churn of the segment mode; float32 halves it).

**Outputs.** The 19-segment float64 decode is **byte-identical to the
exact strand decode** (same SHA-256, `chaindiff.out`: 208,796 chains, 0
differ): with the 8,192-base containment guarantee no chain of this
checkpoint on chr V is cut at a segment seam (on the yeasts at overlap
4,096, 3–4 chains per chromosome were). Float32 differs from float64 by
194 / 199 chains of 208,796 (0.09%; on *S. pombe* 103 of 121,937, 0.08%) —
near-ties of the flat smoke checkpoint, to be re-checked with a fitted
checkpoint before float32 becomes the row of record. The chain counts
(208,796 against 4,965 admitted representatives; AUGUSTUS 3,357 genes)
are those of a 20-step smoke checkpoint and carry no accuracy claim.

**Compute.** Local CPU: A 0.35 CPU-h (three runs), AUGUSTUS 0.30 CPU-h;
cluster CPU-hours 0, GPU-hours 0. No held-out species touched.

Files: `run_A.sh`, `config.json`, `leakage_check.out`, `sources_md5.txt`,
`A/` and `A-f32/` (`measure_*.json`, stdout, `/usr/bin/time -v`, `run.log`),
`augustus/` (run script, `run.log`, `/usr/bin/time -v`, stderr),
`gff3_sha256.txt`, `chaindiff.out` (produced with
`../pombe-normalization/chaindiff.py`).
