# score-canonical/: the fit-cpu-v3 checkpoint decoded with the opt-in restricted-support splice mask, 2026-09-22T16:25Z

The section 3.4 **decoder** ablation. It restricts the support of the
accepted grammar rather than describing what biology permits: proposal
3.1 gives motifs finite scores that never prohibit a junction, and the
chr V reference itself holds 39 introns outside GT-AG/GC-AG, one of
which unmasked v3 matched and the masked decode cannot reach. The
unmasked decode in `../score-final/` stays the default and the baseline
(stalin-0114). Same weights (`best.pt` sha256
`8a32c93e…`, step 3,000 of 3,000, dev NLL 28.640), same tiling, same
scorer, same core; the only difference from `../score-final/` is
`--canonical-splice`, which makes a donor emission `-inf` wherever the
first two intron bases are a concrete dinucleotide other than `GT`/`GC`
and an acceptor `-inf` wherever the last two are not `AG`. Ambiguous
bases and window edges are **not** masked, so the mask can only remove a
site on positive evidence. No refitting: the checkpoint is the one that
produced the rows in `../score-final/`.

`run_canonical.sh` is `../score-final/run_final.sh` with the flag added
and the output directory renamed; substituting both back and diffing
against `run_final.sh` is empty, so the gate, the failure propagation and
every other workload argument are unchanged. All four workloads exit 0
(`run_canonical.out`), the gate passed, and `summarize.py` (copied
unchanged) reproduces every number below from the committed JSON and
GFF3. Measurement commit `77c04bc`; the fit ran at `6c7aac1` and
`../score-final/` was measured at `85d534c`. The only code difference
between `85d534c` and `77c04bc` inside `model/` is the mask itself, which
is off unless the flag is given (`git diff 85d534c 77c04bc -- model/`).

## Cost (one core, `taskset -c 2`, `/usr/bin/time -v`, float64, 19 segments)

| chromosome | CPU-s (user+sys) | CPU-s/Mb masked | unmasked v3 | decode | peak RSS |
|---|---|---|---|---|---|
| S. cerevisiae chr I | 3.39 | **12.14** | 12.38 | 1.47 | 0.96 GiB |
| C. elegans chr V | 181.53 | **8.65** | 8.88 | 100.90 | 2.07 GiB |

In this single run the masked decode is cheaper on chr V by **2.64%**
end to end (stage sum 8.8798 → 8.6454 CPU-s/Mb), with the **decode stage
alone 4.34%** lower (105.478 → 100.896 CPU-s); chr I moves −1.97%. The
mask deletes transitions but also builds the mask on the bias path and
leaves the scan kernels unchanged, so these are observations from one
run per chromosome, not an established speedup (stalin-0114). The
section 5 CPU verdict is unchanged and **no positive CPU allowance for B
follows from this run either**.

## Accuracy, unmasked v3 against the same checkpoint masked

| metric | chr I v3 | chr I masked | chr V v3 | chr V masked |
|---|---|---|---|---|
| nucleotide F1 | 0.891 | 0.890 | **0.598** | 0.585 |
| nt sensitivity / precision | 0.912 / 0.870 | 0.912 / 0.869 | 0.533 / 0.680 | 0.503 / 0.698 |
| locus F1 | 0.786 | 0.778 | 0.603 | 0.590 |
| fusion / split | 3 / 0 | 3 / 0 | 312 / 910 | **226** / 964 |
| exact transcripts (sens) | 64 (0.681) | **65** (0.691) | 70 (0.014) | **101** (0.020) |
| exon F1 (exact) | 0.582 | 0.588 | 0.030 | **0.062** |
| donor F1 | 0 | 0 | 0.070 | **0.117** |
| acceptor F1 | 0 | 0 | 0.081 | **0.158** |
| GT-AG introns TP / FP / FN | 0 / 0 / 3 | 0 / 12 / 3 | 433 / 470 / 22,187 | **1,070** / 1,970 / 21,550 |
| GT-AG intron precision | – | 0 | 0.480 | 0.352 |
| GC-AG introns TP / FP / FN | 0 / 0 / 0 | 0 / 3 / 0 | 1 / 29 / 164 | 1 / 601 / 164 |
| reference introns outside GT-AG/GC-AG matched | 0 of 0 | 0 of 0 | 1 of 39 | **0** of 39 |
| predicted introns outside GT-AG/GC-AG | 16 | **0** | 3,998 | **0** |
| predicted introns (median length) | 16 (81.5 b) | 15 (63 b) | 4,931 (615 b) | 3,642 (550 b) |
| predicted introns above reference q90 = 695 b | – | – | 2,347 (47.60%) | 1,676 (46.02%) |

The transcript sensitivity column is the scorer's own field, exact
transcripts over the scored reference loci (chr V: 4,995). `../score-final/`
and a-pilot section 3.3 divide the same counts by all 6,766 reference
transcripts instead, which turns chr V's 70 and 101 into 0.010 and 0.015.
The counts are the same numbers either way.

The intron rows use the scorer's unit: CDS rows of the committed GFF3
grouped by sequence, strand and `Parent`, introns of at least 20 b
deduplicated by sequence/start/end/strand, which reproduces
`splice.predicted_introns` exactly (4,931 unmasked, 3,642 masked). The
decile cuts `[45, 46, 48, 51, 56, 95, 199, 370, 695]` are the scorer's
own, over unique reference introns.

What the mask does, stated as what the table supports:

- **Splice placement improves by about a factor of two, and not only by
  deletion.** Donor F1 0.070 → 0.117 and acceptor F1 0.081 → 0.158 on
  chr V; correct GT-AG introns 433 → 1,070. The improvement is not only
  the removal of the 3,998 predictions outside GT-AG/GC-AG: the number
  of *correct* introns rises, so masking moves probability onto real
  sites rather than merely deleting predictions.
- **The restriction has its own costs.** Inside the retained classes,
  GT-AG false positives rise 470 → 1,970 and GT-AG precision falls
  0.480 → 0.352; GC-AG keeps its single true positive while its false
  positives rise 29 → 601; and the one reference intron outside both
  classes that unmasked v3 matched is now unreachable. The net gain is
  real; the excluded junctions are not thereby shown to be impossible.
- **Gene structure improves.** Exact transcripts on chr V rise 70 → 101
  (+44%, equal to the v2 checkpoint's count), exact exon F1 doubles
  0.030 → 0.062, and fusions fall 312 → 226.
- **Intron duration is still miscalibrated.** The masked decode predicts
  3,642 unique introns with median 550 b, 46.02% of them above the
  reference q90 of 695 b (unmasked v3: 4,931, 615 b, 47.60%; reference
  q50 is 56 b). The mask selects among supported sites and does not
  address the length distribution (engels-0110, stalin-0114).
- **Nucleotide F1 falls slightly** on chr V, 0.598 → 0.585: sensitivity
  0.533 → 0.503 against precision 0.680 → 0.698, and predicted CDS falls
  4.40 Mb → 4.04 Mb. Forbidden splice sites cost the decoder coding bases
  it used to claim through illegal introns; what it keeps is more often
  right. Locus F1 moves with it (0.603 → 0.590).
- **chr I is unchanged**, as expected where there is almost nothing to
  splice: nucleotide F1 0.891 → 0.890, the 16 predicted introns outside
  GT-AG/GC-AG replaced by 15 inside them (none correct either way),
  64 → 65 exact transcripts.
- **A still misses the accuracy target.** Exact-transcript sensitivity on
  chr V is 0.020 per scored locus, 0.015 per reference transcript.
  Nothing here accepts T-human-014. 21,550 reference GT-AG introns are
  still missed, so the remaining failure is *which* supported site the
  encoder scores highest — the next decoder increment, and the reason a
  mask alone is not the fix. Under the accepted unmasked support that
  question still includes the sites this ablation removes.

The mask is inference-only. A hard `-inf` in the chain-loss numerator
would make a reference intron with an unusual dinucleotide unreachable
and the loss infinite, so the fitting paths do not take it
(`model.a.pooled.motif_bias`, `tests/test_a_pooled.py CanonicalSpliceMask`,
and `measure` refuses `--canonical-splice` outside the chromosome
profile so no unmasked row can be reported as masked).

## Host conditions

One core (core 2), `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, the same
24-core Intel Core Ultra 9 285K as v1–v3; load average 0.31 at launch, no
other work of mine on the machine. An idle core is an operator
precondition, not something the gate can establish (engels-0108); the
cost columns are comparable with `../score-final/`'s on that basis.

A first launch of this same script was stopped by hand during the chr V
measurement, because the measured row did not yet carry the mask fields;
the wrapper reported `measure NC_003283.11 exit=143` and did not score,
which is the stalin-0113 failure propagation working on a real failure.
Its outputs were deleted, not reused: every file here comes from the
single clean run recorded in `run_canonical.out`.

## Files

`run_canonical.sh`, `run_canonical.out`, `declaration.yaml` (the
score-final declaration plus the mask), the two seqid lists,
`measure_*.json` (each with `canonical_splice: true`), `score_*.json`,
`/usr/bin/time -v` reports for all four workloads, the two predicted
GFF3s (chr V gzipped), `sha256.txt`, and `summarize.py`.

`sha256.txt` also lists this README. Its entry is refreshed whenever the
prose is corrected in a later tick (2026-09-22T17:12Z: the
restricted-support framing, the per-class costs, the duration rows and
the cost split, from engels-0110 and stalin-0114). Every run output
listed there is byte-identical to the clean run; only the README line
moves.
