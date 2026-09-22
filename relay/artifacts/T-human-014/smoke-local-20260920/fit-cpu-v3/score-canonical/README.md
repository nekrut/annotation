# score-canonical/: the fit-cpu-v3 checkpoint decoded with the hard splice mask, 2026-09-22T16:25Z

The section 3.4 **decoder** ablation. Same weights (`best.pt` sha256
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

The mask is free: it removes transitions rather than adding work, and the
masked decode is about 4% *cheaper* on chr V (decode 105.48 → 100.90
CPU-s). The section 5 CPU verdict is unchanged and **no positive CPU
allowance for B follows from this run either**.

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
| non-canonical predicted introns | 16 | **0** | 3,986 | **0** |
| predicted introns (median length) | 16 (81.5 b) | 15 (63 b) | 4,931 (615 b) | 3,642 (550 b) |

The transcript sensitivity column is the scorer's own field, exact
transcripts over the scored reference loci (chr V: 4,995). `../score-final/`
and a-pilot section 3.3 divide the same counts by all 6,766 reference
transcripts instead, which turns chr V's 70 and 101 into 0.010 and 0.015.
The counts are the same numbers either way.

What the mask does, stated as what the table supports:

- **Splice placement improves by about a factor of two at no cost.**
  Donor F1 0.070 → 0.117 and acceptor F1 0.081 → 0.158 on chr V; correct
  GT-AG introns 433 → 1,070. The improvement is not only the removal of
  the 3,986 non-canonical predictions: the number of *correct* introns
  rises, so masking the illegal sites moves probability onto real ones
  rather than merely deleting predictions.
- **Gene structure improves.** Exact transcripts on chr V rise 70 → 101
  (+44%, equal to the v2 checkpoint's count but now with v3's plausible
  intron lengths), exact exon F1 doubles 0.030 → 0.062, and fusions fall
  312 → 226.
- **Nucleotide F1 falls slightly** on chr V, 0.598 → 0.585: sensitivity
  0.533 → 0.503 against precision 0.680 → 0.698, and predicted CDS falls
  4.40 Mb → 4.04 Mb. Forbidden splice sites cost the decoder coding bases
  it used to claim through illegal introns; what it keeps is more often
  right. Locus F1 moves with it (0.603 → 0.590).
- **chr I is unchanged**, as expected where there is almost nothing to
  splice: nucleotide F1 0.891 → 0.890, 16 non-canonical introns replaced
  by 15 canonical ones, 64 → 65 exact transcripts.
- **A still misses the accuracy target.** Exact-transcript sensitivity on
  chr V is 0.020 per scored locus, 0.015 per reference transcript.
  Nothing here accepts T-human-014. 21,550 reference GT-AG introns are
  still missed, so the remaining failure is not
  legality but *which* legal site the encoder scores highest — the next
  decoder increment, and the reason a mask alone is not the fix.

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
