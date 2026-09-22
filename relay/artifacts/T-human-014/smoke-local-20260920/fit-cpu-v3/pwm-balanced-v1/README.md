# pwm-balanced-v1/: the source-balanced splice-site PWM (section 3.6 increment 3)

The same donor `[-3, +6)` / acceptor `[-20, +3)` log-odds matrices as
`../pwm-v1/`, estimated over the same 67,325 admitted train junctions
of the same development-excluded split, with one difference: every
**source** carries the same weight instead of every **junction**.

`model.a.splicepwm.fit_grouped` averages the per-source smoothed column
frequencies and the per-source ACGT backgrounds, unweighted, so neither
the matrix nor its null model is decided by whichever species has the
most introns. `model.a.train pwm --balanced` fits it; the JSON records
`weighting: balanced` and, per source, the `train_sites` it counted.

**The imbalance this was meant to correct**, from `pwm.out`:

| source | train windows | train junctions |
|---|---:|---:|
| Saccharomyces_cerevisiae | 5,774 | 265 |
| Caenorhabditis_elegans | 14,201 | 67,060 |

253:1. Totals (67,325 junctions over 19,975 train windows, 4,951
development windows excluded) are identical to `../pwm-v1/`, so the two
matrices differ only in weighting.

**Leakage.** As in `../pwm-v1/`: the fit loads the run's own config and
calls `split_windows` with every source's `dev_seqids`, so the two
scored development chromosomes — S. cerevisiae `NC_001133.9` and
C. elegans `NC_003283.11` — contribute no junction and no background
base. The JSON carries the excluded seqids and the pinned GFF3/FASTA
MD5s. No held-out species is touched.

**It did not work.** See `../score-pwm-bal/` and section 3.6 of
`docs/design/a-pilot.md`: balancing costs both chromosomes exactness.
The matrix is kept as the artifact of a negative result, not as a
recommended default.

Files: `splicepwm.json` (sha256 `067b55e3…`, the matrix the two scored
runs pin), `pwm.out` (the fitting log, including the per-source counts),
`pwm_time.txt` and `fit_time.txt` (`/usr/bin/time -v` of the logged
refit and of the original fit), and `sha256.txt`. The original fit was
run without capturing stdout, so `pwm.out` is a refit, logged verbatim:
it writes `splicepwm.rerun.json`, which came out byte-identical to the
`splicepwm.json` kept here and to the digest the two scored runs pin.
