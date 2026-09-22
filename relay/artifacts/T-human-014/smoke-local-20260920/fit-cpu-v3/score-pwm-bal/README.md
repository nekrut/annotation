# score-pwm-bal/: the fit-cpu-v3 checkpoint decoded with the **source-balanced** PWM and the mask, 2026-09-22

Section 3.6, increment 3. Section 3.5 found that the pooled splice-site
PWM doubles splice placement on C. elegans chr V and **hurts** the
nearly intronless S. cerevisiae chr I, and explained the second finding
by the matrix being dominated by nematode junctions. stalin-0117
objected that the runs establish the regression, not its cause. This
run is the direct test: the same increment with the imbalance removed.

`../pwm-balanced-v1/splicepwm.json` (sha256 `067b55e3…`,
`weighting: balanced`) gives each source the same weight. **Nothing else
changes**: the same `best.pt` (step 3,000 of 3,000, dev NLL 28.640),
the same tiling (`--segments 19 --overlap 4096 --window 12288
--dtype float64`), the same scorer, the same pinned core, no refit,
inference only. `run_pwm_balanced.sh masked` is `run_pwm.sh` with the
matrix path and its digest gate changed.

## Result: the imbalance hypothesis fails

| metric | pooled (`../score-pwm/`) | **balanced (here)** |
|---|---|---|
| chr I exact transcripts (of 94) | 55 | **47** |
| chr I exon exact F1 | 0.3959 | **0.3300** |
| chr I GT-AG introns TP/FP (3 reference) | 2/67 | **3/83** |
| chr I nucleotide F1 | 0.8569 | **0.8610** |
| chr I locus F1 | 0.7215 | **0.7308** |
| chr V exact transcripts (of 6,766) | 835 | **706** |
| chr V correct GT-AG introns (of 22,620) | 18,195 | **16,963** |
| chr V donor / acceptor F1 | 0.642 / 0.666 | **0.628 / 0.646** |
| chr V exon exact F1 | 0.5216 | **0.5008** |
| chr V nucleotide F1 | 0.7675 | **0.7579** |
| chr V locus F1 | 0.6386 | **0.6480** |

Yeast chr I — the chromosome balancing was for — loses eight more exact
transcripts and gains sixteen more false GT-AG introns. It does recover
the third reference intron, and its nucleotide and locus F1 rise
slightly, but the exactness metrics the regression was reported in get
worse. chr V pays the expected price of a 253:1 reweighting.

So **species imbalance is not a sufficient explanation** of the chr I
regression. The surviving hypothesis is the one stalin-0117 named: the
strength and calibration of a fixed extended-context term added to
frozen learned emissions. On a chromosome with three introns in 230 kb,
a confident splice score that is not opposed by an equally calibrated
intron-entry cost buys recall with false junctions whatever species it
was counted on. The next experiment is a scalar weight on the PWM term
(and the entry hazard it competes with), swept and reported as a curve.

## Cost, unchanged

chr V 8.99 CPU-s/Mb (user-only 7.9, user+system 9.0), chr I 12.46
(user-only 12.7, user+system 15.2), peak host RSS 2.05 GiB and 0.96
GiB. Far under the accepted 15 CPU-s/Mb and 8 GB, so **no positive CPU
allowance for B follows**, and **A still misses the accuracy target**.

## Provenance

`measure_*.json` record `commit 9ec0424` with `source_dirty: true` and
`source_sha256 3e433c82…`; that digest is byte-identical to the clean
tree at `3fd9032` on `work/T-human-014-lenin`, so the executed source is
exactly what that commit contains. An earlier chr I workload recorded a
different digest because a docstring was edited while it ran; both
masked workloads were re-run afterwards and the pair here shares one
digest. Accuracy was identical across the two runs; cost differed by
0.7% on chr V and chr I, which is the run-to-run spread on this machine.

`../score-pwm-bal-only/` is the same matrix without the hard mask.
`summarize.py` regenerates every number above from this directory.
