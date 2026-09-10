# baselines/kaks

The KA/KS comparative baseline of T-human-010: the test of Nekrutenko,
Makova and Li (2002, *Genome Research* 12:198-202, doi:10.1101/gr.200901)
as a sliding-window coding classifier on the benchmark's alignment windows,
the floor every later model has to beat. Standard library only.

| file | what it is |
|---|---|
| `kaks.py` | pairwise KA/KS by the Nei and Gojobori (1986) pathway method with the Jukes-Cantor correction, the one-sided z-test of dN < dS, and the decision rule (ratio below one and significant). Standard and ciliate genetic codes. `--self-test`: 22 checks |
| `windows.py` | reads a window written by `scripts/data/fetch_window.py`, pairs the reference with one informant row, tests every window of `W` bases in six frames, projects the calls to bases and scores them against the CDS annotation at the nucleotide level (`docs/benchmark.md` section 4.1, strand-aware), with the metrics restricted to aligned bases, sensitivity by CDS segment length, and the count of windows the test could not decide. `--self-test`: 12 checks on a synthetic two-gene window |

```
python3 baselines/kaks/kaks.py
python3 baselines/kaks/windows.py --self-test
python3 scripts/data/fetch_window.py --assembly dm6 --locus chr2L:14615552-14618902 --flank 500 --track multiz124way --out /tmp/win/Adh --name Adh_124
python3 baselines/kaks/windows.py --stem /tmp/win/Adh/Adh_124 --informants droSim2,droYak3,droPse3 --windows 90,150,300 --out /tmp/kaks/Adh
```

## The test, and how this reproduction differs

The original study fitted KA/KS by maximum likelihood in codeml and tested
it with a likelihood-ratio statistic (the reading of the paper in relay
note 20260909T164106Z-stalin-0016; the full text is not open through
Europe PMC or PMC, only the abstract, and the publisher PDF is a page
image). codeml is a C dependency outside the charter's stack, so `kaks.py`
uses the counting estimator and its z-test instead. The decision rule is
the same: a window is coding when KA/KS is below one *and* the difference
is significant. Because a window is called when any of six frames passes,
each frame is tested at `alpha / 6` (`--no-bonferroni` turns that off).

A window the test cannot decide is not a negative. `windows.py` keeps
them apart: `identical` (no difference between the two sequences at all),
`saturated` (synonymous differences at or beyond three quarters of the
synonymous sites, so dS is undefined) and `too_few_codons` (fewer aligned
codon pairs than `--min-codons`, 15 by default, which is also what an
unaligned informant produces). This is the accounting stalin's note asks
for: alignment availability, frame selection and the fitted test are
reported separately.

## Preliminary numbers (two windows, 2026-09-10)

Not the deliverable; the sweep over a gene sample per species pair comes
next. Fly *Adh* window (dm6 124-way, `chr2L:14615052-14619402`, 1,590 CDS
bases on the plus strand, none on the minus) and human *TP53* window (hg38
100-way, `chr17:7667921-7687990`, 1,263 CDS bases on the minus strand).
`aligned` is the fraction of window bases at which the informant has a
base; window status counts are for the run's windows; the four metrics are
the strand-aware nucleotide ones over the whole window.

| Window | Informant | Distance | W | aligned | called / tested / identical / saturated / too few | Sens | Prec | F1 | MCC |
|---|---|---|---|---|---|---|---|---|---|
| Adh | droSim2 | 0.10 | 90 | 0.89 | 0 / 119 / 11 / 0 / 13 | 0.00 | | | |
| Adh | droSim2 | 0.10 | 300 | 0.89 | 3 / 37 / 0 / 0 / 1 | 0.40 | 0.91 | 0.56 | 0.56 |
| Adh | droYak3 | 0.23 | 150 | 0.88 | 16 / 63 / 0 / 0 / 6 | 0.70 | 0.82 | 0.76 | 0.71 |
| Adh | droYak3 | 0.23 | 300 | 0.88 | 13 / 27 / 0 / 0 / 1 | 0.85 | 0.71 | 0.77 | 0.72 |
| Adh | droAna3 | 1.05 | any | 0.00 | 0 / 0 / 0 / 0 / all | 0.00 | | | |
| Adh | droPse3 | 1.59 | 150 | 0.73 | 23 / 50 / 0 / 0 / 12 | 0.62 | 0.52 | 0.56 | 0.46 |
| Adh | droPse3 | 1.59 | 300 | 0.73 | 16 / 22 / 0 / 0 / 3 | 0.66 | 0.42 | 0.51 | 0.39 |
| Adh | anoGam3 | 3.84 | any | 0.00 | 0 / 0 / 0 / 0 / all | 0.00 | | | |
| TP53 | panTro4 | 0.01 | 90 | 0.93 | 0 / 358 / 263 / 0 / 45 | 0.00 | | | |
| TP53 | panTro4 | 0.01 | 300 | 0.93 | 2 / 171 / 14 / 0 / 11 | 0.00 | 0.00 | 0.00 | -0.02 |
| TP53 | mm10 | 0.50 | 150 | 0.30 | 7 / 130 / 0 / 0 / 262 | 0.42 | 0.62 | 0.50 | 0.50 |
| TP53 | mm10 | 0.50 | 300 | 0.30 | 13 / 67 / 0 / 0 / 118 | 0.76 | 0.38 | 0.51 | 0.52 |
| TP53 | galGal4 | 1.17 | 300 | 0.07 | 4 / 21 / 0 / 0 / 173 | 0.32 | 0.51 | 0.40 | 0.39 |
| TP53 | danRer10 | 2.07 | 300 | 0.04 | 0 / 15 / 0 / 0 / 183 | 0.00 | | | |

What the two windows already show, all of it predicted by the paper's
"sufficiently long and a suitable degree of divergence":

- **Too close is blind.** Chimp against human (0.01) leaves 263 of 358
  90-bp windows without a single difference and calls nothing; fly
  *simulans* (0.10) calls nothing at 90 bp and 3 windows at 300 bp.
- **The useful band is roughly 0.2 to 0.6 substitutions per site**, where
  *yakuba* (0.23) reaches F1 0.77 and mouse (0.50) 0.50, the latter capped
  by coverage: mouse has a base at 30% of the human window, so 118 to 463
  windows are untestable and the untestable ones are the non-coding
  majority, which is why precision, not sensitivity, is the weak number.
- **Distant informants lose the alignment before they lose the signal.**
  Chicken (1.17) and zebrafish (2.07) cover 7% and 4% of the human window;
  *ananassae* (1.05) and *Anopheles* (3.84) have no row at *Adh* at all in
  this track. *pseudoobscura* (1.59) still covers 73% of *Adh* but
  precision falls to 0.4 to 0.5, which is the shadow-frame and
  near-saturation regime.
- **Short windows have no power.** At 90 bp (30 codons) nothing is called
  on any pair at `alpha / 6`; the test needs 150 to 300 bp, so exons
  shorter than that can only be found by a window that spills over their
  boundaries, and the paper's own remark that the test cannot locate
  exon/intron boundaries applies. The per-length sensitivity in the JSON
  output is the number to report once the gene sample is in.

## Next

A driver that takes a species-pair table (reference assembly, track,
informant, gene sample drawn with `scripts/data/sample_genes.py`),
fetches the windows into a cache, runs the sweep over window size and
informant distance, and writes `results/` plus the section 4.8 style table
from one command, then the write-up of where the signal fails and what
that says about the inputs the new model needs.
