# baselines/kaks

The KA/KS comparative baseline of T-human-010: the test of Nekrutenko,
Makova and Li (2002, *Genome Research* 12:198-202, doi:10.1101/gr.200901)
as a sliding-window coding classifier on the benchmark's alignment windows,
the floor every later model has to beat. Standard library only.

| file | what it is |
|---|---|
| `kaks.py` | pairwise KA/KS by the Nei and Gojobori (1986) pathway method with the Jukes-Cantor correction, the one-sided z-test of dN < dS, and the decision rule (ratio below one and significant). Standard and ciliate genetic codes. `--self-test`: 22 checks |
| `windows.py` | reads a window written by `scripts/data/fetch_window.py`, pairs the reference with one informant row, tests every window of `W` bases in six frames, projects the calls to bases and scores them against the CDS annotation at the nucleotide level (`docs/benchmark.md` section 4.1, strand-aware), with the metrics restricted to aligned bases, sensitivity by CDS segment length, and the count of windows the test could not decide. `--self-test`: 12 checks on a synthetic two-gene window |
| `run.py` | the single command: for every row of `pairs.tsv` draws the seeded gene sample (`scripts/data/sample_genes.py`), fetches each gene's window into a cache (`scripts/data/fetch_window.py`, skipped when cached), runs the sweep over informants and window sizes, and pools the confusion counts over genes into `results/` |
| `pairs.tsv` | the four species pairs of the sweep: reference assembly, track, chromosome, informants, gene sample parameters |
| `results/` | `summary.tsv` (pair x informant x window, pooled), `per_gene.tsv`, `summary.md`, `manifest.json` (sample headers with the track `dataTime`, SHA-256 of every fetched input, informants absent from a window) |

```
python3 baselines/kaks/run.py --cache /tmp/kaks-cache        # the deliverable: 4 pairs x 8 genes, ~10 min, ~200 MB cache
python3 baselines/kaks/run.py --cache /tmp/kaks-cache --offline --pairs-filter fly --limit-genes 2   # a bounded rerun
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

## Results (4 pairs x 8 genes, 2026-09-10)

`results/summary.md` is the full sweep: 4 references, 19 informants, 4
window sizes, pooled over the 8 genes of each seeded sample
(`sample_genes.py`, seed 20260909: non-overlapping RefSeq genes with
complete CDS and at least 3 exons, transcript length 3-12 kb for human and
mouse, 2-8 kb for fly and worm, 500 bp flank). `results/per_gene.tsv` has
the per-gene rows behind every pooled number, `results/manifest.json` the
track `dataTime` of each sample and the SHA-256 of every fetched file. The
run took 365 s wall clock, of which the fetch was nearly all; the cache is
92 MB.

The best window size per informant, by pooled nucleotide F1
(`docs/benchmark.md` 4.1, strand-aware, over every base of every window):

| pair | informant | distance | aligned | best W | Sens | Prec | F1 | MCC | Sens <100 bp | Sens 400+ bp |
|---|---|---|---|---|---|---|---|---|---|---|
| human chr17 | panTro4 | 0.01 | 0.98 | any | 0.00 | | | | 0.00 | 0.00 |
| human | rheMac3 | 0.08 | 0.92 | 600 | 0.29 | 0.15 | 0.20 | 0.13 | 0.00-0.16 | 0.65 |
| human | canFam3 | 0.33 | 0.61 | 150 | 0.47 | 0.68 | 0.55 | 0.54 | 0.00 | 0.77 |
| human | mm10 | 0.50 | 0.52 | 150 | 0.62 | 0.67 | 0.64 | 0.62 | 0.17-0.26 | 0.83 |
| human | galGal4 | 1.17 | 0.16 | 600 | 0.39 | 0.25 | 0.31 | 0.25 | 0.57-0.87 | 0.24 |
| human | danRer10 | 2.07 | 0.11 | 600 | 0.21 | 0.25 | 0.23 | 0.18 | 0.17-0.28 | 0.04 |
| mouse chr19 | rn6 | 0.17 | 0.94 | 600 | 0.64 | 0.24 | 0.35 | 0.29 | 0.42-0.65 | 1.00 |
| mouse | hg38 | 0.52 | 0.65 | 150 | 0.62 | 0.71 | 0.66 | 0.63 | 0.00-0.17 | 0.81 |
| mouse | canFam4 | 0.58 | 0.60 | 150 | 0.66 | 0.71 | 0.68 | 0.65 | 0.00-0.23 | 0.85 |
| mouse | galGal6 | 1.31 | 0.13 | 300 | 0.27 | 0.42 | 0.33 | 0.28 | 0.14-0.21 | 0.17 |
| fly chr2L | droSim2 | 0.10 | 0.98 | 600 | 0.79 | 0.73 | 0.76 | 0.69 | 0.24-0.59 | 0.87 |
| fly | droYak3 | 0.23 | 0.95 | 300 | 0.85 | 0.78 | 0.81 | 0.76 | 0.00-0.42 | 0.94 |
| fly | droAna3 | 1.05 | 0.80 | 600 | 0.84 | 0.39 | 0.53 | 0.39 | 1.00 | 0.80 |
| fly | droPse3 | 1.59 | 0.67 | 300 | 0.75 | 0.53 | 0.62 | 0.51 | 0.03-0.62 | 0.82 |
| fly | anoGam3 | 3.84 | 0.39 | 600 | 0.58 | 0.37 | 0.45 | 0.26 | 0.13-0.54 | 0.60 |
| worm chrIII | caeJap4 | 1.12 | 0.44 | 600 | 0.31 | 0.22 | 0.26 | 0.02 | 0.00-0.40 | 0.05 |
| worm | caeRem4 | 1.14 | 0.59 | 600 | 0.57 | 0.28 | 0.37 | 0.14 | 0.67-0.84 | 0.53 |
| worm | caeAng2 | 1.25 | 0.48 | 600 | 0.36 | 0.28 | 0.32 | 0.10 | 0.00-0.40 | 0.37 |
| worm | C_briggsae | 1.31 | 0.60 | 600 | 0.60 | 0.29 | 0.39 | 0.16 | 0.58-0.84 | 0.53 |

Distance is patristic from the reference on the track's own tree
(substitutions per site). `aligned` is the fraction of window bases at
which the informant has a base. The last two columns are sensitivity for
CDS bases in segments under 100 bp (the two strata, `<60` and `60-100`)
and in segments of 400 bp or more.

The genes behind the four samples are not alike, which matters for
reading the table. CDS is 13% of the human windows and 19% of the mouse
windows, against 44% (fly) and 43% (worm); 48% of human CDS bases and 46%
of mouse CDS bases lie in segments shorter than 200 bp, against 13% of
fly and 28% of worm CDS bases (computed from the fetched annotations with
`windows.truth_arrays`; the per-pair numbers are in the task log).

**The floor.** For a mammal with an informant in the useful band and a
window of 150 bp, the test reaches nucleotide F1 0.64-0.68 (MCC
0.62-0.65) over whole windows. For fly with *yakuba* it reaches F1 0.81
(MCC 0.76) at 300 bp. Every later model is measured against the same
metric on the same windows (`results/per_gene.tsv` has the window
coordinates), and has to beat these numbers by a margin before its extra
machinery is justified. The worm result (F1 at most 0.39, MCC at most
0.16) is not a floor; it is a failure, explained below.

## Where the signal fails

1. **Too close: no substitutions, no test.** Chimp (0.01) leaves 776 of
   the 1,975 90-bp human windows without a single difference and calls
   nothing at any window size. Rhesus (0.08) reaches F1 0.20. *simulans*
   (0.10) needs 600 bp to reach F1 0.76 and calls nothing useful at 90 or
   150 bp. The test counts synonymous and nonsynonymous differences; with
   only a handful of differences in a window the z-test cannot reject
   dN = dS whatever the ratio. A model that wants to use close relatives
   (which is where the alignments are complete: chimp covers 98% of the
   human windows) must use a signal that does not need substitutions,
   such as the absence of indels that break frame, and must pool
   substitutions across many close informants instead of testing one.

2. **Too far: the alignment goes first, then the frame.** Coverage falls
   with distance in every clade: human-mouse 52%, human-chicken 16%,
   human-zebrafish 11%; mouse-chicken 13%; fly-*Anopheles* 39%. The
   unaligned bases are mostly non-coding (exons are what survives), so
   coverage loss shows up as precision, not sensitivity: at 300 bp,
   human-mouse sensitivity is 0.77 but precision 0.34, because a called
   window spills into flanking intron that the informant does not align.
   Restricted to aligned bases (`aligned F1` in `results/summary.md`) the
   numbers improve only a little, because the spill is inside the aligned
   part too. Beyond about 1 substitution per site, the second failure
   appears: the pairs that are aligned are near saturation at synonymous
   sites, so dS is large and noisy and dN/dS is below one in shadow frames
   nearly as often as in the true frame (fly *ananassae* at 1.05 and
   600 bp: sensitivity 0.84, precision 0.39, MCC 0.39). The worm pairs are all in this regime
   (1.12-1.31 on the 135-way tree; there is no *Caenorhabditis* at 0.2-0.6
   from *elegans* in the track), and their aligned-only MCC is 0.14 or
   below: the test is at noise level.

3. **Short exons: the window is the wrong unit.** At the 150-bp window
   that maximizes F1 for mammals, sensitivity is 0.00-0.26 for CDS segments
   under 100 bp and 0.44-0.63 for 100-200 bp, against 0.77-0.85 for
   segments of 400 bp or more. Nearly half of mammalian CDS bases are in
   segments under 200 bp (48% human, 46% mouse in these samples), so
   half the coding sequence is invisible to the window that works best on
   the other half. Raising the window to 600 bp recovers the short exons
   (sensitivity 0.47-1.00 under 100 bp) but precision collapses to
   0.17-0.23, because the window is now mostly intron. There is no window
   size that resolves both; the sweep's best F1 is a compromise between
   two failure modes, not a good operating point. The paper's own remark
   that the test cannot locate exon boundaries is the same fact.

4. **Compact genomes: introns inside the window break the frame.** Worm
   introns have a median of 65 bp and fly introns 102 bp
   (`benchmark/panel.tsv`), so a 300-bp window in a worm gene usually
   contains an intron, and the intron shifts the reading frame of the
   downstream exon relative to the upstream one. The six-frame test then
   sees a window whose two halves are coding in different frames, and
   neither frame passes. Fly escapes this because 72% of its CDS bases
   are in segments of 400 bp or more (long exons, few introns per gene in
   these samples); worm does not (39%). This, together with the distance
   problem, is why the worm result is a failure rather than a floor.

5. **The stop-codon veto is unused, and precision shows it.** The default
   run does not veto frames with reference stop codons (`--max-ref-stops`
   unset), so a window with dN/dS below one in a frame that contains a
   stop is still called. The shadow-frame false positives of item 2 and
   the spill of item 3 would both be reduced by the veto; it was left off
   so the floor is the 2002 test and nothing more. The synthetic self-test
   shows that shadow frames pass when substitutions sit only at third
   positions.

Not tested here: fast-evolving genes (the sample is a random draw, not a
selection by rate) and conserved non-coding elements as a source of false
positives (the windows are gene loci with 500 bp of flank, not intergenic
space). Both need the benchmark's whole-chromosome runs, which the
`--out` GFF3 of `windows.py` is written for.

## What this says about the inputs the new model needs

- **Many informants at once, spanning distances, with the tree.** No
  single informant is in the useful band for every base: the close ones
  align everywhere and carry no substitutions, the far ones carry
  substitutions and align only at exons. The information is in the sum
  over the tree, which is what an alignment-plus-tree input gives a model
  that treats phylogeny as a coordinate (the HyphAeon pattern in
  `relay/TASK.md`). A pairwise design inherits this section's failure
  modes 1 and 2 by construction.
- **Codon-resolution, not window-resolution, with an explicit frame.**
  The unit of evidence has to be the aligned codon column, and the frame
  has to be a state that a decoder carries across an intron, so that a
  40-bp exon contributes its 13 codons to a coding call instead of being
  averaged into 260 bp of intron. Failure modes 3 and 4 are properties of
  fixed windows, not of the KA/KS signal.
- **Boundary signals in the same input.** The test cannot see splice
  sites; the model must, from the reference sequence itself and from
  where the informants' alignment gaps and frame shifts sit. Alignment
  coverage boundaries are themselves informative (coverage loss is
  concentrated in introns), so the model should see "unaligned" as a
  value, not as missing data.
- **Precision comes from what the test does not use.** Reference stop
  codons in frame, frame-breaking indels in the informants, and the
  requirement that a coding segment be open in one frame from one splice
  site to the next are all cheap and all absent from the floor. Any model
  that has them should beat the precision column easily; the sensitivity
  column on short exons is the harder target.
- **Clade coverage of the alignments is the real limit for worms and
  everything compact.** The 135-way nematode track has no informant at a
  useful distance from *elegans*, and no UCSC alignment exists at all for
  14 of the 20 panel species (`docs/data-sources.md`). The model's
  training plan needs to say where the alignments for those come from,
  or admit a DNA-only path for them.
