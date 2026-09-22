# score-final/: development-chromosome scores under the finished fit-cpu-v1 `best.pt` (step 900)

Run 2026-09-22T04:05Z on core 2, one thread, code `aeed615` (PR #38; the
`model/a` + `model/grammar` source digest is unchanged since the fit's
`4c43819`, see `measure_*.json` `source`). `run_final.sh` runs, per
chromosome, `measure --profile chromosome --window 12288 --overlap 4096
--segments 19 --dtype float64 --checkpoint best.pt --gff-out` under
`/usr/bin/time -v`, then `benchmark/score.py --seqids --genome
--declaration`. All four exits 0 (`run_final.out`). `sha256.txt` pins the
two predicted GFF3s and the host-local `best.pt`
(`181adad7…4cf88b`, the fit's `best_pt.sha256`). The chr V GFF3 is
committed gzipped (796 KB raw).

## Cost (checkpoint-independent apart from the decoded chain count)

| chromosome | bases | tiles | CPU-s | CPU-s/Mb | encoder | decode | preprocess | io | output | peak RSS | process wall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *S. cerevisiae* chr I | 230,218 | 76 | 2.78 | **12.09** | 1.26 | 1.48 | 0.033 | 0.008 | 0.0002 | 0.96 GiB | 3.4 s |
| *C. elegans* chr V | 20,924,180 | 3,420 | 185.4 | **8.86** | 77.6 | 105.0 | 2.33 | 0.44 | 0.015 | 2.05 GiB | 186.1 s |

`score.py`: chr I 0.13 CPU-s / 32 MB; chr V 2.84 CPU-s / 173 MB.

Budget reading (correction after stalin-0103): 15 CPU-s/Mb is bound to the
cost-baseline runner; this machine's portable ceiling is AUGUSTUS
51.41 / 11 = **4.67 CPU-s/Mb** (`../../pombe-normalization/`). Stage sum /
user / user+system: chr I 12.09 / 12.42 / 14.64 = 2.59× / 2.66× / 3.13×
the ceiling; chr V 8.86 / 7.77 / 8.89 = 1.90× / 1.66× / 1.90×.
Machine-normalized (×3.21) the stage sums are 38.8 and 28.4 against 15.
**Both fitted rows miss the CPU budget**; memory passes.

## Accuracy (`score_*.json`; development chromosomes of seen species, declaration in `declaration.yaml`)

| metric | chr I (94 ref transcripts, 3 ref introns) | chr V (6,766 ref transcripts / 4,995 loci, 22,824 ref introns) |
|---|---|---|
| predicted transcripts | **6** | **287** |
| nucleotide sens / prec / F1 / MCC | 0.039 / 0.356 / 0.071 / 0.020 | 0.470 / 0.586 / 0.522 / 0.461 |
| predicted / reference CDS bp | 15,537 / 141,204 | 4,491,933 / 5,601,813 |
| locus TP / FP / FN, fusion, split | 3 / 3 / 91, 2, 0 | 269 / 18 / 4,726, **238**, 121 |
| transcript exact | 0 | 0 |
| exon exact (all) | 0 / 37 FP | 67 TP / 10,481 FP / 28,502 FN |
| predicted introns by `score.py` dinucleotide category | 31: 1 GT-AG / 30 other | 10,261: 586 GT-AG (162 TP + 424 FP) / 43 GC-AG / 10 AT-AC / **9,622 other** (9,675 non-GT-AG) |
| start / stop codon TP | 1 / 0 | 7 / 3 |
| predicted chain span min / median / max | 234 / 11,212 / 15,989 | 147 / **81,233** / **842,820** |
| CDS exons per chain median / max | 6.5 / 10 | 25 / 243 |
| predicted intron length min / median / max | 20 / 845 / 12,028 | 20 / 1,336 / 52,199 |

Reading: the step-900 checkpoint does not produce gene-by-gene structure.
On chr V it emits 287 chains of median 81 kb that run through many
reference genes (238 fusions, 121 splits) joined by long non-canonical
"introns" (median 1.3 kb; 94 % non-GT-AG), which is why nucleotide
sensitivity reaches 0.47 while no transcript, and only 67 of 28,569 CDS
exons, are exact. On the intronless-dominated yeast chr I it emits almost
nothing (6 chains, 4 % nucleotide sensitivity). The interim step-300
checkpoint (`../score-dry/`) over-predicted instead (662 chains on chr I,
0.96 sensitivity / 0.46 precision). Working hypothesis (engels-0099's
draw replay; not established as the sole cause): both are consistent with
a model whose intron-versus-intergenic emission is barely supervised. In
the fit's training windows every non-CDS base is intron except the 10-base
flanks and the sampled background tiles; the seeded draws through step 900
contained 168 background draws (344,064 bases) plus 140,640 flank bases =
484,704 intergenic of 15,628,525 sampled bases (3.1 %; 785,532 of
26,259,772 = 3.0 % over all 1,500 steps), against 8.76 Mb CDS and 6.39 Mb
intron. (The 548 loaded tiles are the train-plus-dev inventory, 106 of
them dev.) The duration factor cannot stop the fusions on its own
(stalin-0102: −305 nats at 1,000 bases for phase 0 at step 300) because
emissions over a long intergenic stretch dominate it. The revision to test
first, before any B allowance (task goal, last bullet), is the section 3.6
adjacent-gene / wide-flank increment (a-pilot section 2) plus a background
share matched to the genome composition, scored the same way; no decoder
or scorer defect was found by the reviews.
