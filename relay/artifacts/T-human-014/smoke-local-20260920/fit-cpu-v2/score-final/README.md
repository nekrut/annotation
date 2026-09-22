# score-final/: development-chromosome scores under the finished fit-cpu-v2 `best.pt` (step 900)

Run 2026-09-22T07:27Z on core 2, one thread, code `4d766b0` (PR #38). The
`model/a` + `model/grammar` source digest (`56e472d1…`, `measure_*.json`
`source`) differs from the fit's `61a9480` digest (`429fea43…`) by the
docstring-only edits of `4d766b0` (engels-0101 and stalin-0105 verified
identical ASTs after stripping docstrings). `run_final.sh` is v1's with
the v2 checkpoint, config and output paths: per chromosome `measure
--profile chromosome --window 12288 --overlap 4096 --segments 19 --dtype
float64 --checkpoint best.pt --gff-out` under `/usr/bin/time -v`, then
`benchmark/score.py --seqids --genome --declaration`. All four exits 0
(`run_final.out`). `sha256.txt` pins the two predicted GFF3s and the
host-local `best.pt` (`99f773ed…2a8cf6`, the fit's `best_pt.sha256`). The
chr V GFF3 is committed gzipped (3.36 MB raw, 13,410 chains).
`summarize.py DIR` reproduces every number below from the JSON and GFF3
records (run it on `../../fit-cpu-v1/score-final` for the v1 column).

## Cost (checkpoint-independent apart from the decoded chain count)

| chromosome | bases | tiles | CPU-s | CPU-s/Mb | encoder | decode | preprocess | io | output | peak RSS | process wall |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *S. cerevisiae* chr I | 230,218 | 76 | 2.89 | **12.54** | 1.32 | 1.53 | 0.034 | 0.008 | 0.001 | 0.96 GiB | 3.5 s |
| *C. elegans* chr V | 20,924,180 | 3,420 | 185.8 | **8.88** | 76.9 | 106.1 | 2.26 | 0.44 | 0.054 | 2.05 GiB | 186.6 s |

v1 rows: 12.09 and 8.86 CPU-s/Mb. `score.py`: chr I 0.14 CPU-s / 32 MB;
chr V 4.36 CPU-s / 175 MB (13,410 predicted transcripts against 287).
`/usr/bin/time` user + system: chr I 3.49 s, chr V 186.4 s. The cost is
unchanged by the checkpoint (decode 1.53 / 106.1 s traces back 301 /
13,460 chains against v1's 45 / 324 at the same cost, so traceback is
not the decode's cost). The CPU verdict of `../../fit-cpu-v1/score-final/`
stands: both rows exceed this machine's portable ceiling of 4.67 CPU-s/Mb
(2.7× and 1.9×; machine-normalized 40.3 and 28.5 against 15); memory
passes.

## Accuracy (`score_*.json`; development chromosomes of seen species, declaration in `declaration.yaml`; v1 in parentheses)

| metric | chr I (94 ref transcripts, 3 ref introns) | chr V (6,766 ref transcripts / 4,995 loci, 22,824 ref introns) |
|---|---|---|
| predicted transcripts | **214** (6) | **13,410** (287) |
| nucleotide sens / prec / F1 / MCC | **0.931 / 0.742 / 0.825 / 0.746** (0.039 / 0.356 / 0.071 / 0.020) | **0.531 / 0.546 / 0.538 / 0.468** (0.470 / 0.586 / 0.522 / 0.461) |
| predicted / reference CDS bp | 177,189 / 141,204 | 5,448,717 / 5,601,813 |
| locus TP / FP / FN, fusion, split | 86 / 128 / 8, 2, 1 (3 / 3 / 91, 2, 0) | 3,804 / 9,606 / 1,191, **4**, **1,827** (269 / 18 / 4,726, 238, 121) |
| locus sens / prec | 0.915 / 0.402 | 0.762 / 0.284 |
| transcript exact (sens / prec) | **69** (0.734 / 0.322) (0) | **101** (0.020 / 0.008) (0) |
| exact CDS exons TP / FP / FN | 69 / 173 / 28 (0 / 37 / 97) | 239 / 15,461 / 28,330 (67 / 10,481 / 28,502) |
| exact exons by type TP (FP) single / initial / internal / terminal | 69 (123) / 0 (22) / 0 (6) / 0 (22) | 107 (11,428) / 37 (1,838) / 12 (403) / 57 (1,818) |
| predicted introns by `score.py` dinucleotide category | 28: 1 GT-AG / 27 other (31: 1 / 30) | 2,290: 165 GT-AG (68 TP + 97 FP) / 7 GC-AG / 4 AT-AC / **2,114 other** (10,261: 586 / 43 / 10 / 9,622) |
| donor / acceptor site TP (F1) | 0 / 0 | 168 (0.014) / 306 (0.025) (613 / 751) |
| start / stop codon TP (FP) | 77 (137) / 79 (135) (1 / 0) | 1,412 (11,996) / 1,377 (12,029) (7 / 3) |
| predicted chain span min / median / max | 117 / 347 / 4,416 (234 / 11,212 / 15,989) | 96 / **294** / 12,343 (147 / 81,233 / 842,820) |
| chains by CDS exon count 1 / 2 / 3 / 4+ | 197 / 12 / 4 / 1 | 11,535 / 1,561 / 242 / 72 |
| predicted intron length min / median / max | 20 / 21 / 403 | 20 / **34** / 6,785 (20 / 1,336 / 52,199) |
| predicted intron lengths 20–30 / 31–50 / 51–100 / 101–500 / >500 | 25 / 1 / 0 / 2 / 0 | 1,078 / 486 / 425 / 277 / 24 |
| predicted chains with CDS < 300 bases | 92 of 214 | 6,827 of 13,410 |

Reading, in two parts.

**The fusions are sharply reduced.** The v1 fit's whole-chromosome fault —
287 chains of median 81 kb on chr V running through many genes (238
fusions), and almost nothing on chr I — is nearly absent under the v2
checkpoint. Chr V chains have median span 294 bases and 4 fusions (chr I
2); chr I now has 86 of 94
reference loci hit (nucleotide sensitivity 0.93, MCC 0.75), 69 reference
transcripts of the 91 single-exon ones reproduced exactly with their
start and stop codons, and 77 / 79 exact start / stop codons against
1 / 0. This is a comparison of two loader configurations, not a
single-variable experiment: the v2 intervention (`context` 512 plus 800
background draws per species; `../README.md`) raises the sampled `U`
share from 3.1 % to 30.6 % of the bases seen through the selected step
900 (stalin-0106 replay: 484,704 of 15,628,525 against 6,486,097 of
21,177,237 bases), but it also changes the context geometry of every
window, the background pool, which chains are drawn, the length
exclusions, and the development windows on which the checkpoint was
selected (stalin-0104). The common-chromosome result supports the v1
working hypothesis — that the fusions came from a model whose
intron-versus-`U` emission was barely supervised — more strongly than
the draw replay alone did; it does not isolate the `U` fraction as the
sole cause, and it does not show that the supervision defect is fully
removed (4 + 2 fusions remain).

**Candidate A still misses the accuracy target in this fit, on splicing
and on precision.** The v2 checkpoint predicts almost no real intron:
on chr V it emits 2,290 introns of median 34 bases (1,078 at 20–30 bases,
`min_intron` is 20) of which 2,125 are non-GT-AG (2,114 `other` + 7 GC-AG
+ 4 AT-AC) and 68 match a reference intron, giving donor / acceptor F1 0.014 / 0.025 against 22,824 reference
introns; 11,535 of 13,410 chains are single-exon. Multi-exon genes are
therefore emitted as several single-exon chains (1,827 splits, 239 exact
CDS exons of 28,569), so transcript-level sensitivity on chr V is 0.020
even though locus sensitivity is 0.76. Precision is the second miss: `score.py`
reports 9,606 of 13,410 chr V predicted loci and 128 of 214 chr I loci
as false positives, and about half of the predicted chains carry under
300 CDS bases. A locus FP in `score.py:loci` is a prediction left
unmatched after one-to-one greedy matching on shared same-strand CDS
bases, so it counts two different things; splitting the FP by whether
the predicted locus shares any same-strand reference CDS base (read-only
replay of the scorer's matching on the checksummed inputs, stalin-0106,
reproduced here):

| v2 chromosome | matched loci | FP with no same-strand reference CDS overlap | FP overlapping reference CDS but unmatched (extra fragments of a matched gene) | reported locus FP |
|---|---:|---:|---:|---:|
| *S. cerevisiae* chr I | 86 | 127 | 1 | 128 |
| *C. elegans* chr V | 3,804 | 6,376 | 3,230 | 9,606 |

On chr V a third of the FP loci are fragments of genes the model already
hits (the same fragmentation the 1,827 splits measure), and two thirds
share no same-strand CDS base with any reference gene; that subset is
the candidate set for short open reading frames the emissions accept as
coding once the `U` state is available, but it has not been checked
against full gene spans (UTRs) or opposite-strand genes, so "no reference
gene" is not established for it. The tiny predicted introns are
consistent with one reading — the decoder using an intron of minimal
length as the cheapest way to bridge a frame break inside what it scores
as coding, rather than a learned splice event — which is a hypothesis
from the length distribution, not a measurement of the mechanism; the motif tables that moved in the canonical direction
during v1 (stalin-0102) have not become decisive splice-site emissions in
either bounded fit (1,500 steps of batch 8, ~0.6 of one pass, both fits
best at step 900). Whether that is the training budget or the encoder's
splice-site receptive field is not decided by this run.

Consequences for the task goal (last bullet): the CPU regime still misses
as in a-pilot 3.2 and 3.3, so no positive CPU allowance for B; the
accuracy failure has moved from a supervision-balance defect (sharply
reduced under the loader revision) to a splice-site / precision defect. The revision to
test next, before any B allowance, is a longer fit under the same loader
(the GPU grant of lenin-0083 exists for exactly this), scored the same
way; a decoder change (an intron duration floor above `min_intron`, or a
short-ORF penalty) is not proposed until a fit that covers more than one
pass shows whether the splice-site emissions sharpen on their own.
