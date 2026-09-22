# score-pwm/: the fit-cpu-v3 checkpoint decoded with the splice mask **and** the splice-site PWM, 2026-09-22T19:4xZ

The section 3.5 **decoder** increment 2, the revision section 3.4 asked
for. The mask of `../score-canonical/` decides *whether* a junction is
in the support; nothing in the decoder decided *which* of the surviving
sites to prefer, because after the mask the learned 16-entry
dinucleotide table scores every legal donor identically. That is
visible in `../score-canonical/`: masking doubled true GT-AG introns
(433 → 1,070) but quadrupled the false ones (470 → 1,970, precision
0.480 → 0.352).

`model/a/splicepwm.py` supplies the missing ranking as a position-weight
matrix: donor offsets `[-3, +6)` around the junction (three exon bases
and six intron bases) and acceptor offsets `[-20, +3)` (twenty intron
bases, long enough to hold a polypyrimidine tract, and three exon
bases), scored as natural-log odds against the base composition of the
same windows, with a pseudocount of 1. It is added to the `donor` and
`acceptor` emission rows exactly as the dinucleotide bias is, so both
kernels and both reference decoders are unchanged, and it is
**inference only**: the chain-loss numerator is untouched, so a
reference intron with an atypical junction stays reachable and the
encoder is not fitted against a bias estimated from its own labels.

**No refitting.** The weights are the same `best.pt` (step 3,000 of
3,000, dev NLL 28.640) that produced `../score-final/` and
`../score-canonical/`; the only differences from `../score-final/` are
`--canonical-splice` and `--splice-pwm`.

## The matrix and its leakage rules

`../pwm-v1/splicepwm.json` (sha256 `a964e3e6…`), fitted by
`model.a.train pwm`, which loads the **fit's own config** and then calls
`split_windows` with every source's `dev_seqids`, so the two development
chromosomes this run is scored on — S. cerevisiae `NC_001133.9` and
C. elegans `NC_003283.11` — contribute no junction and no background
base. It counts **67,325 donors and 67,325 acceptors over 19,975 train
windows**, with 4,951 development windows excluded; those two window
counts are exactly the `train_windows` and `dev_windows` of
`../run_manifest.json`, so the matrix was estimated on precisely the set
the encoder was fitted on. The JSON records the spans, the site counts,
and per source the excluded dev seqids with the pinned GFF3/FASTA MD5s.
`run_pwm.sh` re-checks the matrix's sha256 before scoring, so a rerun
cannot silently score a refitted matrix.

The estimated consensus is the textbook one and was not put there by
hand: donor `MAG|GTAAGT` (`A` at −2, `G` at −1, `GT` at 0/+1 at 1.99 and
1.94 bits, then `AAGT`), acceptor a T-rich tract at −8…−4 followed by
`CAG|` (`AG` at −2/−1 at 1.99 bits each). Both species pool into one
matrix; see the caveat below.

## Reproduction

`run_pwm.sh` (at `../run_pwm.sh`, shared with `../score-pwm-only/`) is
`../score-final/run_final.sh` with the two flags added, the output
directory selected by its `$1`, and the PWM digest gate appended; the
completion gate, the failure propagation and every other workload
argument are unchanged. All four workloads exit 0 (`run_pwm.out`), and
`summarize.py` (copied unchanged) reproduces every number below from the
committed JSON and GFF3.

Provenance. The measurement rows record `commit 72eecf5` with
`source_dirty: true` and `source_sha256
e65ea9b181a87fde436faa75a60d9064d94d77139e1afbabbd3d05e777345423`,
because the run happened before the code was committed. That digest is
byte-identical to the one `_source_provenance()` reports on the clean
tree at **`c823d97`**, so the executed source is exactly what `c823d97`
contains — the dirty flag is honest, not a gap.

## Cost (one core, `taskset -c 2`, `/usr/bin/time -v`, float64, 19 segments)

| chromosome | CPU-s (user+sys) | CPU-s/Mb mask+PWM | mask only | unmasked v3 | decode | peak RSS |
|---|---|---|---|---|---|---|
| S. cerevisiae chr I | 2.87 | **12.47** | 12.14 | 12.38 | 1.55 | 0.96 GiB |
| C. elegans chr V | 187.92 | **8.98** | 8.65 | 8.88 | 106.85 | 2.05 GiB |

The PWM costs **+3.88%** end to end on chr V against the masked decode
(8.6454 → 8.9810 CPU-s/Mb) and **+2.76%** on chr I; against the
unmasked default it is +1.14% on chr V and +0.74% on chr I. These are
single-run observations on one pinned core, not repeated timings. The
CPU verdict is unchanged: every configuration stays far under the
accepted **15 CPU-s/Mb** budget, and **no positive CPU allowance for B
follows** from this run.

## Accuracy, C. elegans chr V (`NC_003283.11`, 20.92 Mb, 6,766 reference transcripts)

| metric | unmasked v3 | mask | **mask+PWM** |
|---|---|---|---|
| nucleotide F1 | 0.5978 | 0.5847 | **0.7675** |
| nucleotide sensitivity | 0.5334 | 0.5033 | **0.9477** |
| nucleotide precision | 0.6799 | 0.6976 | 0.6450 |
| exon exact F1 (all) | 0.0300 | 0.0622 | **0.5216** |
| exon exact F1 (internal) | 0.0152 | 0.0485 | **0.6055** |
| donor F1 | 0.0705 | 0.1173 | **0.6420** |
| acceptor F1 | 0.0812 | 0.1577 | **0.6662** |
| acceptor precision | 0.2247 | 0.5626 | 0.5253 |
| correct GT-AG introns (TP/FP) | 433/470 | 1,070/1,970 | **18,195**/20,131 |
| exact transcripts | 70 | 101 | **835** |
| transcript F1 | 0.0126 | 0.0180 | **0.1334** |
| locus F1 | 0.6029 | 0.5900 | **0.6386** |
| locus fusion / split | 312/910 | 226/964 | 735/229 |

This is the largest single change any decoder increment has produced.
Splice placement is no longer the dominant failure on chr V: donor and
acceptor F1 go from 0.12/0.16 to 0.64/0.67, correct GT-AG introns from
1,070 to 18,195 of the 22,620 reference GT-AG introns, and exact
transcripts from 101 to 835. Sensitivity, not precision, is what moved:
nucleotide sensitivity 0.503 → 0.948 at a precision cost of 0.698 →
0.645. The fusion/split trade also inverts — splits fall 964 → 229
while fusions rise 226 → 735 — so the decoder now over-joins where it
used to over-split.

**Intron duration also improves but is still not calibrated.** Under the
scorer's convention the predicted chr V introns go from 3,642 with a
median of 550 b to 38,730 with a median of **160 b**, against a
reference chr V q50 of 56 b; the fraction above the reference q90 of
695 b falls 46.02% → **17.17%**. The overshoot named in section 3.3
item 5 is reduced, not closed, and the duration law itself was not
touched by this change.

**A still misses the accuracy target.** 835 exact transcripts of 6,766
reference transcripts is a per-reference-transcript exact sensitivity of
0.123 (transcript F1 0.1334), and 4,425 reference GT-AG introns are
still missed. This run does not accept T-human-014.

## Accuracy, S. cerevisiae chr I (`NC_001133.9`) — the PWM **hurts** here

| metric | unmasked v3 | mask | **mask+PWM** |
|---|---|---|---|
| nucleotide F1 | 0.8906 | 0.8902 | **0.8569** |
| exon exact F1 (all) | 0.5818 | 0.5882 | **0.3959** |
| exact transcripts (of 94) | 64 | 65 | **55** |
| locus F1 | 0.7861 | 0.7783 | **0.7215** |
| predicted transcripts | 107 | 109 | 125 |
| GT-AG introns (TP/FP) | 0/0 | 0/12 | **2/67** |

Yeast chr I is nearly intronless, and the pooled matrix is dominated by
C. elegans: of the 67,325 junctions counted, the overwhelming majority
are nematode, so a *C. elegans* splice-site score is being applied to a
*S. cerevisiae* chromosome. It finds sites that are not there — 67 false
GT-AG introns against 2 true ones — and splits single-exon genes, which
is exactly why exact transcripts fall 65 → 55 and locus F1 0.778 →
0.722. The gain on chr V and this loss on chr I are the same mechanism
seen from two intron-density regimes.

This is a **species-independence** finding, and it is the charter's
central concern, not an incidental one: a decoder term fitted by pooling
junctions across clades is not clade-neutral. The next revision should
weight or condition the matrix per species (or per intron density)
rather than pool it, and re-measure both chromosomes. Nothing here
argues for keeping the pooled matrix as a default.

## What this run does not establish

- It is one checkpoint, two development chromosomes, two species, one
  run per configuration. No held-out species was touched.
- The columns are treated as independent, the standard weight-matrix
  assumption; it is false at a real branch point, so the score is a
  ranking device, not a likelihood.
- The mask contributes almost nothing once the PWM is present
  (see `../score-pwm-only/`): the PWM alone gives donor F1 0.6403 vs
  0.6420 with the mask on chr V. The mask is not what produced the gain.
- No fitting of the encoder, no GPU, no cluster, no held-out access.
