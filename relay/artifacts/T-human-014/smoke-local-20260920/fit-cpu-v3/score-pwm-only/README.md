# score-pwm-only/: the splice-site PWM **without** the hard mask, 2026-09-22T19:5xZ

The control for `../score-pwm/`. Same weights, same tiling, same
scorer, same pinned core, same matrix (`../pwm-v1/splicepwm.json`,
sha256 `a964e3e6…`); the only difference from `../score-pwm/` is that
`--canonical-splice` is **not** given, so the decode keeps the accepted
finite-score support of proposal 3.1 — the PWM reorders junctions, and
nothing forbids one.

It exists to answer one question: how much of the section 3.5 gain is
the mask and how much is the PWM?

## Answer: essentially all of it is the PWM

C. elegans chr V (`NC_003283.11`):

| metric | mask only | **PWM only** | mask+PWM |
|---|---|---|---|
| nucleotide F1 | 0.5847 | **0.7661** | 0.7675 |
| exon exact F1 (all) | 0.0622 | **0.5199** | 0.5216 |
| donor F1 | 0.1173 | **0.6403** | 0.6420 |
| acceptor F1 | 0.1577 | **0.6640** | 0.6662 |
| exact transcripts | 101 | **830** | 835 |
| locus F1 | 0.5900 | **0.6383** | 0.6386 |
| non-canonical (`other`) intron FP | 0 | 246 | 0 |
| CPU-s/Mb | 8.65 | 9.12 | 8.98 |

S. cerevisiae chr I (`NC_001133.9`): PWM-only and mask+PWM agree on
every accuracy metric reported in `../score-pwm/` to four decimals
(nucleotide F1 0.8569, exon exact F1 0.3959, 55 exact transcripts,
locus F1 0.7215), so the regression there is the PWM's, not the mask's.

Adding the mask on top of the PWM moves chr V donor F1 by 0.0017 and
exact transcripts by 5 of 6,766. What it does buy is the 246
non-canonical false introns the PWM alone still emits, removed at the
cost of the reference's own unusual junctions staying unreachable — the
restricted-support trade stalin-0114 described, now at a much smaller
scale than in `../score-canonical/`.

**Consequence for the default.** The mask was justified in section 3.4
by a large placement gain that the PWM turns out to supply on its own.
The unmasked finite-score decode remains the accepted default, and on
this evidence the PWM is the increment worth keeping; the mask is a
small, optional cleanup on top of it. Cost is the other way round: the
PWM alone is the most expensive of the four configurations on chr V
(9.1174 CPU-s/Mb vs 8.9810 with the mask), because the mask removes
transitions the scan would otherwise consider. All four stay far under
the 15 CPU-s/Mb budget and **no positive CPU allowance for B follows**.

The caveats of `../score-pwm/README.md` — one checkpoint, two
development chromosomes, single runs, independent columns, the yeast
regression, no held-out species — apply here unchanged.
