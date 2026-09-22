# score-pwm-bal-only/: the source-balanced PWM **without** the hard mask, 2026-09-22

The control for `../score-pwm-bal/`, exactly as `../score-pwm-only/` is
the control for `../score-pwm/`: same `best.pt`, same tiling, same
scorer, same pinned core, the same balanced matrix
(`../pwm-balanced-v1/splicepwm.json`, sha256 `067b55e3…`), and the
accepted finite-score support of proposal 3.1 instead of
`--canonical-splice`.

It answers whether section 3.5's conclusion — that the PWM, not the
section 3.4 mask, carries the increment — survives the reweighting. It
does:

| metric | masked (`../score-pwm-bal/`) | unmasked (here) |
|---|---|---|
| chr I exact transcripts (of 94) | 47 | 46 |
| chr I exon exact F1 | 0.3300 | 0.3172 |
| chr V exact transcripts (of 6,766) | 706 | 675 |
| chr V exon exact F1 | 0.5008 | 0.4923 |
| chr V correct GT-AG introns | 16,963 | 16,928 |
| chr V donor / acceptor F1 | 0.628 / 0.646 | 0.621 / 0.638 |

The mask is worth 0.0085 exon-exact F1 on chr V and 0.0129 on chr I
(one transcript), against 0.0018 on chr V with the pooled matrix
(0.52164 masked vs 0.51987 unmasked): a larger cleanup than pooling
needed, still an order of magnitude below what the PWM itself moves.
An earlier revision of this paragraph said the mask was worth nothing
on chr I's exactness, which the table above contradicts (stalin-0118).
The unmasked decode stays the accepted default.

Cost: chr V 9.08 CPU-s/Mb (user-only 8.0, user+system 9.1), chr I 12.50
(12.8 / 15.0). As in section 3.5, the unmasked configuration is the most
expensive of the four, because the mask removes transitions the scan
would otherwise consider.

Provenance as in `../score-pwm-bal/`: `commit 9ec0424`,
`source_sha256 3e433c82…`, byte-identical to the clean tree at
`3fd9032` on `work/T-human-014-lenin`. `summarize.py` regenerates every
number above from this directory.
