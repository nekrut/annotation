# scripts/data

Data-access tools for the gene prediction project (task T-human-008). The
inventory of sources, sizes, licences and rate limits is in
`docs/data-sources.md`. Standard library only; Python 3.11.

## fetch_window.py

Fetch one locus window with its multiple alignment, tree, reference
sequence, annotation and conservation scores. Run `--help` for the
argument list; the module docstring describes the two backends (UCSC API
plus HTTP range reads into the uncompressed MAF, and Ensembl Compara REST)
and the output files. Nothing is written inside the repository; point
`--out` at a scratch directory.

Demonstration runs on 2026-09-09 from a cloud runner with no cached data
(the definition of done asks for one human and one non-mammal locus; the
other rows exercise the remaining code paths). Sizes are what the runner
downloaded; the SHA-256 prefix is of the `.maf` written, so a rerun can be
compared. Alignments change when the sources update, so a differing hash
is a prompt to look, not an error.

| Run | Assembly | Window (0-based, half-open, flank 500 included) | Alignment | Blocks | Sources | Transcripts | Download / requests | MAF bytes | MAF sha256[:16] |
|---|---|---|---|---|---|---|---|---|---|
| HBB_470 | hg38 | chr11:5224964-5229895 | multiz470way (bigMaf via API) | 241 | 461 | 1 | 15.8 MB / 7 | 15,359,758 | 7a85a363d5b240c7 |
| HBB_241 | hg38, locus given as `NC_000011.10` | chr11:5224964-5229895 | cactus241wayBM (bigMaf via API; RefSeq name resolved through chromAlias) | 2517 | 240 | 1 | 68.6 MB / 9 | 67,775,774 | 0753bfed06589869 |
| Adh_124 | dm6 | chr2L:14615052-14619402 | multiz124way (wigMaf index + Range reads) | 352 | 89 | 12 | 4.6 MB / 11 | 3,398,732 | e4f6282ff8e9043b |
| TP53_100 | hg38 | chr17:7667921-7687990 | multiz100way (wigMaf index + Range reads) | 1021 | 100 | 28 | 9.1 MB / 14 | 7,252,181 | f6a925b3be471e0d |
| ATP5PO_kg | hg38 | chr21:33901026-33915853 | multiz100way (wigMaf index + Range reads), `--annotation knownGene` | 998 | 100 | 53 | 9.8 MB / 15 | 7,975,753 | c65620e6d9224180 |
| GAPDH_sauropsids | gallus_gallus (GRCg7b) | 1:76900098-76906736 | Ensembl EPO sauropsids, with 5 ancestral rows | 2 | 11 | 5 | 0.1 MB / 9 | 54,580 | 6e24cc945abb6fd8 |
| Gapdh_mammals | mus_musculus (GRCm39) | 6:125134789-125144018 | Ensembl EPO 44 mammals, mouse as query species | 1 | 3 | 65 | 0.2 MB / 6 | 32,258 | 3dc17314f67d2554 |
| Pf_none | GCF_000002765.6 (GenArk hub) | NC_004325.2:100000-120000 | `--track none`: RefSeq annotation and sequence only | 0 | 0 | 4 | 0.0 MB / 3 | 90 | 2d35a98bb3d5e6a8 |

The commands, in the same order:

```
python3 scripts/data/fetch_window.py --assembly hg38 --locus chr11:5225464-5229395 --flank 500 --track multiz470way --out /tmp/win/HBB --name HBB_470
python3 scripts/data/fetch_window.py --assembly hg38 --locus NC_000011.10:5225464-5229395 --flank 500 --track cactus241wayBM --out /tmp/win/HBB241 --name HBB_241
python3 scripts/data/fetch_window.py --assembly dm6 --locus chr2L:14615552-14618902 --flank 500 --track multiz124way --out /tmp/win/Adh --name Adh_124
python3 scripts/data/fetch_window.py --source ensembl --assembly gallus_gallus --locus 1:76900598-76906236 --flank 500 --species-set sauropsids --out /tmp/win/GAPDH --name GAPDH_sauropsids
python3 scripts/data/fetch_window.py --source ensembl --assembly mus_musculus --locus 6:125135289-125143518 --flank 500 --species-set mammals --out /tmp/win/mm --name Gapdh_mammals
python3 scripts/data/fetch_window.py --assembly GCF_000002765.6 --locus NC_004325.2:100000-120000 --track none --out /tmp/win/Pf --name Pf_none
```

Notes from the runs:

- Every source in the fetched blocks was a leaf of the tree fetched next
  to the track (461 of 461 for the 470-way, 240 of 240 for the 241-way
  Cactus, 89 of 89 for the fly 124-way).
- The fly window took 358 requests when each block was fetched by its own
  range request; coalescing the contiguous index offsets into 1 MiB reads
  brought it to 11 requests, which is the current behaviour.
- Ensembl's `mammals` EPO block at mouse `Gapdh` holds only mouse, rat and
  their ancestor; the same locus under `murinae` returns 43 rows. Which
  species set answers depends on how Compara partitioned the region, so a
  training pipeline on Ensembl should query the largest set that contains
  the reference and fall back.
- Mouse `Hbb-bs` (7:103475730-103482989) has no block at all in the
  `mammals` EPO set and 43 plus 15 rows in `murinae`; duplicated loci are
  where EPO coverage is patchiest.
- The per-species coverage tables (`coverage` in each manifest) are the
  source of the numbers in `docs/data-sources.md` section 6.

## coverage_by_distance.py

Aggregates the per-informant coverage tables from many `fetch_window.py`
manifests, binned by patristic distance from the reference on the track's
tree (`--bins` gives the upper edges in substitutions per site). Reports
base-weighted coverage per annotation class and the share of
informant-window pairs above 0.5; `--tsv` writes every pair. The numbers in
`docs/data-sources.md` section 6.1 come from:

```
python3 scripts/data/coverage_by_distance.py /tmp/cov/fly --reference dm6 --markdown --bins 0.1,0.25,0.5,1.0,2.0
python3 scripts/data/coverage_by_distance.py /tmp/cov/human --reference hg38 --markdown --bins 0.1,0.25,0.5,1.0,2.0
python3 scripts/data/coverage_by_distance.py /tmp/cov/mouse --reference mm39 --markdown --bins 0.1,0.25,0.5,1.0,2.0
python3 scripts/data/coverage_by_distance.py /tmp/cov/worm --reference ce11 --markdown --bins 0.1,0.25,0.5,1.0,2.0
```

A source that is not a leaf of the tree is reported on stderr and skipped
for the bins, but a tree leaf with no coverage entry counts as an
all-unaligned informant; that is how the mm39 35-way's GenArk assembly
(`GCF_003668045.3` in the MAF, `GCF_003668045v3` in the tree) scored zero
until `fetch_window.py` learned the naming convention (`maf_source()`), so
check the warning before trusting a bin.

## cut_windows.py

Turns one fetched locus into fixed-length training examples (`.npz` plus a
JSON sidecar) following the convention in `docs/data-sources.md` section
6.3: reference one-hot codes, per-base labels, CDS frame, strand, boundary
marks, one row per informant with explicit `gap` and `unaligned` states, an
insertion-length channel, tree distances, conservation. Standard library
only; the `.npz` is written by hand and reads with `numpy.load`.
`--drop-species` removes held-out informants (the benchmark leakage rule)
and, on Ensembl EPO, the ancestral rows of any clade containing one;
`--drop-ancestors` removes all ancestral rows; `--transcript-types`
chooses which transcripts paint labels (default: the benchmark's truth
rule, no pseudogenes or Ig/TCR segments); `--isoforms` chooses which
isoforms of a locus paint (`union`, the default; `longest-cds`, one per
locus; or `representative` with `--representatives FILE`, for example a
MANE Select list, falling back to `longest-cds` where a locus has none
listed; the sidecar records the policy and every isoform dropped, with
its locus and reason, and under `distinct_sites` the distinct start
codons, stop codons, donors and acceptors over all isoforms against those
the painted isoforms keep, so the policy's cost in splice sites is on
record); `--partial-ends` decides whether a CDS end is complete and so
gets a start or stop mark (`both`, the default: the annotation's
declaration where the fetcher recorded one, from RefSeq's `cdsStartStat`
columns, GENCODE's `cds_start_NF` / `cds_end_NF` tags or a non-zero first
coding-exon frame, else the reference sequence, ATG first and a stop last
under `--genetic-code` or `--stop-codons`; `declared`; `sequence`; or
`none`, the pre-0.5 behaviour of marking every end); the sidecar's
`cds_ends` records the declared and sequence evidence, the transcripts
judged incomplete at each end, a declared frame offset, and every
disagreement between the two signals; `--reference-anchored` declares
a track unknown to the built-in class table to be reference-anchored;
`--duplicate-rows` chooses the copy when one species has several rows in
one block, as Cactus exports (`identity`, the default: the copy with most
bases identical to the reference in that block; or `first`), and the
sidecar's `block_selection` counts the discarded copies, the positions
covered by overlapping blocks with what first-wins lost there, and the
minus-strand reference blocks flipped into forward coordinates (section
6.3); `--both-strands` adds the reverse-complement example. `--self-test` checks
labels, frames, boundaries, informant codes, insertions, distances, the
reverse complement, the transcript filter, the alignment-class table,
label counts on a padded window, an EPO-shaped two-block window with
ancestral rows, a gene nested on the opposite strand inside another gene's
intron (strand and frame follow the nested gene), byte-identical
archives from repeated runs, and the three isoform policies on a locus
with two isoforms plus an unnamed locus clustered by span, the
distinct-site accounting under each policy, and the incomplete-end rules
on a fixture with real codons (declared, undeclared, declared against the
sequence, a declared frame on the minus strand, the four modes, genetic
code 6, a stop-excluded convention and a CDS reaching past the window),
overlapping blocks with a lost informant base, a minus-strand reference
row and two copies of one species under both policies (58 checks). The `.npz` members carry a fixed
timestamp, so an example's checksum depends only on its inputs; the whole
cutter, the fetcher, `coverage_by_distance.py` and `tree_composition.py`
were checked to give identical output under `PYTHONHASHSEED` 0, 1, 2 and 42
on the fly `Adh` window (a stratified count that depends on set iteration
order is the defect lenin found in the benchmark scorer, relay note
20260909T133356Z-lenin-0013).

```
python3 scripts/data/cut_windows.py --self-test
python3 scripts/data/cut_windows.py --stem /tmp/win/Adh/Adh_124 --out /tmp/ex --length 2048 --stride 1024 --drop-species apiMel4 --both-strands
python3 scripts/data/cut_windows.py --stem /tmp/win/GAPDH/GAPDH_sauropsids --out /tmp/ex --length 0 --both-strands --drop-species taeniopygia_guttata
python3 scripts/data/fetch_window.py --assembly hg38 --locus chr21:33901026-33915853 --track multiz100way --annotation knownGene --out /tmp/win/ATP5PO --name ATP5PO_kg
python3 scripts/data/cut_windows.py --stem /tmp/win/ATP5PO/ATP5PO_kg --out /tmp/ex --length 0
```

The last two commands are the incomplete-end demonstration (section 6.3
of the document): a GENCODE window with four `cds_start_NF` and one
`cds_end_NF` transcript among 45 coding ones, where the sequence rule
agrees with every tag, the default cut paints one start mark against
five under `--partial-ends none`, and the tagged 3' end lies outside the
window so only the declaration catches it.

The third command on the chicken `GAPDH` window (re-fetched 2026-09-09
after `fetch_window.py` started recording the species set; 14 requests,
92,661 bytes): the 17-member sauropsid EPO set gives 16 leaf rows in sorted
order minus the reference, the manifest lists 6 ancestral rows, and
dropping zebra finch removes its row and the 5 ancestors whose clade
contains it (`Ggal-Mgal-Pmaj-Scan-Tgut[5]`, `Ggal-Pmaj-Tgut[3]`,
`Pmaj-Scan-Tgut[3]`, `Pmaj-Tgut[2]`, `Scan-Tgut[2]`), leaving `K` = 16 with
`Ggal-Mgal[2]` as the one ancestor, `alignment_class: jointly_inferred`
from the table, `dropped_rows_only: true`, and identical `label_counts`
(1,275 CDS, 387 UTR, 4,334 intron, 642 intergenic) on the forward and
reverse-complement sidecars. Tree distances are 0.038 to the chicken-turkey
ancestor, 0.079 to turkey, 0.168 to great tit and 0.175 to canary; the
other members do not align in this window and are all-unaligned rows.

## sample_genes.py and the sampling runs (section 6.1 of the document)

`docs/data-sources.md` section 6.1 aggregates 12 fly, 10 human, 10 mouse
and 10 worm windows. The fly and human genes were drawn with a fixed seed
from `ncbiRefSeqCurated` on dm6 chr2L and hg38 chr11 (NM_ transcripts,
complete CDS, at least 3 exons, transcript length 2 to 8 kb for fly and
3 to 12 kb for human, one per gene, non-overlapping), so the human sample
is biased to short genes; the list is in the document. That draw was done
in an interactive session; `sample_genes.py` is the same procedure as a
script (longest qualifying NM_ per gene, genes overlapping any other
transcript of the track dropped, `random.Random(seed).shuffle`, first `n`),
and the mouse and worm samples were drawn with it on 2026-09-09:

```
python3 scripts/data/sample_genes.py --genome mm39 --chrom chr19 --length 3000-12000 --n 10 --seed 20260909
python3 scripts/data/sample_genes.py --genome ce11 --chrom chrIII --length 2000-8000 --n 10 --seed 20260909
```

Each printed line is `gene`, `transcript`, `chrom:start-end`, `exons`; the
locus goes straight to `fetch_window.py --locus ... --flank 500`
(`--track multiz35way` on mm39, `--track multiz135way` on ce11,
`--conservation none` since the coverage run does not use it). The header
line records the track's `dataTime`, which is what the draw depends on. Midway through the run `hgdownload.soe.ucsc.edu` reset
every connection for several minutes while the API host and
`hgdownload2.soe.ucsc.edu` kept serving; the fetcher now rotates through
the mirrors on connection errors and `--download-host` sets the primary,
and the remaining windows were fetched with
`--download-host https://hgdownload2.soe.ucsc.edu`. Requests that fell back
carry `mirror_for` in the manifest log.

## tree_composition.py

Summarises a track's Newick tree: leaves and distinct species per clade,
patristic distance from the reference (min, median, max), and how many
leaves fall within 0.5 and beyond 1.0 substitutions per site, the
non-coding alignability horizon of `docs/data-sources.md` section 6.2.
`scripts/data/hg38.470way.orders.tsv` maps the 470 leaves of
`hg38.470way.scientificNames.nh` to their order and family (GBIF
species-match API, 2026-09-09; one leaf, the Hawaiian monk seal, filled by
hand after a connection reset). The composition table in section 2.1 of the
document is:

```
curl -sSO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/hg38.470way.scientificNames.nh
python3 scripts/data/tree_composition.py hg38.470way.scientificNames.nh --reference Homo_sapiens \
    --clades scripts/data/hg38.470way.orders.tsv --markdown
```

`--lookup-gbif <tsv>` rebuilds the clade table from GBIF (about ten
minutes for 470 leaves at one request per leaf). Without `--clades` the
script prints one row and the distance histogram, which on
`hg38.100way.nh` with `--reference hg38` reproduces the informant counts
of the section 6.1 human table.

## pairwise_codons.py

Exports reference-versus-informant codon pairs from a fetched window for
the KA/KS baseline (T-human-010), one sequential PHYLIP per informant,
with a `.pairs.tsv` and sidecar recording per informant how many codons the
CDS has, how many the pairwise rule keeps, how many a complete-case rule
over all rows would keep, the drop reasons (unaligned, gap, ambiguous,
insertion inside a codon, informant stop), identity and tree distance.
Written because codeml's pairwise mode deletes every column any input row
gaps (relay note 20260909T164106Z-stalin-0016); `docs/data-sources.md`
section 6.4 has the measurement. Reuses the MAF, tree and isoform code of
`cut_windows.py`, so it needs that file next to it.

```
python3 scripts/data/pairwise_codons.py --self-test
python3 scripts/data/pairwise_codons.py --stem /tmp/win/Adh/Adh_124 --out /tmp/pairs --isoforms longest-cds --drop-species apiMel4 --write-multi
python3 scripts/data/pairwise_codons.py --stem /tmp/win/Adh/Adh_124 --out /tmp/cc --transcript NM_001032098.2 --informants droSim2,droSec1,droEre2,droYak3 --complete-case
```

On the Adh window (2026-09-09) the first command writes 159 files (two
coding loci, Adh and Adhr; the enclosing outspread transcript has CDS
outside the window and is skipped with a message), byte-identical under
`PYTHONHASHSEED` 0 and 42. Self-test: 17 checks on a synthetic window built
from stalin's three-row example plus a minus-strand spliced gene, an
insertion inside a codon, an informant stop and a tree leaf absent from
every block.

## Politeness

The fetcher spaces requests at least 0.34 s apart by default (`--pause`),
retries 429 and 5xx with backoff, sends a descriptive User-Agent, and logs
every request with its byte count into the manifest. UCSC asks for about
one request per second on the API and NCBI for at most three per second
on E-utilities; Ensembl's limit was 55,000 requests per hour on the day of
the runs.
