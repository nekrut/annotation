# benchmark/

The evaluation panel for the charter's gene prediction work. The design
document is [`docs/benchmark.md`](../docs/benchmark.md); this directory holds
the manifest and the scripts that produce and check it.

No genome or annotation data is committed here. `panel.tsv` records the
accession, the annotation release, and the MD5 of the annotation file, so a
download is verifiable later; the data itself is fetched into a directory
outside the repository.

| file | what it is |
|---|---|
| `panel.tsv` | 20 species: accessions, assembly and annotation metadata, intron/exon statistics, split assignment, annotation MD5 |
| `taxonomy.tsv` | NCBI Taxonomy lineage per species, cached 2026-09-09 |
| `fetch.py` | download and checksum-verify genomes/annotations from the NCBI FTP mirror |
| `annotation_stats.py` | recompute every statistic in `panel.tsv` from a GFF3 |
| `leakage_check.py` | enforce the held-out phylogenetic distance rule; check declared informants |
| `score.py` | score one predicted GFF3 against the reference; emits the `docs/benchmark.md` section 4 metrics as JSON |
| `degrade.py` | make a seeded, known-perturbation copy of a reference: the control run that an identity run cannot be |
| `report.py` | join scored species into the section 4.8 table and its two aggregates |
| `validation/` | scored runs behind `docs/benchmark.md` section 6 |

Python 3.11, standard library only. Nothing to install.

```
# what would be downloaded
python3 benchmark/fetch.py --all --what gff --dry-run

# one species, ~2 MB, into a scratch directory
python3 benchmark/fetch.py --species Saccharomyces_cerevisiae --what gff --dest /tmp/panel

# the statistics behind panel.tsv
python3 benchmark/annotation_stats.py /tmp/panel/Saccharomyces_cerevisiae/*_genomic.gff.gz

# the held-out rules
python3 benchmark/leakage_check.py

# the scorer against its fixtures
python3 benchmark/score.py --self-test

# score a prediction (the declaration is docs/benchmark.md section 3.3)
python3 benchmark/score.py \
    --reference /tmp/panel/Saccharomyces_cerevisiae/*_genomic.gff.gz \
    --prediction predicted.gff3 --species Saccharomyces_cerevisiae \
    --declaration my-run.yaml --genome /tmp/panel/*/*_genomic.fna.gz \
    --out results/Saccharomyces_cerevisiae.json

# the section 4.8 table
python3 benchmark/report.py results/*.json --markdown

# the control run: degrade a reference by a known amount and score it
python3 benchmark/degrade.py --self-test
python3 benchmark/degrade.py \
    --reference /tmp/panel/Saccharomyces_cerevisiae/*_genomic.gff.gz \
    --out /tmp/sc.degraded.gff3 --summary /tmp/sc.degrade.json
```

`score.py --genome` is optional and only adds the splice dinucleotide and
local-GC stratifications; everything else is computed from annotation alone.

`score.py --stop-outside-cds` is needed for AUGUSTUS, BRAKER, GeneMark and
SNAP output, whose CDS excludes the stop codon while every panel reference
includes it. The scorer detects the convention from the prediction's
`stop_codon` features and warns when the flag disagrees, because without the
correction every terminal exon, single-exon gene, stop codon and exact
transcript match is scored 3 bp short while the nucleotide score stays
healthy. See `validation/` and `docs/benchmark.md` section 4.5.

`degrade.py` deletes 10% of transcripts and moves the *downstream boundary in
transcription order* of 10% of CDS segments 3 bp further downstream. Which
field that is depends on the strand -- `end` on `+`, `start` on `-` -- so the
perturbation lands on stop codons and donor sites and never on start codons or
acceptor sites. That asymmetry is the test: a scorer whose splice-site or
codon assignment is not strand-aware cannot reproduce it. See
`docs/benchmark.md` section 6.

`fetch.py --what` accepts `gff`, `fasta`, `protein`, `cds`. The whole panel is
318 MB of gzipped GFF; the FASTAs are considerably larger.
