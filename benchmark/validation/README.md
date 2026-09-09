# benchmark/validation/

Scored runs of real gene predictors, plus the panel-wide degraded-copy
control in `degraded/`, kept because they are the evidence for
the claims in [`docs/benchmark.md`](../../docs/benchmark.md) §6 and because
they are what found the defects in `score.py` and `report.py` that
section records.

No genome, annotation or prediction data is here: the declarations say what
was run and the JSON is the scorer's own output. Each JSON records the
SHA-256 of the declaration beside it, so the pairing is checkable.

| declaration | result | what it is |
|---|---|---|
| `augustus-Saccharomyces_cerevisiae.yaml` | `.json` | AUGUSTUS on *S. cerevisiae* with its own parameters |
| `augustus-Schizosaccharomyces_pombe.yaml` | `.json` | AUGUSTUS on *S. pombe* with its own parameters |
| `augustus-Schizosaccharomyces_pombe-crossparam.yaml` | `.json` | AUGUSTUS on *S. pombe* with the *S. cerevisiae* parameters |
| `helixer-Takifugu_rubripes.yaml` | `.json` | Helixer 0.3.7 on *T. rubripes*, vertebrate model, GPU |
| `helixer-Neurospora_crassa.yaml` | `.json` | Helixer 0.3.7 on *N. crassa*, fungi model, GPU |
| `helixer-Saccharomyces_cerevisiae.yaml` | `.json` | Helixer 0.3.7 on *S. cerevisiae*, fungi model, GPU |
| `gencode50-Homo_sapiens.yaml` | `.json` | GENCODE 50 scored as a submission against the human RefSeq reference |
| `augustus-partial-Tetrahymena_thermophila.yaml` | `.json` | AUGUSTUS `--genemodel=partial` over 1,158 scaffolds, ids concatenated without renaming, genetic code 6 |
| `augustus-partial-Apis_mellifera.yaml` | `.json` | AUGUSTUS `--genemodel=partial` with `honeybee1`, which excludes the stop codon from the CDS |
| `tiberius-Takifugu_rubripes.yaml` | `.json` | Tiberius 2.0.7 on *T. rubripes*, `vertebrates` model, GPU — the GFF3 of the run |
| `tiberius-Takifugu_rubripes.yaml` | `-gtf.json` | the **GTF** of the same invocation, scored separately; every metric equals the GFF3 result |
| `degraded/` | the panel-wide control: all 20 references scored against a known degraded copy of themselves ([README](degraded/README.md)) |

The three Helixer runs are the second tool and the second output shape: no
`stop_codon` features at all, UTRs present, and two species (*T. rubripes*,
*S. cerevisiae*) that the shipped checkpoints' own training lists name,
against one (*N. crassa*) that they do not. The *S. cerevisiae* one is also
the direct comparison with the AUGUSTUS run on the same genome, which is what
settles §7 item 10: two unrelated tools over-predict yeast introns by 2 to
2.6x. Every JSON here was regenerated with the
current scorer, so all five carry `stop_codon_convention_from_genome` and
`stop_codon_convention_source`; the AUGUSTUS numbers are unchanged to five
decimals by that regeneration.

The two Tiberius results are one prediction in two dialects. Tiberius writes
GTF and GFF3 from one invocation (`--out tiberius.gtf tiberius.gff3`), so the
pair is a regression on the scorer itself: the two JSONs differ in exactly
five fields, all of which report what the *file* states about itself — the
input filename, `transcripts_with_stop_codon_feature`, `stop_inside_cds`,
`stop_codon_convention_detected` and `predicted_partial_source`. Every scored
metric is identical. They are kept as a pair because the GTF is what found
the `_attr` defect in §6: before the fix the scorer read that file as 241,333
single-exon transcripts and still reported nucleotide F1 0.93801 and exon F1
0.89110, the same two values the correct parse gives.
*T. rubripes* is also the one vertebrate on this panel that is genuinely held
out of a shipped checkpoint — `model_cfg/vertebrates.yaml`'s own training list
does not contain it and its test set names it — so this is the only §6
vertebrate row that is a measurement rather than an upper bound.

The third AUGUSTUS run is the ablation: the same genome and the same tool, one parameter
set away. Exon F1 falls from 0.774 to 0.296 and donor F1 from 0.854 to 0.175
while nucleotide F1 only falls from 0.955 to 0.868.

## Reproducing

AUGUSTUS is not a dependency of this benchmark and is not installed by it.
The runs above used bioconda `augustus-3.5.0-pl5321h5653ebf_10`:

```
micromamba create -p /tmp/aug -c conda-forge -c bioconda augustus
export AUGUSTUS_CONFIG_PATH=/tmp/aug/config

python3 benchmark/fetch.py --species Saccharomyces_cerevisiae --what gff --dest /tmp/panel
python3 benchmark/fetch.py --species Saccharomyces_cerevisiae --what fasta --dest /tmp/panel
zcat /tmp/panel/Saccharomyces_cerevisiae/*_genomic.fna.gz > /tmp/sc.fna

# one invocation per sequence, in parallel; AUGUSTUS restarts its gene
# numbering at g1 in every invocation, so the ids MUST be made unique before
# the parts are concatenated -- see score.py's
# predicted_conflicting_transcript_ids, which exists because of this.
mkdir parts && cd parts
awk '/^>/{n=substr($1,2); f=n".fa"} {print > f}' /tmp/sc.fna
ls *.fa | xargs -P 17 -I{} sh -c \
  '/tmp/aug/bin/augustus --species=saccharomyces_cerevisiae_S288C \
     --gff3=on --genemodel=complete --UTR=off {} > {}.gff3'
(echo '##gff-version 3'; for f in *.fa.gff3; do n=${f%.fa.gff3}; \
   grep -v '^#' $f | sed "s/ID=/ID=${n}_/; s/Parent=/Parent=${n}_/"; done) > /tmp/sc_augustus.gff3
cd ..

# AUGUSTUS ships stopCodonExcludedFromCDS=true, so --stop-outside-cds is
# required; score.py detects the convention and warns if it is omitted.
python3 benchmark/score.py \
    --reference /tmp/panel/Saccharomyces_cerevisiae/*_genomic.gff.gz \
    --prediction /tmp/sc_augustus.gff3 --species Saccharomyces_cerevisiae \
    --declaration benchmark/validation/augustus-Saccharomyces_cerevisiae.yaml \
    --genome /tmp/panel/Saccharomyces_cerevisiae/*_genomic.fna.gz \
    --stop-outside-cds --out /tmp/sc.json
```

*S. pombe* is the same with `--species=schizosaccharomyces_pombe`, and the
ablation is the same *S. pombe* FASTA with
`--species=saccharomyces_cerevisiae_S288C`.

Helixer is not a dependency either. The runs above used the published
container image and one consumer GPU; the model files are fetched once into a
mounted directory so the image stays read-only:

```
docker pull gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8
mkdir -p /tmp/helixer_models
docker run --rm -v /tmp/helixer_models:/home/helixer_user/.local/share/Helixer \
    gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8 \
    fetch_helixer_models.py --lineage vertebrate   # and --lineage fungi

zcat /tmp/panel/Takifugu_rubripes/*_genomic.fna.gz > /tmp/panel/fugu.fna
docker run --rm --gpus all -v /tmp/panel:/data \
    -v /tmp/helixer_models:/home/helixer_user/.local/share/Helixer \
    gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8 \
    bash -lc 'Helixer.py --lineage vertebrate --fasta-path /data/fugu.fna \
        --species Takifugu_rubripes --gff-output-path /data/fugu_helixer.gff3 \
        --temporary-dir /data/htmp --batch-size 8'

# No --stop-outside-cds: Helixer emits no stop_codon features and includes
# the stop in its CDS.  --genome is what lets score.py establish that rather
# than assume it (docs/benchmark.md 4.5).
python3 benchmark/score.py \
    --reference /tmp/panel/Takifugu_rubripes/*_genomic.gff.gz \
    --prediction /tmp/panel/fugu_helixer.gff3 --species Takifugu_rubripes \
    --declaration benchmark/validation/helixer-Takifugu_rubripes.yaml \
    --genome /tmp/panel/fugu.fna --out /tmp/fugu.json
```

`--batch-size 8` is not a default: the shipped default of 32 exhausts 16 GB of
GPU memory at `--subsequence-length 213840`. The image's TensorFlow 2.15.1
carries no cubins for compute capability 12.0, so on an RTX 5080 every kernel
is JIT-compiled from PTX by the driver on first use; mounting a cache
directory at `/home/helixer_user/.nv` keeps that cost to the first run.
*N. crassa* and *S. cerevisiae* are the same with `--lineage fungi`.

These JSON files were regenerated on 2026-09-09 after §4 stopped scoring
pseudogene and gene-fragment CDS rows as truth and §4.5 stopped charging
incomplete CDS ends; every F1 moved by at most 0.002, and the reference
transcript counts fell by the 6 and 32 pseudogenes the two yeasts have.

They were regenerated again the same day after §4.4 made an exact chain match
outrank shared CDS bases in the within-locus isoform pairing. Only the
`transcript` block moved, and only in the two runs against a reference with
several isoforms per locus: *T. rubripes* F1 0.23782 -> 0.25248 (tp 5,564 ->
5,907) and *N. crassa* 0.68609 -> 0.68947 (tp 6,889 -> 6,923). The three
AUGUSTUS runs and the *S. cerevisiae* Helixer run are byte-identical apart
from their regeneration timestamps.

All eleven were regenerated once more the same day when `splice.short_gaps`
gained `motif_by_class`, `motif_unresolved_windows` and `contexts_by_class`
(§4.3). The change is additive: with the three new keys removed, every
regenerated JSON is identical to the one it replaced.

`gencode50-Homo_sapiens` is not a predictor run at all: it is one human
annotation scored against another, and it is here because it is the only
submission so far with more than one isoform per locus (370,476 scored chains
against the reference's 131,442). That is what exercises the §4.4
false-positive branch on real data, and it found two defects: the §4 biotype
filter read only RefSeq's `gene_biotype` spelling, and
`predicted_transcripts_not_scored` pooled the filter's drops with off-panel
sequences. Neither touches the six runs above, which declare no biotype; the
JSON here is the output after both fixes. Reproducing it needs no tool:

```
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gff3.gz
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
# rewrite column 1 from the UCSC-style / GenBank names in the GFF3 to the
# RefSeq accessions the reference uses, using columns 5, 7 and 10 of the
# assembly report; drop rows on names the report does not map (KI270721.1,
# KI270734.1).  Nothing else about the file is changed.
python3 benchmark/fetch.py --species Homo_sapiens --what fasta --dest /tmp/panel
python3 benchmark/score.py \
    --reference /tmp/panel/Homo_sapiens/*_genomic.gff.gz \
    --prediction gencode.v50.refseqnames.gff3.gz --species Homo_sapiens \
    --declaration benchmark/validation/gencode50-Homo_sapiens.yaml \
    --genome /tmp/panel/Homo_sapiens/*_genomic.fna.gz \
    --out /tmp/gencode.json
```

`--genome` here is the scorer's largest input to date: 972,898,531 bytes of
gzip, 705 FASTA records, 3.1 Gb of sequence streamed once, 88 s and 2.12 GB
against 43 s and 1.06 GB without it. It changes no metric — it adds the §4.3 dinucleotide and
local-GC strata and a genome-read second opinion on the stop-codon
convention.

These are not benchmark results for AUGUSTUS. Both species are in AUGUSTUS's
own training set — `heldout_seen_in_pretraining: yes` in two of the three
declarations — so the first two rows are an upper bound, not a measurement.
They are here to exercise the scorer, and the third row is here because it is
the one comparison among them that the leakage rules do allow.

The two `augustus-partial-` runs are the only submissions here that declare
their own incomplete gene ends, and they were added because §4.5's
predicted-partial handling had until then only ever run on fixtures built
around the *reference* convention (`start_range=`/`end_range=`), which no
predictor writes. They are not measurements — both species are in AUGUSTUS's
training set — and they are here for three things they found:

- *T. thermophila*'s 1,158 scaffolds give 214 colliding gene ids, and merely
  *counting* them (the earlier fix) still welded their chains: 218 scored
  chains instead of 10,355 and nucleotide F1 0.020 instead of 0.420. Ids are
  now keyed by `(sequence, strand)`. The control is the same AUGUSTUS run with
  `--uniqueGeneId=true`, which scores identically block for block.
- `predicted_partial_5prime`/`_3prime` could not fire on any prediction, since
  they read the reference attribute. They now read the predictor's own
  statement — an omitted `start_codon`/`stop_codon` feature — which on
  *T. thermophila* removes 208 false-positive starts and 22 false-positive
  stops.
- *A. mellifera* is the run where `honeybee1` puts the stop outside the CDS,
  so the 3 bp extension runs and the skip for a 3'-partial chain fires on real
  data: 13 chains, 8 false-positive stops and 2 false-positive terminal exons.

All nine JSONs here were regenerated once more on 2026-09-09 after §4.3 made
its stratification units explicit, added the local-GC table for acceptors and
started counting the sites that sit in introns of more than one dinucleotide
class. **No existing field moved in any of them**: the only differences are
the new `splice.strata_units`, the new `by_local_gc.acceptor` block (the old
donor-only table is now `by_local_gc.donor`, value for value) and the two new
`sites_multiple_dinuc_classes` counts. The `gencode50-Homo_sapiens` run is
the exception, and only because it is the first one given `--genome`: it
gains the two strata and the genome-read stop-codon convention, which agrees
with the `stop_codon` features the file already carried. Its metrics are
identical to five decimals either way.

All eight JSONs here were regenerated with the scorer carrying those changes.
The three AUGUSTUS and three Helixer results are unchanged in every value
apart from the two new `codon` fields; only the GENCODE run's codon blocks
move, start F1 0.637 to 0.750 and stop F1 0.360 to 0.408, because GENCODE
marks 13,203 CDS starts and 19,535 CDS ends as incomplete by omitting the
codon feature — which its own `cds_start_NF`/`cds_end_NF` tags independently
confirm at 99.7% and 99.2% agreement.

Reproducing the two partial runs (`--uniqueGeneId=true` is what a submission
should do; the *T. thermophila* file here deliberately omits it):

```
python3 benchmark/fetch.py --species Tetrahymena_thermophila --what fasta --dest /tmp/panel
python3 benchmark/fetch.py --species Tetrahymena_thermophila --what gff   --dest /tmp/panel
zcat /tmp/panel/Tetrahymena_thermophila/*_genomic.fna.gz > /tmp/tt.fna
mkdir parts && cd parts && awk '/^>/{n=substr($1,2); f=n".fa"} {print > f}' /tmp/tt.fna && cd ..
ls parts/*.fa | xargs -P 22 -I{} sh -c \
  '/tmp/aug/bin/augustus --species=tetrahymena --gff3=on --genemodel=partial \
     --UTR=off --softmasking=0 {} > {}.gff3'
(echo '##gff-version 3'; cat parts/*.fa.gff3 | grep -v '^#') > /tmp/tt_augustus.gff3

# --genetic-code 6: TAA and TAG are glutamine in the Tetrahymena nucleus, so
# the genome probe must not read them as stops (docs/benchmark.md 2.1, 4.5).
# No --stop-outside-cds: the `tetrahymena` parameter set includes the stop.
python3 benchmark/score.py \
    --reference /tmp/panel/Tetrahymena_thermophila/*_genomic.gff.gz \
    --prediction /tmp/tt_augustus.gff3 --species Tetrahymena_thermophila \
    --declaration benchmark/validation/augustus-partial-Tetrahymena_thermophila.yaml \
    --genome /tmp/panel/Tetrahymena_thermophila/*_genomic.fna.gz \
    --genetic-code 6 --out /tmp/tt.json
```

*A. mellifera* is the same with `--species=honeybee1 --softmasking=1
--uniqueGeneId=true`, the default genetic code, and `--stop-outside-cds`,
which AUGUSTUS itself asks for on stderr: `honeybee1` is one of the 44 species
parameter sets of 166 that set `stopCodonExcludedFromCDS true`.

Reproducing the Tiberius run. Tiberius is not a dependency of this benchmark
either; the run used the published container, so nothing is installed on the
host:

```
python3 benchmark/fetch.py --species Takifugu_rubripes --what fasta --dest /tmp/panel
python3 benchmark/fetch.py --species Takifugu_rubripes --what gff   --dest /tmp/panel
zcat /tmp/panel/Takifugu_rubripes/*_genomic.fna.gz > /tmp/tibrun/fugu.fa

# --out takes both dialects; writing both is what makes the GTF/GFF3
# equality in docs/benchmark.md 6 checkable on real output.
# batch_size is auto-computed from VRAM (4 on a 16 GB card); the model_cfg
# fixes seq_len at 400,050.
docker run --gpus all --rm -v /tmp/tibrun:/data \
  larsgabriel23/tiberius@sha256:2c3bddda32cc621b805de40dc0395942cb5dd8a5766b1fe6c98fba845740f9bd \
  tiberius --genome /data/fugu.fa --model_cfg vertebrates \
           --out /data/tiberius_fugu.gtf /data/tiberius_fugu.gff3

# score each dialect separately; the two results must agree everywhere except
# the fields that report what the file states about its own convention.
for d in gtf gff3; do
  out=benchmark/validation/tiberius-Takifugu_rubripes.json
  [ $d = gtf ] && out=benchmark/validation/tiberius-Takifugu_rubripes-gtf.json
  python3 benchmark/score.py \
      --reference /tmp/panel/Takifugu_rubripes/*_genomic.gff.gz \
      --prediction /tmp/tibrun/tiberius_fugu.$d --species Takifugu_rubripes \
      --declaration benchmark/validation/tiberius-Takifugu_rubripes.yaml \
      --genome /tmp/panel/Takifugu_rubripes/*_genomic.fna.gz --out $out
done
```

No `--stop-outside-cds`: Tiberius includes the stop codon in the CDS, and the
scorer confirms it from both routes (the GTF's `stop_codon` features and the
genome probe, which finds 23,946 of 23,948 chains ending on a stop).

This run is also what `splice.short_gaps` (§4.3) was added for. Tiberius emits
413 CDS gaps below the 20 bp intron floor on fugu; the reference has 1,124.
All 413 predicted gaps satisfy the donor and acceptor motif masks — 312 on
their own bases, 101 by borrowing an exon base at each end, and those 101 are
exactly the two one-base contexts the masks allow, `A|G|T` (72) and `A|G|C`
(29). Of the 1,124 reference gaps, 7 pass the same combined test — the same
two borrowed contexts, `A|G|C` ×4 and `A|G|T` ×3 — and 1,117 fail it; every
one is 1 or 2 bp and none is a multiple of three. Helixer on the same genome
emits none.

Two things that earlier wording here got wrong, both raised by engels in
`relay/messages/20260909T222452Z-engels-0022.md`:

- **These are two aggregates, not an intersection.** Both sides are counted
  from their own gaps and the report keeps no coordinates, so it cannot say
  whether any predicted gap sits at a reference gap. The distributions differ
  strongly; they are not disjoint in motif class, since both sides contain
  `A|G|T` and `A|G|C` gaps.
- **Failing the combined test is not "satisfying neither mask."** The test is
  a conjunction, so its failure covers three cases, and `motif_by_class` now
  reports them apart. Of the 1,117 fugu reference gaps that fail, **104 carry
  a GT/GC donor with no AG after it**, **61 carry an AG with no donor before
  it**, and 952 fail both windows. `motif_unresolved_windows` is 0 on every
  side of every committed run, so no gap on this panel sits over an `N`: the
  class counts are measurements and not silent ambiguity.

The floor is therefore not hiding a disagreement about micro-introns: it is
separating an annotation's frameshift encoding from a decoder that constrains
splice-site composition without constraining intron duration. It also puts a
number on a half-masked decoder: a donor-only mask would admit 104 of the
1,124 fugu reference frameshift steps (9.3%) as intron starts, against 7 for
the conjunction.

The `heldout_seen_in_pretraining: no` in that declaration is a set operation,
not a judgement, and it can be rechecked without the container:

```
git clone --filter=blob:none --no-checkout https://github.com/Gaius-Augustus/Tiberius /tmp/tib
git -C /tmp/tib checkout c6d92f2fa15cee0bc845141395216cea89ec89ec -- model_cfg
```

`training_species` is a bracketed whitespace-separated bare sequence in
`fungi.yaml`, `insecta.yaml` and the plant/algal configs and a `-` list with
trailing accession comments in `vertebrates.yaml`, so a parser written for one
shape returns nothing useful on the other. Across all nine configs, five panel
species appear in some training list (*M. musculus*, *A. mellifera*,
*S. cerevisiae*, *S. pombe*, *N. crassa*) and fifteen appear in none.
