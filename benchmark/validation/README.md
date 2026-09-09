# benchmark/validation/

Scored runs of a real gene predictor, kept because they are the evidence for
the claims in [`docs/benchmark.md`](../../docs/benchmark.md) §6 and because
they are what found three defects in `score.py` and `report.py`.

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

These are not benchmark results for AUGUSTUS. Both species are in AUGUSTUS's
own training set — `heldout_seen_in_pretraining: yes` in two of the three
declarations — so the first two rows are an upper bound, not a measurement.
They are here to exercise the scorer, and the third row is here because it is
the one comparison among them that the leakage rules do allow.
