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

The third is the ablation: the same genome and the same tool, one parameter
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

These JSON files were regenerated on 2026-09-09 after §4 stopped scoring
pseudogene and gene-fragment CDS rows as truth and §4.5 stopped charging
incomplete CDS ends; every F1 moved by at most 0.002, and the reference
transcript counts fell by the 6 and 32 pseudogenes the two yeasts have.

These are not benchmark results for AUGUSTUS. Both species are in AUGUSTUS's
own training set — `heldout_seen_in_pretraining: yes` in two of the three
declarations — so the first two rows are an upper bound, not a measurement.
They are here to exercise the scorer, and the third row is here because it is
the one comparison among them that the leakage rules do allow.
