Body of the `alert` to `human` for the two cost-baseline rows (the
resubmission of marx-0033 in the `relay/templates/compute-request.md` form,
as the decision human-0013 asks). Posted once `gagarin` has said hello.

## Task

T-human-009 (done; the two rows feed `docs/cost-baseline/measured.tsv` and
section 3 of `docs/cost-baseline.md` through a follow-up pull request from
`work/T-human-009-marx`). The document currently carries them as
"not attempted" (EGAPx) and "not measurable on this runner" (Tiberius).

## What to run

Two runs, three Slurm jobs, one attempt each, from branch
`work/T-human-009-marx` at commit `3ce9cb4` (scripts under
`docs/cost-baseline/cluster/`, described in its `README.md`).

### Run 1: EGAPx v1.0 on *Ciona intestinalis*

- Command or job script: `docs/cost-baseline/cluster/egapx_ciona.sbatch <checkout> <scratch>/egapx`,
  through `cluster-run.sh submit ... --cpus 20 --mem 240 --gpus 0 --time 12:00:00 --name egapx-ciona`.
  The script clones `ncbi/egapx` at tag `v1.0`, runs `ui/egapx.py input.yaml -e singularity`
  with an input of genome and taxid 7719 only (no RNA-seq, proteins from
  EGAPx's own taxon-matched NCBI sets), caps Nextflow CPU and memory labels
  to the allotment, and times the driver under GNU time with a Nextflow trace.
- Container or environment: `docker://ncbi/egapx:1.0` pulled by EGAPx itself
  through singularity or apptainer; Nextflow 23.10 or later, Python 3.11 or
  later and PyYAML on the node.
- Inputs (verified against NCBI MD5 by `benchmark/fetch.py`; SHA-256 recorded after download):
  - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/327/825/GCF_018327825.1_ASM1832782v2/GCF_018327825.1_ASM1832782v2_genomic.fna.gz, 42,164,729 bytes, md5 c9eacbcc32a5182a496cfa7e057b9aa2
  - https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/327/825/GCF_018327825.1_ASM1832782v2/GCF_018327825.1_ASM1832782v2_genomic.gff.gz, 7,079,141 bytes, md5 5acd9a693eea48d5a944fa7f7273cf21
  - EGAPx support files and protein sets, fetched by EGAPx on first run (several GB; not checksummed here).
- Outputs to return: `<scratch>/egapx/return/` (all text, under 1 MB):
  `summary.tsv`, `trace.txt`, `egapx.time`, `resource_caps.txt`,
  `egapx-Ciona_intestinalis.json`, `egapx-Ciona_intestinalis.yaml`,
  `inputs.sha256`, `egapx_commit.txt`, `versions.txt`, `hardware.txt`,
  `disk_usage.txt`, the stdout and stderr tails. The prediction GFF (about
  10 MB compressed) stays in scratch.
- Scoring or summary command: run by the script:
  `summarize_trace.py` for the cost row and
  `benchmark/score.py --reference <gff.gz> --prediction out/complete.genomic.gff --species Ciona_intestinalis --declaration ... --genome <fna.gz>`
  for accuracy.

### Run 2: Tiberius 2.0.7 on eleven panel species, two jobs

- Command or job script: `docs/cost-baseline/cluster/tiberius_panel.sbatch <checkout> <scratch>/tiberius <species...>`,
  through `cluster-run.sh submit ... --cpus 8 --mem 128 --gpus 1 --time 24:00:00`, twice:
  `--name tiberius-ten` with
  `Saccharomyces_cerevisiae Schizosaccharomyces_pombe Neurospora_crassa Drosophila_melanogaster Apis_mellifera Takifugu_rubripes Gallus_gallus Danio_rerio Xenopus_tropicalis Mus_musculus`
  and `--name tiberius-human` with `Homo_sapiens`. Same scratch directory is
  fine (per-species subdirectories; the image and weights are shared).
  Each species: `tiberius.py --genome ... --out tiberius.gff3 --model_cfg <vertebrates|insecta|fungi> --batch_size 8 --seq_len 400050 --no_softmasking`
  under GNU time, GPU memory sampled every 5 s, then scored.
- Container or environment:
  `docker://larsgabriel23/tiberius@sha256:2c3bddda32cc621b805de40dc0395942cb5dd8a5766b1fe6c98fba845740f9bd`
  (the image behind `benchmark/validation/tiberius-Takifugu_rubripes.yaml`),
  run with `apptainer exec --nv --writable-tmpfs`. Weights: vertebrates in the
  image; insecta and fungi fetched from bioinf.uni-greifswald.de on first use
  (tarballs also fetched to scratch and checksummed).
- Inputs: the 22 files in `docs/cost-baseline/cluster/inputs.tsv` for the
  eleven species (genome FASTA and RefSeq GFF each, 3.25 GB compressed in
  total; URL, size and NCBI MD5 per file there), verified by `benchmark/fetch.py`.
- Outputs to return: `<scratch>/tiberius/return/`: `summary.tsv` (one row per
  species), per-species `tiberius-<sp>.time`, `.stderr.tail`, `.json` and
  `.yaml`, `inputs.sha256`, `image.sha256`, `weights.sha256`, `gpu.txt`,
  `hardware.txt`. Predictions stay in scratch.
- Scoring or summary command: run by the script per species:
  `benchmark/score.py --reference <gff.gz> --prediction tiberius.gff3 --species <sp> --declaration ... --genome <fna.gz>`.

## Resources

| job | CPUs | memory (GB) | GPUs | wall clock | scratch (GB) |
|---|---|---|---|---|---|
| egapx-ciona | 20 | 240 | 0 | 12:00 | 150 |
| tiberius-ten | 8 | 128 | 1 | 24:00 | 60 (shared) |
| tiberius-human | 8 | 128 | 1 | 24:00 | (same) |

- Estimated total: EGAPx 40 to 100 CPU-hours, 0 GPU-hours (README: the
  144 Mb fly with one RNA-seq run took 71 CPU-h; this run has no reads and a
  140 Mb genome, but 20 CPUs instead of 32 and no fan-out across instances,
  so 3 to 8 hours wall). Tiberius 6 to 16 GPU-hours across both jobs (10.60
  Gb at the paper's 1.9 GPU-s/Mb on an A100 is 5.6 h; the A5000's dense
  BF16 tensor throughput is about a third of the A100's, so up to 2.8x
  that), of which the human job is 1.6 to 4.6 GPU-hours; CPU-hours at most
  the 8-core allotment times wall, in practice about twice the GPU-hours.
  Both jobs fit 24 hours with margin; if the ten-species job is cut off, the
  species already finished are still returned.

## Why the laptop budget is not enough

EGAPx states a 32-CPU, 256 GB minimum and this runner has 4 cores and 15 GB;
Tiberius needs a GPU with at least 8 GB and this runner has none, and the
panel is 10.6 Gb.

## Where the results go

`docs/cost-baseline/measured.tsv` (one EGAPx row, eleven Tiberius rows) and
`docs/cost-baseline.md` section 3 (new subsections replacing the two
"documented only" cells in the section 1 table), with actual CPU-hours and
GPU-hours recorded in the T-human-009 log, through a pull request from
`work/T-human-009-marx`.
