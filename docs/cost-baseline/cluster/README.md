# Cluster runs for the cost baseline (T-human-009)

`docs/cost-baseline.md` section 5.4 names two rows the document could not
measure on a laptop: EGAPx end to end and Tiberius on the benchmark panel.
This directory holds the job scripts, input manifest and submission
declarations for those two runs, submitted by `gagarin` through
`relay/bin/cluster-run.sh` after a coordinator `decision` on the compute
request that cites this directory and its commit. The scripts write only
under their scratch directory; everything that comes back to the relay is
text under `<scratch>/return/`.

| file | what it is |
|---|---|
| `inputs.tsv` | the twelve assemblies (genome FASTA and RefSeq GFF each): NCBI URL, size in bytes, and the MD5 NCBI publishes in `md5checksums.txt`, read 2026-09-15. `benchmark/fetch.py` verifies these on download; the jobs record SHA-256 after download in `return/inputs.sha256` |
| `egapx_ciona.sbatch` | EGAPx v1.0 on *Ciona intestinalis* (140 Mb, `heldout`), genome and taxid only, one node, Nextflow local executor under `-e singularity`; per-process CPU and memory labels capped to the node; driver timed with GNU time, workers through the Nextflow trace; scored with `benchmark/score.py` |
| `tiberius_panel.sbatch` | Tiberius 2.0.7 (the image digest pinned in `benchmark/validation/`) on the panel species the `vertebrates`, `insecta` and `fungi` models cover, one GPU, batch 8, chunk 400,050, softmasking off; each species timed with GNU time and `nvidia-smi` sampled every 5 s; scored per species |
| `summarize_trace.py` | folds a Nextflow trace and a GNU time file into one `measured.tsv`-style row (CPU-s from `realtime x %cpu` over tasks and from user+sys on the driver; peak RSS as the largest task) |
| `declarations/` | the section 3.3 declaration for each run, `heldout_seen_in_pretraining` taken from the model configs' literal `training_species` lists at Tiberius commit `c6d92f2` (yes for *M. musculus*, *A. mellifera* and the three fungi; no for the other six), and for EGAPx marked `yes` because its taxon-matched protein sets are not filtered for the target |

## Species and models

| species | Gb | Tiberius model | in its training list |
|---|---|---|---|
| Homo_sapiens | 3.10 | vertebrates | no |
| Mus_musculus | 2.73 | vertebrates | yes |
| Xenopus_tropicalis | 1.45 | vertebrates | no |
| Danio_rerio | 1.45 | vertebrates | no |
| Gallus_gallus | 1.05 | vertebrates | no |
| Takifugu_rubripes | 0.38 | vertebrates | no |
| Apis_mellifera | 0.23 | insecta | yes |
| Drosophila_melanogaster | 0.14 | insecta | no |
| Neurospora_crassa | 0.04 | fungi | yes |
| Schizosaccharomyces_pombe | 0.01 | fungi | yes |
| Saccharomyces_cerevisiae | 0.01 | fungi | yes |

Total 10.60 Gb (`genome_bp` in `benchmark/panel.tsv`). The eleven are split
into two jobs so the human genome does not share one 24-hour window with the
other ten. *Ciona intestinalis* (0.14 Gb) is the EGAPx genome.

## What comes back

- EGAPx: `summary.tsv` (wall, CPU-s by trace and by GNU time, peak task RSS,
  CPU-s per Mb, task counts), `trace.txt`, `egapx.time`, `resource_caps.txt`
  (every label the script lowered), `egapx-Ciona_intestinalis.json` from the
  scorer and its declaration, `inputs.sha256`, commits and versions. The
  prediction GFF stays in scratch (about 10 MB compressed).
- Tiberius: `summary.tsv` with one row per species (wall, user and system
  CPU-s, peak RSS, peak GPU memory, gene count, batch size), a GNU time file
  and a stderr tail per species, `tiberius-<species>.json` and declaration per
  scored species, `inputs.sha256`, `image.sha256`, `weights.sha256`, `gpu.txt`.

The rows go into `docs/cost-baseline/measured.tsv` and section 3 of
`docs/cost-baseline.md`, with the actual CPU-hours and GPU-hours in the task
log, as the charter's Compute paragraph asks.

## Known risks, one attempt each

- Neither script has run on the target cluster; this runner has no Slurm,
  no GPU and no container runtime. Each fails loudly with its logs in
  `return/` rather than partially succeeding in silence.
- First attempt, 2026-09-15 (Slurm job 95577, `relay/artifacts/T-gagarin-001/
  egapx-ciona/`): the EGAPx job died 24 s in, before EGAPx started, because
  node03's Python 3.14 has no `ensurepip` and `python3 -m venv` refused to
  build the venv. `egapx_ciona.sbatch` step 2 now falls back to a venv built
  `--without-pip` with pip bootstrapped from get-pip.py, then to the PATH
  interpreter if it already imports `yaml`, then to `pip install --user`; a
  half-made venv from an earlier attempt is rebuilt, and the interpreter used
  is written to `return/versions.txt`. Each path was exercised on the writing
  runner with a shim standing in for node03's interpreter; the real node has
  still not run EGAPx.
- EGAPx's minimum machine is 32 CPUs and 256 GB; the cluster node has 20
  CPUs. The script caps Nextflow labels at the allotment, so wall clock will
  be longer than the README's AWS figures and any process that truly needs
  more than the node has will fail and say so in the trace.
- EGAPx downloads several GB of support files and taxon protein sets on the
  first run; the compute node needs outbound network access, as does the
  Tiberius job for the NCBI genomes and the weights.
- Tiberius on the A5000 (24 GB, GA102) runs at batch 8 by the README's
  RTX 3090 guidance; a memory error on the largest chunk would show in the
  per-species stderr tail and the summary row's exit status, and the loop
  continues with the next species.
