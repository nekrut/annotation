---
id: T-gagarin-001
title: run: egapx-ciona, tiberius-ten, tiberius-human for T-human-009
status: blocked
owner: null
created_by: gagarin
created: 2026-09-15T14:46:09Z
lease_until: null
depends_on: [T-human-009]
touches: []
pr: null
---

## Goal

Execute the cluster grant in decision `20260915T143249Z-human-0016` (human,
2026-09-15) for marx's compute request `20260915T015342Z-marx-0042`
(T-human-009): EGAPx v1.0 on *Ciona intestinalis* and Tiberius 2.0.7 on eleven
panel species, scripts from `work/T-human-009-marx` at `acd5db5`
(`docs/cost-baseline/cluster/`). Full request text: `relay/messages/20260915T015342Z-marx-0042.md`.

Granted limits (decision: "all three jobs at the limits requested in marx-0042,
in the order egapx-ciona, tiberius-ten, tiberius-human; the two Tiberius jobs
may be queued at the same time"):

| job | CPUs | memory (GB) | GPUs | wall clock | scratch (GB) |
|---|---|---|---|---|---|
| egapx-ciona | 20 | 240 | 0 | 12:00 | 150 |
| tiberius-ten | 8 | 128 | 1 | 24:00 | 60 (shared) |
| tiberius-human | 8 | 128 | 1 | 24:00 | (same dir) |

Decision conditions for gagarin: scratch under `RELAY_CLUSTER_SCRATCH` on
`/data`; apptainer for containers; stage a JDK 17+ and the Nextflow launcher
under scratch, user-level, with `JAVA_HOME` and `PATH` set inside the job;
if staging fails, block and ask; Tiberius jobs proceed regardless; return
summary files and usage with a `note` on T-human-009.

Scratch: `/data/relay-scratch/T-gagarin-001/` (`checkout/` = 0be4b1d since 16:44Z, acd5db5 before,
`tools/jdk` = Temurin 21.0.12.1, `tools/nextflow` = Nextflow 24.10.5 dist,
`egapx/`, `tiberius/`, `tiberius-human/`). Slurm stdout/stderr: `/data/relay-scratch/T-gagarin-001/<jobid>.{out,err}`.

## Definition of done

Each of the three jobs has finished (or failed) once; its return files under
5 MB, log tails, and usage are in `relay/artifacts/T-gagarin-001/<label>/`;
one `note` per finished job posted to marx on T-human-009; any resubmission
question posted to human. Task released `done`, or `blocked` on a failure
awaiting a decision.

## Log
- 2026-09-15 gagarin: claimed.
- 2026-09-15T14:46Z gagarin: staged Temurin JDK 21.0.12.1 and Nextflow 24.10.5 (self-contained dist) under `tools/`; `nextflow -version` verified on head with `NXF_HOME=egapx/nxf`. Cloned acd5db5 to `checkout/`.
- 2026-09-15T14:47Z gagarin: `sbatch --test-only` showed `--mem=240G` (245760 MiB) is unschedulable: node03 has RealMemory 245000 MiB and no node is larger. Submitted egapx-ciona at 239 GB (largest schedulable, within the 240 GB grant) rather than not at all. 128 GB + 1 GPU fits only node03 (node02 has 119000 MiB), so the two Tiberius jobs serialize on node03 as written.
- 2026-09-15T14:47Z gagarin: submitted job 95577 egapx-ciona (20 CPU, 239 GB, 0 GPU, 12:00:00), job 95578 tiberius-ten (8 CPU, 128 GB, 1 GPU, 24:00:00), job 95579 tiberius-human (same). `JAVA_HOME` and `PATH` (jdk, nextflow) exported into the jobs via sbatch's environment propagation; scripts unmodified.
- 2026-09-15T14:47Z gagarin: job 95577 egapx-ciona FAILED on node03 after 24 s. Inputs fetched and MD5-verified, nextflow/apptainer/GNU time checks passed (versions.txt: Nextflow 24.10.5, apptainer 1.4.5, Python 3.14.4), then `python3 -m venv venv` failed: ensurepip not available (python3.14-venv not installed on node03). EGAPx never started. Return files, log tails and usage in `relay/artifacts/T-gagarin-001/egapx-ciona/`. Usage: elapsed 00:00:24, 20 CPUs x 239 GB allocated, 0 GPU; no RSS/CPU-time (died before the timed section; `sacct` is disabled on this cluster, so `cluster-run.sh status` reports `unknown` and usage comes from scontrol while a job lives and from the scripts' GNU time files afterwards). Broken `egapx/venv` must be removed before any resubmission (the script skips venv creation if the dir exists). Question to human posted on whether to stage a user-level Python with venv and resubmit.
- 2026-09-15T14:48Z gagarin: job 95578 tiberius-ten RUNNING on node03 since 14:47Z; job 95579 tiberius-human PENDING (Resources: needs node03). Lease renewed.
- 2026-09-15T15:43Z gagarin: job 95578 tiberius-ten RUNNING on node03 (55 min elapsed of 24 h; 8 CPU, 128 GB, 1 GPU). 5 of 10 species done, all exit 0 (S. cerevisiae 74 s, S. pombe 78 s, N. crassa 156 s, D. melanogaster 707 s, A. mellifera 655 s; peak RSS 5.5 to 26.8 GB; GPU 22.3 GB); T. rubripes in progress. Job 95579 tiberius-human still PENDING (Resources, waits for node03). `cluster-run.sh status` reports `unknown` for all three (sacct disabled); state taken from squeue/scontrol. Lease renewed.
- 2026-09-15T15:43Z gagarin: read marx-0050 (note): egapx_ciona.sbatch rewritten at 0be4b1d on work/T-human-009-marx to not need ensurepip, and 116 GB on node02 is acceptable to the requester for tiberius-human. No new decision from human on gagarin-0003, so nothing resubmitted; egapx-ciona stays failed pending a decision. Push via SSH failed this tick (no ssh-agent in the runner shell); pushed over HTTPS with the gh credential helper instead.
- 2026-09-15T16:44Z gagarin: human-0017 (answer to gagarin-0003, within the grant of decision human-0016) authorised both deviations: resubmit egapx-ciona from work/T-human-009-marx at 0be4b1d at the same limits, and cancel 95579 to resubmit tiberius-human at 116 GB on node02. Deviation 1: `checkout/` advanced acd5db5 -> 0be4b1d (whole-tree diff is only `docs/cost-baseline/cluster/{README.md,REQUEST.md,egapx_ciona.sbatch}`, so the running tiberius-ten job, which reads `benchmark/` from the same checkout, is unaffected). Removed the half-made `egapx/venv`. Submitted job 95580 egapx-ciona (20 CPU, 239 GB, 0 GPU, 12:00:00, `<checkout> egapx/`), JAVA_HOME and PATH (jdk, nextflow) and NXF_HOME propagated as before. It is PENDING (Resources) behind tiberius-ten: 239 GB only fits node03, which holds 128 GB for 95578, so EGAPx starts when the ten-species job ends.
- 2026-09-15T16:44Z gagarin: Deviation 2: cancelled 95579 (was PENDING, never started) and submitted job 95581 tiberius-human (8 CPU, 116 GB, 1 GPU, 24:00:00). The request's "(same dir)" was not followed: `tiberius_panel.sbatch` truncates `return/summary.tsv` at start, so a concurrent run in `tiberius/` would have wiped the ten-species rows. It runs in `tiberius-human/` with the 14.7 GB image and the three weight tarballs hardlinked from `tiberius/` (script skips the pull when the files exist; checksums are recomputed into its own `return/`). Backfilled onto node02 at 16:45Z (AllocTRES cpu=8,mem=116G, node02's GPU); its `return/` is being written under `/data`, so the NFS mount on node02 is live (ssh to node02 is not possible from this runner shell; the nfsd client list on head also shows node02 connected).
- 2026-09-15T16:46Z gagarin: job 95578 tiberius-ten RUNNING on node03 (1 h 58 min elapsed). 6 of 10 species done, all exit 0 (T. rubripes 1263 s, peak RSS 27.1 GB, GPU 22.3 GB); G. gallus in progress since 15:59Z; D. rerio, X. tropicalis, M. musculus remain. Partial return files (all species finished so far, 232 KB) copied to `relay/artifacts/T-gagarin-001/tiberius-ten/` so marx can start folding; the directory is refreshed in place on later ticks. Lease renewed to 18:46Z.
- 2026-09-15T17:43Z gagarin: job 95578 tiberius-ten RUNNING on node03 (2 h 56 min elapsed of 24 h; AllocTRES cpu=8,mem=128G,node03 GPU). 7 of 10 species done, all exit 0; G. gallus finished 16:52Z (3156 s wall, 3300 s user CPU, peak RSS 46.1 GB, GPU 22.3 GB, 19058 genes). D. rerio in progress since 16:53Z; X. tropicalis and M. musculus remain. `relay/artifacts/T-gagarin-001/tiberius-ten/` refreshed in place with the 7-species return files (272 KB) and a scontrol snapshot. `cluster-run.sh status` still reports `unknown` (sacct disabled).
- 2026-09-15T17:43Z gagarin: job 95581 tiberius-human RUNNING on node02 since 16:45Z (58 min elapsed; AllocTRES cpu=8,mem=116G, node02 GPU). Image and weights verified by checksum, GRCh38.p14 fasta and gff fetched and MD5-verified by 16:49Z; the Tiberius run on H. sapiens is in progress, no per-species output yet. Slurm stderr shows the job could not chdir to the submit directory on node02 and fell back to /tmp (same message as 95578 on node03; home is not mounted on compute nodes); the script uses absolute paths under /data so this had no effect. Provenance files copied to `relay/artifacts/T-gagarin-001/tiberius-human/` (refreshed in place).
- 2026-09-15T17:43Z gagarin: job 95580 egapx-ciona PENDING (Resources) on node03 behind 95578; Slurm's estimated start is 2026-09-16T10:47Z (worst case, the tiberius-ten wall-clock limit). No new decision or question for gagarin in the inbox (marx-0052, lenin-0046, engels-0037, stalin-0044 concern PR #31 only). Lease renewed to 2026-09-15T19:43:00Z.
- 2026-09-15T18:42Z gagarin: job 95578 tiberius-ten RUNNING on node03 (3 h 55 min elapsed of 24 h; AllocTRES cpu=8,mem=128G,node03 GPU). 8 of 10 species done, all exit 0; D. rerio finished 18:05Z (4232 s wall, 4440 s user CPU, peak RSS 38.3 GB, GPU 22.3 GB, 29799 genes). X. tropicalis in progress since 18:07Z; M. musculus remains. `relay/artifacts/T-gagarin-001/tiberius-ten/` refreshed in place with the 8-species return files (308 KB) and a scontrol snapshot.
- 2026-09-15T18:42Z gagarin: job 95581 tiberius-human RUNNING on node02 (1 h 57 min elapsed; AllocTRES cpu=8,mem=116G, node02 GPU). H. sapiens run in progress: `tiberius.gff3` in the work dir is 28.7 MB and still growing, no summary row yet. Provenance files and scontrol snapshot refreshed in `relay/artifacts/T-gagarin-001/tiberius-human/`.
- 2026-09-15T18:42Z gagarin: job 95580 egapx-ciona PENDING (Resources) behind 95578 on node03; Slurm estimate unchanged (worst case 2026-09-16T10:47Z). Inbox: marx-0053, lenin-0047, engels-0038, stalin-0045 are PR #31 notes, no decision or question for gagarin. No note posted this tick (nothing finished; artifacts refreshed in place). git over SSH still fails in this runner shell; pulled and pushing over HTTPS with the gh credential helper. Lease renewed to 2026-09-15T20:43:00Z.
- 2026-09-15T19:43Z gagarin: job 95581 tiberius-human COMPLETED, exit 0, on node02 (8 CPU, 116 GB, 1 GPU). Job ran 16:45:06Z to 19:25:17Z; Tiberius step wall 9186.3 s (2:33:06), user CPU 9590.82 s, sys 1945.13 s, peak RSS 59034196 kB (57.7 GB), peak GPU memory 22323 MiB, 28476 genes; score nucleotide F1 0.916, exon F1 0.825, transcript F1 0.649, locus F1 0.852. All return files (84 KB), Slurm stdout/stderr tails, `usage.txt` and the scontrol snapshot copied to `relay/artifacts/T-gagarin-001/tiberius-human/`; `tiberius.gff3` (45.7 MB) left in scratch. Slurm's only stderr is the chdir-to-/tmp notice (no effect). Note gagarin-0006 posted to marx on T-human-009. No sacct on this cluster (`cluster-run.sh status` reports `unknown`).
- 2026-09-15T19:43Z gagarin: job 95578 tiberius-ten RUNNING on node03 (4 h 55 min elapsed of 24 h; AllocTRES cpu=8,mem=128G,node03 GPU). 9 of 10 species done, all exit 0; X. tropicalis finished 19:19Z (4219 s wall, 4466 s user CPU, peak RSS 50.0 GB, GPU 22.3 GB, 25463 genes). M. musculus (last) in progress since 19:20Z, inputs MD5-verified. `relay/artifacts/T-gagarin-001/tiberius-ten/` refreshed in place with the 9-species return files and a scontrol snapshot.
- 2026-09-15T19:43Z gagarin: job 95580 egapx-ciona PENDING (Resources) behind 95578 on node03; Slurm estimate unchanged (worst case 2026-09-16T10:47Z). Inbox: lenin-0048, marx-0054, engels-0039, stalin-0046 are PR #31 notes and reviews, no decision or question for gagarin. git over SSH still fails in this runner shell; pulled and pushing over HTTPS with the gh credential helper. Task stays in_progress (two jobs outstanding). Lease renewed to 2026-09-15T21:43:00Z.
- 2026-09-15T20:42Z gagarin: job 95578 tiberius-ten RUNNING on node03 (5 h 55 min elapsed of 24 h; AllocTRES cpu=8,mem=128G,node03 GPU). Still 9 of 10 species done, all exit 0. M. musculus (last, GRCm39, 2.7 GB fasta) in progress since 19:20Z: `tiberius.gff3` in its work dir is 24.0 MB and was last written 20:37Z, no summary row yet; by the X. tropicalis rate (1451 Mb in 4219 s) it should finish around 21:40Z. `relay/artifacts/T-gagarin-001/tiberius-ten/` refreshed in place (9-species return files unchanged, scontrol snapshot updated).
- 2026-09-15T20:42Z gagarin: job 95580 egapx-ciona PENDING (Resources), SchedNodeList=node03 behind 95578; Slurm estimate unchanged (worst case 2026-09-16T10:47Z, the tiberius-ten limit; realistically shortly after M. musculus finishes). Inbox: marx-0055, lenin-0049, engels-0040, stalin-0047 are PR #31 notes and reviews, no decision or question for gagarin. No note posted this tick (nothing finished). git over SSH still fails in this runner shell; pulled and pushing over HTTPS with the gh credential helper. Task stays in_progress. Lease renewed to 2026-09-15T22:43:00Z.
- 2026-09-15T21:45Z gagarin: job 95578 tiberius-ten COMPLETED, exit 0 for all 10 species, on node03 (8 CPU, 128 GB, 1 GPU). Job ran 14:47:34Z to 21:25:18Z (6 h 37 min 44 s); M. musculus finished 21:23Z (7268.7 s wall, 7653.82 s user CPU, peak RSS 45.2 GB, GPU 22.3 GB, 25321 genes). Totals: tiberius wall 21807.8 s, CPU user 23196.32 s, sys 4365.29 s, peak RSS 50029 MB (X. tropicalis), GPU 22323 MiB every species, 173735 genes over 7502.5 Mb. All return files, `scores.tsv` (per-species F1 extracted from the json), Slurm stdout/stderr tails and `usage.txt` in `relay/artifacts/T-gagarin-001/tiberius-ten/` (408 KB); the 20:42Z scontrol snapshot is the last one (job left Slurm before this tick; sacct disabled). GFF3 predictions (1.8 to 48 MB) left in scratch. Note gagarin-0007 posted to marx on T-human-009.
- 2026-09-15T21:45Z gagarin: job 95580 egapx-ciona FAILED on node03, 21:25:19Z to 21:25:38Z (19 s of 20 CPU, 239 GB, 0 GPU). The 0be4b1d venv bootstrap worked (venv --without-pip, get-pip.py, PyYAML 6.0.3); EGAPx v1.0 exited 1 in 0.16 s: "Local cache directory .../egapx/cache does not exist" (ui/egapx.py line 3118 requires the -lc path to exist; the script passes -lc "$SCRATCH/cache" but never mkdirs it). Nextflow never started; summary.tsv holds a zero row. Return files, slurm-95580 tails and a two-attempt `usage.txt` in `relay/artifacts/T-gagarin-001/egapx-ciona/` (refreshed in place, 95577 files kept). Usage: driver user 0.12 s, sys 0.02 s, max RSS 35.9 MB; no scontrol snapshot (job started and ended between ticks). Question gagarin-0008 to human: resubmit from a fixed commit (recommended) or pre-create the directory and resubmit 0be4b1d. Nothing submitted.
- 2026-09-15T21:45Z gagarin: definition of done reached for the Tiberius half (both jobs done, notes posted); EGAPx failed twice and awaits a decision. Releasing the task `blocked`. Scratch under /data/relay-scratch/T-gagarin-001/ kept intact for a resubmission. git over SSH still fails in this runner shell; pulled and pushing over HTTPS with the gh credential helper.
- 2026-09-15 gagarin: -> blocked.
