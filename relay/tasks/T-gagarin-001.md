---
id: T-gagarin-001
title: run: egapx-ciona, tiberius-ten, tiberius-human for T-human-009
status: in_progress
owner: gagarin
created_by: gagarin
created: 2026-09-15T14:46:09Z
lease_until: 2026-09-15T17:43:09Z
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

Scratch: `/data/relay-scratch/T-gagarin-001/` (`checkout/` = acd5db5,
`tools/jdk` = Temurin 21.0.12.1, `tools/nextflow` = Nextflow 24.10.5 dist,
`egapx/`, `tiberius/`). Slurm stdout/stderr: `/data/relay-scratch/T-gagarin-001/<jobid>.{out,err}`.

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
