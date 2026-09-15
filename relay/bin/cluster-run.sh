#!/usr/bin/env bash
# Submit or poll a Slurm job on behalf of a relay task. Usage:
#   cluster-run.sh submit --task T-gagarin-NNN --decision MSG-ID \
#       --cpus N --mem GB --gpus N --time HH:MM:SS --name LABEL SCRIPT [ARGS...]
#   cluster-run.sh status --task T-gagarin-NNN
#
# Every submission must name the relay task that owns it and the coordinator
# decision that granted it. All four resource limits are required and are
# passed to sbatch as hard limits; Slurm kills the job when it exceeds them.
# The wrapper refuses anything over the cluster ceilings below (override
# per machine with RELAY_CLUSTER_MAX_*). Job ids go to
# relay/artifacts/<task>/jobs.tsv so later ticks can poll them, and Slurm
# stdout/stderr go to $RELAY_CLUSTER_SCRATCH/<task>/ (never into the repo).
set -euo pipefail

max_cpus="${RELAY_CLUSTER_MAX_CPUS:-20}"
max_mem="${RELAY_CLUSTER_MAX_MEM_GB:-256}"
max_gpus="${RELAY_CLUSTER_MAX_GPUS:-2}"
max_hours="${RELAY_CLUSTER_MAX_HOURS:-24}"
scratch="${RELAY_CLUSTER_SCRATCH:-$HOME/relay-scratch}"

here="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$here"

die() { echo "cluster-run.sh: $*" >&2; exit 2; }

mode="${1:-}"; shift || true
task="" decision="" cpus="" mem="" gpus="" walltime="" label=""
while [ $# -gt 0 ]; do
  case "$1" in
    --task) task="$2"; shift 2 ;;
    --decision) decision="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --mem) mem="$2"; shift 2 ;;
    --gpus) gpus="$2"; shift 2 ;;
    --time) walltime="$2"; shift 2 ;;
    --name) label="$2"; shift 2 ;;
    --) shift; break ;;
    -*) die "unknown option $1" ;;
    *) break ;;
  esac
done

[[ "$task" =~ ^T-[a-z][a-z0-9-]*-[0-9]{3,}$ ]] || die "--task must be a relay task id"
[ -f "relay/tasks/$task.md" ] || die "no such task file relay/tasks/$task.md"
artdir="relay/artifacts/$task"
jobs="$artdir/jobs.tsv"

case "$mode" in
  submit)
    [ -n "$decision" ] || die "--decision MSG-ID is required (the coordinator's grant)"
    ls relay/messages/ | grep -q "^${decision}\.md$" || die "decision $decision is not a message in relay/messages/"
    grep -q '^type: decision$' "relay/messages/$decision.md" || die "$decision is not a decision"
    [ -n "$cpus" ] && [ -n "$mem" ] && [ -n "$gpus" ] && [ -n "$walltime" ] \
      || die "--cpus, --mem, --gpus and --time are all required"
    [ "$cpus" -ge 1 ] && [ "$cpus" -le "$max_cpus" ] || die "cpus $cpus outside 1..$max_cpus"
    [ "$mem" -ge 1 ] && [ "$mem" -le "$max_mem" ] || die "mem $mem GB outside 1..$max_mem"
    [ "$gpus" -ge 0 ] && [ "$gpus" -le "$max_gpus" ] || die "gpus $gpus outside 0..$max_gpus"
    [[ "$walltime" =~ ^([0-9]+):([0-9]{2}):([0-9]{2})$ ]] || die "--time must be HH:MM:SS"
    hours="${BASH_REMATCH[1]}"
    [ "$hours" -lt "$max_hours" ] || [ "$hours" -eq "$max_hours" -a "${BASH_REMATCH[2]}${BASH_REMATCH[3]}" = "0000" ] \
      || die "time $walltime exceeds ceiling ${max_hours}h"
    script="${1:-}"; [ -n "$script" ] || die "SCRIPT is required"
    [ -f "$script" ] || die "script $script not found"
    label="${label:-$task}"

    mkdir -p "$artdir" "$scratch/$task"
    [ -f "$jobs" ] || printf 'job_id\tsubmitted\tdecision\tlabel\tcpus\tmem_gb\tgpus\ttime\tscript\n' >"$jobs"

    gres=()
    [ "$gpus" -gt 0 ] && gres=(--gres="gpu:$gpus")
    job_id="$(sbatch --parsable \
      --job-name="$label" \
      --cpus-per-task="$cpus" --mem="${mem}G" "${gres[@]}" --time="$walltime" \
      --output="$scratch/$task/%j.out" --error="$scratch/$task/%j.err" \
      --comment="relay $task $decision" \
      "$script" "${@:2}")"
    job_id="${job_id%%;*}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$job_id" "$(date -u +%FT%TZ)" "$decision" "$label" "$cpus" "$mem" "$gpus" "$walltime" "$script" >>"$jobs"
    echo "submitted $job_id ($label) for $task under $decision; logs in $scratch/$task/"
    ;;
  status)
    [ -f "$jobs" ] || die "no jobs recorded for $task"
    # One line per job: id, state, elapsed, max RSS, CPU time. Records the
    # numbers the charter asks agents to log after a run.
    tail -n +2 "$jobs" | cut -f1 | while read -r id; do
      sacct -j "$id" --noheader --parsable2 \
        --format=JobID,State,Elapsed,MaxRSS,TotalCPU,AllocTRES 2>/dev/null \
        | grep -v '\.batch\|\.extern' || echo "$id|unknown|||"
    done
    ;;
  *)
    die "usage: cluster-run.sh submit|status ... (see header)"
    ;;
esac
