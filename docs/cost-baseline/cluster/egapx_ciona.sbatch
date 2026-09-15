#!/bin/bash
# EGAPx v1.0 end to end on Ciona intestinalis (GCF_018327825.1, 140 Mb,
# `heldout` in benchmark/panel.tsv), proteins only, no RNA-seq, on one node.
# The measured row this produces goes into docs/cost-baseline/measured.tsv
# (docs/cost-baseline.md section 3), replacing the "not attempted" EGAPx cell.
#
# Submitted by gagarin through relay/bin/cluster-run.sh, which sets the Slurm
# limits (--cpus-per-task, --mem, --time, --output); nothing here overrides them.
#
#   cluster-run.sh submit --task T-gagarin-NNN --decision <id> \
#       --cpus 20 --mem 240 --gpus 0 --time 12:00:00 --name egapx-ciona \
#       docs/cost-baseline/cluster/egapx_ciona.sbatch <repo checkout> [scratch dir]
#
# Arguments: $1 = path of a checkout of nekrut/annotation at the commit named
# in the request (needs benchmark/fetch.py and benchmark/score.py);
# $2 = scratch directory (default $RELAY_CLUSTER_SCRATCH/T-human-009/egapx).
# Needs on PATH: python3 >= 3.11, git, nextflow (23.10 or later), and
# singularity or apptainer, plus /usr/bin/time (GNU time) and network access
# for the NCBI download, the egapx clone, the ncbi/egapx:1.0 image and EGAPx's
# support files. Everything that comes back to the relay is written to
# $SCRATCH/return/ and is small (text only).
set -euo pipefail

REPO=${1:?usage: egapx_ciona.sbatch <repo checkout> [scratch dir]}
SCRATCH=${2:-${RELAY_CLUSTER_SCRATCH:-$HOME/relay-scratch}/T-human-009/egapx}
CPUS=${SLURM_CPUS_PER_TASK:-20}
MEM_GB=${SLURM_MEM_PER_NODE:-245760}; MEM_GB=$((MEM_GB / 1024))
SPECIES=Ciona_intestinalis
TAXID=7719
EGAPX_TAG=v1.0
RET="$SCRATCH/return"
log() { echo "[$(date -u +%FT%TZ)] $*"; }
die() { log "FAIL: $*"; exit 1; }

mkdir -p "$SCRATCH" "$RET"; cd "$SCRATCH"
log "host=$(hostname) cpus=$CPUS mem_gb=$MEM_GB scratch=$SCRATCH repo=$REPO"
[ -f "$REPO/benchmark/fetch.py" ] || die "no benchmark/fetch.py under $REPO"
git -C "$REPO" rev-parse HEAD > "$RET/repo_commit.txt"
command -v nextflow >/dev/null || die "nextflow not on PATH"
command -v /usr/bin/time >/dev/null || die "GNU time not at /usr/bin/time"
# EGAPx's -e singularity calls the `singularity` binary; apptainer ships one
# under that name on most installs, and a shim covers the rest.
if ! command -v singularity >/dev/null; then
  command -v apptainer >/dev/null || die "neither singularity nor apptainer on PATH"
  mkdir -p "$SCRATCH/bin"; ln -sf "$(command -v apptainer)" "$SCRATCH/bin/singularity"
  export PATH="$SCRATCH/bin:$PATH"
fi
{ nextflow -version 2>&1 | grep -i version; singularity --version; python3 --version; } > "$RET/versions.txt" 2>&1 || true
{ grep -m1 'model name' /proc/cpuinfo; nproc; free -g | head -2; } > "$RET/hardware.txt" 2>&1 || true

# 1. Inputs, verified against NCBI's md5checksums.txt by fetch.py; sha256 recorded.
log "fetching $SPECIES"
python3 "$REPO/benchmark/fetch.py" --species "$SPECIES" --what fasta,gff --dest "$SCRATCH/panel"
FNA_GZ=$(ls "$SCRATCH/panel/$SPECIES"/*_genomic.fna.gz)
GFF_GZ=$(ls "$SCRATCH/panel/$SPECIES"/*_genomic.gff.gz)
sha256sum "$FNA_GZ" "$GFF_GZ" > "$RET/inputs.sha256"
FNA="${FNA_GZ%.gz}"; [ -s "$FNA" ] || gunzip -kf "$FNA_GZ"
GENOME_BP=$(grep -v '^>' "$FNA" | tr -d '\n' | wc -c)
log "genome_bp=$GENOME_BP"

# 2. EGAPx at a pinned tag, in a venv with its Python requirement (PyYAML).
if [ ! -d egapx ]; then git clone -q --depth 1 --branch "$EGAPX_TAG" https://github.com/ncbi/egapx.git egapx; fi
git -C egapx rev-parse HEAD > "$RET/egapx_commit.txt"
[ -d venv ] || python3 -m venv venv
# shellcheck disable=SC1091
. venv/bin/activate
pip install -q pyyaml
[ -f egapx/ui/requirements.txt ] && pip install -q -r egapx/ui/requirements.txt || true

# 3. Cap per-process CPU and memory requests so the Nextflow local executor
# does not refuse them on a 20-CPU node (EGAPx's labels go up to 32 CPUs and
# 128 GB; the README says 200 GB for large genomes). Every edit is logged.
MEM_CAP=$((MEM_GB - 10))
python3 - "$CPUS" "$MEM_CAP" egapx > "$RET/resource_caps.txt" <<'PY'
import re, sys, pathlib
cpus, mem, root = int(sys.argv[1]), int(sys.argv[2]), pathlib.Path(sys.argv[3])
cpu_re = re.compile(r"(\bcpus\s*[=:]\s*)(\d+)")
mem_re = re.compile(r"(\bmemory\s*[=:]\s*['\"]?)(\d+)(\s*\.?\s*GB)", re.I)
for p in list(root.rglob("*.config")) + list(root.rglob("*.nf")):
    s = p.read_text(); t = s
    t = cpu_re.sub(lambda m: m.group(1) + (str(cpus) if int(m.group(2)) > cpus else m.group(2)), t)
    t = mem_re.sub(lambda m: m.group(1) + (str(mem) if int(m.group(2)) > mem else m.group(2)) + m.group(3), t)
    if t != s:
        p.write_text(t)
        for a, b in zip(s.splitlines(), t.splitlines()):
            if a != b: print(f"{p}: {a.strip()} -> {b.strip()}")
PY

# 4. Input YAML: genome and taxid only. Proteins come from EGAPx's own
# taxon-matched NCBI sets (miniprot default); no reads, by design of the
# request (the panel row has no RNA-seq and the cost of interest is the
# pipeline floor).
cat > input.yaml <<YAML
genome: $FNA
taxid: $TAXID
YAML
cp input.yaml "$RET/input.yaml"

# 5. Nextflow trace through a private NXF_HOME so nothing outside scratch changes.
export NXF_HOME="$SCRATCH/nxf"; mkdir -p "$NXF_HOME"
cat > "$NXF_HOME/config" <<NXF
trace { enabled = true; overwrite = true; file = '$SCRATCH/trace.txt'
        fields = 'task_id,hash,name,status,exit,submit,start,complete,duration,realtime,%cpu,peak_rss,peak_vmem,rchar,wchar,cpus,memory' }
report { enabled = true; overwrite = true; file = '$SCRATCH/report.html' }
NXF
export NXF_SINGULARITY_CACHEDIR="$SCRATCH/sif"; mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# 6. The run, timed on the driver process. CPU time of the workers comes from
# the trace (they are children of Nextflow, so GNU time's user+sys covers them
# too; both are reported).
log "starting EGAPx"
T0=$(date +%s)
set +e
/usr/bin/time -v -o "$RET/egapx.time" \
  python3 egapx/ui/egapx.py input.yaml -e singularity -w "$SCRATCH/work" -o "$SCRATCH/out" -lc "$SCRATCH/cache" \
  > "$SCRATCH/egapx.stdout" 2> "$SCRATCH/egapx.stderr"
STATUS=$?
set -e
T1=$(date +%s)
log "EGAPx exit=$STATUS wall_s=$((T1 - T0))"
tail -n 200 "$SCRATCH/egapx.stderr" > "$RET/egapx.stderr.tail"
tail -n 200 "$SCRATCH/egapx.stdout" > "$RET/egapx.stdout.tail"
cp "$SCRATCH/trace.txt" "$RET/trace.txt" 2>/dev/null || true
ls -la "$SCRATCH/out" > "$RET/out_listing.txt" 2>/dev/null || true
du -sh "$SCRATCH/work" "$SCRATCH/cache" "$SCRATCH/out" > "$RET/disk_usage.txt" 2>/dev/null || true

# 7. Summary row (CPU-s from the trace and from GNU time, peak RSS, per-Mb).
python3 "$REPO/docs/cost-baseline/cluster/summarize_trace.py" \
  --trace "$SCRATCH/trace.txt" --time "$RET/egapx.time" --wall-s "$((T1 - T0))" \
  --genome-bp "$GENOME_BP" --species "$SPECIES" --exit "$STATUS" > "$RET/summary.tsv" || true
cat "$RET/summary.tsv" || true
[ "$STATUS" -eq 0 ] || die "EGAPx failed; logs in $RET"

# 8. Score with the benchmark's own scorer under a declaration (docs/benchmark.md 3.3).
PRED=$(find "$SCRATCH/out" -maxdepth 2 -name 'complete.genomic.gff*' | head -1)
[ -n "$PRED" ] || die "no complete.genomic.gff under $SCRATCH/out"
DECL="$RET/egapx-$SPECIES.yaml"
sed "s|^hardware:.*|hardware: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'), $CPUS cores allotted, $MEM_GB GB, no GPU|; s|^model:.*|model: EGAPx $EGAPX_TAG (git $(cat "$RET/egapx_commit.txt"), docker ncbi/egapx:1.0 via singularity, nextflow $(nextflow -version 2>&1 | grep -io 'version [0-9.]*' | head -1)), input.yaml with genome and taxid only|" \
  "$REPO/docs/cost-baseline/cluster/declarations/egapx-$SPECIES.yaml" > "$DECL"
python3 "$REPO/benchmark/score.py" --reference "$GFF_GZ" --prediction "$PRED" --species "$SPECIES" \
  --declaration "$DECL" --genome "$FNA_GZ" --out "$RET/egapx-$SPECIES.json"
gzip -c "$PRED" > "$SCRATCH/egapx-$SPECIES.gff.gz"
ls -la "$SCRATCH/egapx-$SPECIES.gff.gz" >> "$RET/out_listing.txt"
log "done; return files in $RET"
