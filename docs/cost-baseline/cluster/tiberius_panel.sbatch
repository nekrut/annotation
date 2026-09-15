#!/bin/bash
# Tiberius 2.0.7 (the image pinned in benchmark/validation/) on the panel
# species its `vertebrates`, `insecta` and `fungi` models cover, one GPU,
# each species timed under GNU time with GPU memory sampled, then scored with
# benchmark/score.py under a declaration. Produces the measured Tiberius rows
# for docs/cost-baseline/measured.tsv (docs/cost-baseline.md section 3).
#
# Submitted by gagarin through relay/bin/cluster-run.sh, which sets the Slurm
# limits; the request asks for two jobs so the human genome does not share a
# 24 h window with the other ten:
#
#   cluster-run.sh submit ... --cpus 8 --mem 128 --gpus 1 --time 24:00:00 --name tiberius-ten \
#       docs/cost-baseline/cluster/tiberius_panel.sbatch <repo> <scratch> \
#       Saccharomyces_cerevisiae Schizosaccharomyces_pombe Neurospora_crassa \
#       Drosophila_melanogaster Apis_mellifera Takifugu_rubripes Gallus_gallus \
#       Danio_rerio Xenopus_tropicalis Mus_musculus
#   cluster-run.sh submit ... --cpus 8 --mem 128 --gpus 1 --time 24:00:00 --name tiberius-human \
#       docs/cost-baseline/cluster/tiberius_panel.sbatch <repo> <scratch> Homo_sapiens
#
# Arguments: $1 = checkout of nekrut/annotation at the commit named in the
# request; $2 = scratch directory; $3.. = species (panel.tsv names), run in
# the order given. Needs: python3 >= 3.11, apptainer or singularity with
# --nv, nvidia-smi, /usr/bin/time, network for NCBI, Docker Hub and the
# Greifswald weights server. Each species writes its return files as soon as
# it finishes, so a job cut off by the wall clock still returns the rest.
set -euo pipefail

REPO=${1:?usage: tiberius_panel.sbatch <repo checkout> <scratch dir> <species...>}
SCRATCH=${2:?scratch dir}
shift 2
IMAGE="docker://larsgabriel23/tiberius@sha256:2c3bddda32cc621b805de40dc0395942cb5dd8a5766b1fe6c98fba845740f9bd"
SIF="$SCRATCH/tiberius_2.0.7.sif"
BATCH=${TIBERIUS_BATCH:-8}     # README: RTX 3090 25 GB -> 8; the A5000 has 24 GB
SEQ_LEN=400050                  # as in benchmark/validation/tiberius-Takifugu_rubripes.yaml
RET="$SCRATCH/return"
log() { echo "[$(date -u +%FT%TZ)] $*"; }
die() { log "FAIL: $*"; exit 1; }
[ $# -gt 0 ] || die "no species given"
model_for() {
  case "$1" in
    Homo_sapiens|Mus_musculus|Xenopus_tropicalis|Danio_rerio|Gallus_gallus|Takifugu_rubripes) echo vertebrates ;;
    Apis_mellifera|Drosophila_melanogaster) echo insecta ;;
    Neurospora_crassa|Schizosaccharomyces_pombe|Saccharomyces_cerevisiae) echo fungi ;;
    *) return 1 ;;
  esac
}

mkdir -p "$SCRATCH" "$RET"; cd "$SCRATCH"
CTR=$(command -v apptainer || command -v singularity || true)
[ -n "$CTR" ] || die "neither apptainer nor singularity on PATH"
command -v nvidia-smi >/dev/null || die "nvidia-smi not on PATH"
command -v /usr/bin/time >/dev/null || die "GNU time not at /usr/bin/time"
git -C "$REPO" rev-parse HEAD > "$RET/repo_commit.txt"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv > "$RET/gpu.txt"
{ grep -m1 'model name' /proc/cpuinfo; nproc; free -g | head -2; "$CTR" --version; } > "$RET/hardware.txt" 2>&1 || true
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader | head -1)
CPU_NAME=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')
for sp in "$@"; do model_for "$sp" >/dev/null || die "no Tiberius model mapped for $sp"; done

# Image, pulled once by digest.
[ -s "$SIF" ] || "$CTR" pull "$SIF" "$IMAGE"
sha256sum "$SIF" > "$RET/image.sha256"
# Weights: the vertebrates set ships in the image; insecta and fungi are
# fetched by Tiberius on first use from the URLs in its model_cfg. A writable
# tmpfs overlay lets that download land inside the container; the tarballs
# are also fetched here once so their checksums are on record.
mkdir -p weights
for u in https://bioinf.uni-greifswald.de/bioinf/tiberius/models/vertebrates_weights.tar.gz \
         https://bioinf.uni-greifswald.de/bioinf/tiberius/models/insecta_unmasked_v2.tar.gz \
         https://bioinf.uni-greifswald.de/bioinf/tiberius/models/fungi.tar.gz; do
  f="weights/$(basename "$u")"; [ -s "$f" ] || curl -sSL -o "$f" "$u" || log "warning: could not fetch $u"
done
sha256sum weights/*.tar.gz > "$RET/weights.sha256" 2>/dev/null || true

printf 'species\tmodel_cfg\tgenome_mb\texit\twall_s\tcpu_user_s\tcpu_sys_s\tpeak_rss_mb\tpeak_gpu_mem_mb\tgenes\tbatch_size\n' > "$RET/summary.tsv"
for SP in "$@"; do
  CFG=$(model_for "$SP")
  W="$SCRATCH/$SP"; mkdir -p "$W"
  log "== $SP ($CFG)"
  python3 "$REPO/benchmark/fetch.py" --species "$SP" --what fasta,gff --dest "$SCRATCH/panel"
  FNA_GZ=$(ls "$SCRATCH/panel/$SP"/*_genomic.fna.gz); GFF_GZ=$(ls "$SCRATCH/panel/$SP"/*_genomic.gff.gz)
  sha256sum "$FNA_GZ" "$GFF_GZ" >> "$RET/inputs.sha256"
  FNA="$W/genome.fna"; [ -s "$FNA" ] || zcat "$FNA_GZ" > "$FNA"
  GENOME_BP=$(grep -v '^>' "$FNA" | tr -d '\n' | wc -c)
  # GPU memory sampler, every 5 s, killed when the run ends.
  nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 5 > "$W/gpu_samples.csv" &
  SMI=$!
  T0=$(date +%s.%N)
  set +e
  /usr/bin/time -v -o "$W/tiberius.time" \
    "$CTR" exec --nv --writable-tmpfs -B "$SCRATCH:$SCRATCH" "$SIF" \
      python3 /opt/Tiberius/tiberius.py --genome "$FNA" --out "$W/tiberius.gff3" \
        --model_cfg "$CFG" --batch_size "$BATCH" --seq_len "$SEQ_LEN" --no_softmasking \
      > "$W/tiberius.stdout" 2> "$W/tiberius.stderr"
  STATUS=$?
  set -e
  T1=$(date +%s.%N)
  kill "$SMI" 2>/dev/null || true; wait "$SMI" 2>/dev/null || true
  WALL=$(python3 -c "print(f'{$T1 - $T0:.1f}')")
  USER_S=$(awk -F': ' '/User time/{print $2}' "$W/tiberius.time"); SYS_S=$(awk -F': ' '/System time/{print $2}' "$W/tiberius.time")
  RSS_MB=$(awk -F': ' '/Maximum resident/{printf "%d", $2/1024}' "$W/tiberius.time")
  GPU_MB=$(awk -F', ' 'NR>1{gsub(/ MiB/,"",$2); if($2+0>m)m=$2+0} END{print m+0}' "$W/gpu_samples.csv")
  GENES=$( [ -s "$W/tiberius.gff3" ] && grep -c $'\tgene\t' "$W/tiberius.gff3" || echo 0 )
  printf '%s\t%s\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SP" "$CFG" "$(python3 -c "print($GENOME_BP/1e6)")" "$STATUS" "$WALL" "$USER_S" "$SYS_S" "$RSS_MB" "$GPU_MB" "$GENES" "$BATCH" >> "$RET/summary.tsv"
  tail -n 100 "$W/tiberius.stderr" > "$RET/tiberius-$SP.stderr.tail"; cp "$W/tiberius.time" "$RET/tiberius-$SP.time"
  log "$SP exit=$STATUS wall_s=$WALL user_s=$USER_S rss_mb=$RSS_MB gpu_mb=$GPU_MB genes=$GENES"
  [ "$STATUS" -eq 0 ] || { log "skipping score for $SP"; continue; }
  DECL="$RET/tiberius-$SP.yaml"
  sed "s|^hardware:.*|hardware: $CPU_NAME, ${SLURM_CPUS_PER_TASK:-8} cores allotted, ${SLURM_MEM_PER_NODE:-?} MB, 1x $GPU_NAME $GPU_MEM|; s|batch size [0-9]*|batch size $BATCH|" \
    "$REPO/docs/cost-baseline/cluster/declarations/tiberius-$SP.yaml" > "$DECL"
  EXTRA=""; [ "$SP" = Tetrahymena_thermophila ] && EXTRA="--genetic-code 6"
  python3 "$REPO/benchmark/score.py" --reference "$GFF_GZ" --prediction "$W/tiberius.gff3" --species "$SP" \
    --declaration "$DECL" --genome "$FNA_GZ" $EXTRA --out "$RET/tiberius-$SP.json" || log "scoring failed for $SP"
  gzip -kf "$W/tiberius.gff3"
  rm -f "$FNA"
done
cat "$RET/summary.tsv"
log "done; return files in $RET"
