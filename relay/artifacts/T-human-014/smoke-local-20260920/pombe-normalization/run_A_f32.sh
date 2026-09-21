#!/bin/bash
# Candidate A, frozen smoke checkpoint, S. pombe nuclear chromosomes, float32 decode, one thread pinned to core 2.
set -u
cd /tmp/wt014
O=/tmp/lenin-scratch/pombe/A-f32
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "start=$(date -u +%FT%TZ) commit=$(git rev-parse HEAD) core=2 cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > $O/run.log
for SEG in 19 1; do
for SEQ in NC_003424.3 NC_003423.3 NC_003421.2; do
    N=${SEQ}_seg${SEG}_ov4096_w12288_f32
    taskset -c 2 /usr/bin/time -v -o $O/${N}_time.txt \
      /tmp/lenin-venv/bin/python -m model.a.train measure --config /tmp/lenin-scratch/pombe/config.json \
      --species Schizosaccharomyces_pombe --seqid $SEQ --checkpoint /tmp/lenin-scratch/pooledrun/best.pt \
      --profile chromosome --window 12288 --overlap 4096 --segments $SEG --dtype float32 \
      --gff-out $O/$N.gff3 --json-out $O/measure_$N.json > $O/$N.out 2>&1
    echo "$N exit=$? $(date -u +%FT%TZ)" >> $O/run.log
  done
done
echo "end=$(date -u +%FT%TZ)" >> $O/run.log
