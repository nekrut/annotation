#!/bin/bash
# Candidate A, frozen smoke checkpoint, S. pombe nuclear chromosomes, one thread pinned to core 2.
set -u
cd /tmp/wt014
O=/tmp/lenin-scratch/pombe/A
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "start=$(date -u +%FT%TZ) commit=$(git rev-parse HEAD) core=2 cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > $O/run.log
for SEQ in NC_003424.3 NC_003423.3 NC_003421.2; do
  for SEG in 1 19; do
    N=${SEQ}_seg${SEG}_ov4096_w12288
    taskset -c 2 /usr/bin/time -v -o $O/${N}_time.txt \
      /tmp/lenin-venv/bin/python -m model.a.train measure --config /tmp/lenin-scratch/pombe/config.json \
      --species Schizosaccharomyces_pombe --seqid $SEQ --checkpoint /tmp/lenin-scratch/pooledrun/best.pt \
      --profile chromosome --window 12288 --overlap 4096 --segments $SEG \
      --gff-out $O/$N.gff3 --json-out $O/measure_$N.json > $O/$N.out 2>&1
    echo "$N exit=$? $(date -u +%FT%TZ)" >> $O/run.log
  done
done
echo "end=$(date -u +%FT%TZ)" >> $O/run.log
