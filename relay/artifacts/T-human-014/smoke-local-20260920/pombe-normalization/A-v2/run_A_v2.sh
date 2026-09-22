#!/bin/bash
# Candidate A, fit-cpu-v2 best.pt (step 900), S. pombe nuclear chromosomes, 19 segments per strand, one thread pinned to core 2.
set -u
cd /tmp/wt014
O=/tmp/lenin-scratch/pombe/A-v2
CK=/tmp/lenin-scratch/fit/fit-cpu-v2/best.pt
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "start=$(date -u +%FT%TZ) commit=$(git rev-parse HEAD) core=2 checkpoint_sha256=$(sha256sum $CK | cut -d' ' -f1) cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > $O/run.log
for SEQ in NC_003424.3 NC_003423.3 NC_003421.2; do
  SEG=19; N=${SEQ}_seg${SEG}_ov4096_w12288
  taskset -c 2 /usr/bin/time -v -o $O/${N}_time.txt \
    /tmp/lenin-venv/bin/python -m model.a.train measure --config /tmp/lenin-scratch/pombe/config.json \
    --species Schizosaccharomyces_pombe --seqid $SEQ --checkpoint $CK \
    --profile chromosome --window 12288 --overlap 4096 --segments $SEG \
    --gff-out $O/$N.gff3 --json-out $O/measure_$N.json > $O/$N.out 2>&1
  echo "$N exit=$? $(date -u +%FT%TZ)" >> $O/run.log
done
sha256sum $O/*.gff3 > $O/gff3_sha256.txt
echo "end=$(date -u +%FT%TZ)" >> $O/run.log
