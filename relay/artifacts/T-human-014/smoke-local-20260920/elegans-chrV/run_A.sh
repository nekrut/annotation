#!/bin/bash
# Candidate A, frozen smoke checkpoint, C. elegans chr V (NC_003283.11, 20,924,180 bases), one thread pinned to core 2.
# Seam overlap 16,384 (8,192-base chain containment) for metazoan gene lengths; float64 rows of record, then float32.
set -u
cd /tmp/wt014
S=/tmp/lenin-scratch/elegans
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "start=$(date -u +%FT%TZ) commit=$(git rev-parse HEAD) core=2 cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > $S/A/run.log
SEQ=NC_003283.11
for RUN in "19 float64 A" "1 float64 A" "19 float32 A-f32"; do
  set -- $RUN; SEG=$1; DT=$2; O=$S/$3
  N=${SEQ}_seg${SEG}_ov16384_w12288_${DT}
  taskset -c 2 /usr/bin/time -v -o $O/${N}_time.txt \
    /tmp/lenin-venv/bin/python -m model.a.train measure --config $S/config.json \
    --species Caenorhabditis_elegans --seqid $SEQ --checkpoint /tmp/lenin-scratch/pooledrun/best.pt \
    --profile chromosome --window 12288 --overlap 16384 --segments $SEG --dtype $DT \
    --gff-out $O/$N.gff3 --json-out $O/measure_$N.json > $O/$N.out 2>&1
  echo "$N exit=$? $(date -u +%FT%TZ)" >> $S/A/run.log
done
echo "end=$(date -u +%FT%TZ)" >> $S/A/run.log
