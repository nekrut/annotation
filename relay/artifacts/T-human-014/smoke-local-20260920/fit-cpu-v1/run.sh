#!/bin/bash
# Candidate A, bounded local CPU fit v1: S. cerevisiae (dev chr I) + C. elegans (dev chr V), one thread pinned to core 2.
set -u
cd /tmp/wt014
S=/tmp/lenin-scratch/fit/fit-cpu-v1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "start=$(date -u +%FT%TZ) commit=$(git rev-parse HEAD) core=2 cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > $S/run.log
/tmp/lenin-venv/bin/python benchmark/leakage_check.py > $S/leakage_check.out 2>&1
echo "leakage_check exit=$?" >> $S/run.log
md5sum /tmp/lenin-scratch/sources/Saccharomyces_cerevisiae/*.gz /tmp/lenin-scratch/sources/Caenorhabditis_elegans/*.gz > $S/source_md5.txt
taskset -c 2 /usr/bin/time -v -o $S/train_time.txt \
  /tmp/lenin-venv/bin/python -u -m model.a.train train --config $S/config.json > $S/train.out 2>&1
echo "train exit=$? end=$(date -u +%FT%TZ)" >> $S/run.log
