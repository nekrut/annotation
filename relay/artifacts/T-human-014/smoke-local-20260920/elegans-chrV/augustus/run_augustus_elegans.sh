#!/bin/bash
# AUGUSTUS 3.5.0 (bioconda) on C. elegans chr V, one process pinned to core 4 (runs alongside the A runs on core 2),
# the marx-0026 command with --species=caenorhabditis.
set -u
cd /tmp/lenin-scratch/elegans/augustus
AUG=/tmp/lenin-scratch/augustus-env/bin/augustus
export AUGUSTUS_CONFIG_PATH=/tmp/lenin-scratch/augustus-env/config
echo "start=$(date -u +%FT%TZ) core=4 version=$($AUG --version 2>&1 | head -1) cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > run.log
taskset -c 4 /usr/bin/time -v -o augustus_time.txt $AUG --species=caenorhabditis --gff3=on --UTR=off chrV.fna > augustus.gff3 2> augustus.err
echo "exit=$? end=$(date -u +%FT%TZ)" >> run.log
echo "genes=$(grep -c $'\tgene\t' augustus.gff3)" >> run.log
grep -E "User time|System time|Elapsed|Maximum resident" augustus_time.txt >> run.log
