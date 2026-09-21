#!/bin/bash
# AUGUSTUS 3.5.0 (bioconda) on the three S. pombe nuclear chromosomes, one process, pinned to core 2,
# the marx-0026 command (docs/cost-baseline.md section 3.2 row 1). Waits for the A runs to finish first.
set -u
cd /tmp/lenin-scratch/pombe/augustus
while ! grep -q '^end=' /tmp/lenin-scratch/pombe/A/run.log 2>/dev/null; do sleep 10; done
AUG=/tmp/lenin-scratch/augustus-env/bin/augustus
export AUGUSTUS_CONFIG_PATH=/tmp/lenin-scratch/augustus-env/config
echo "start=$(date -u +%FT%TZ) core=2 version=$($AUG --version 2>&1 | head -1) cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)" > run.log
taskset -c 2 /usr/bin/time -v -o augustus_time.txt $AUG --species=schizosaccharomyces_pombe --gff3=on --UTR=off nuclear.fna > augustus.gff3 2> augustus.err
echo "exit=$? end=$(date -u +%FT%TZ)" >> run.log
echo "genes=$(grep -c $'\tgene\t' augustus.gff3)" >> run.log
grep -E "User time|System time|Elapsed|Maximum resident" augustus_time.txt >> run.log
