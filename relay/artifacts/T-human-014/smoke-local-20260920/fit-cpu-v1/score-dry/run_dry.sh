#!/bin/bash
# Scoring-pipeline dry run: interim best.pt (step 300 of fit-cpu-v1) on S. cerevisiae chr I; core 3 (the fit holds core 2).
set -u
cd /tmp/wt014
S=/tmp/lenin-scratch/fit/score-dry
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
SEQ=NC_001133.9; N=${SEQ}_seg19_ov4096_w12288_float64
taskset -c 3 /usr/bin/time -v -o $S/${N}_time.txt \
  /tmp/lenin-venv/bin/python -m model.a.train measure --config /tmp/lenin-scratch/fit/fit-cpu-v1/config.json \
  --species Saccharomyces_cerevisiae --seqid $SEQ --checkpoint $S/best_interim.pt \
  --profile chromosome --window 12288 --overlap 4096 --segments 19 --dtype float64 \
  --gff-out $S/$N.gff3 --json-out $S/measure_$N.json > $S/$N.out 2>&1
echo "measure exit=$?"
taskset -c 3 /usr/bin/time -v -o $S/score_time.txt /tmp/lenin-venv/bin/python benchmark/score.py \
  --reference /tmp/lenin-scratch/sources/Saccharomyces_cerevisiae/GCF_000146045.2_R64_genomic.gff.gz \
  --genome /tmp/lenin-scratch/sources/Saccharomyces_cerevisiae/GCF_000146045.2_R64_genomic.fna.gz \
  --prediction $S/$N.gff3 --species Saccharomyces_cerevisiae --declaration $S/declaration.yaml \
  --seqids $S/seqids_chrI.txt --out $S/score_$N.json > $S/score.out 2>&1
echo "score exit=$?"
