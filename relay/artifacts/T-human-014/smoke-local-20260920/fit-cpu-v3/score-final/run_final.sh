#!/bin/bash
# Final scoring: fit-cpu-v3 best.pt on S. cerevisiae chr I and C. elegans chr V; core 2, after the fit has exited 0.
set -u
cd /tmp/wt014
S=/tmp/lenin-scratch/fit/score-final-v3
F=/tmp/lenin-scratch/fit/fit-cpu-v3
CK=$F/best.pt
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# Completion gate (engels-0108 P3): taskset only sets affinity. It does not wait
# for the fit, does not check that the fit succeeded, and does not make the core
# idle. Require the fit's own recorded status before scoring; an idle core stays
# an operator precondition, recorded in the README when the run happens.
grep -q '^train exit=0 ' $F/run.log || { echo "gate: no 'train exit=0' in $F/run.log"; exit 1; }
grep -q 'Exit status: 0' $F/train_time.txt || { echo "gate: /usr/bin/time did not record Exit status: 0"; exit 1; }
if pgrep -f 'model[.]a[.]train train --config .*fit-cpu-v3' > /dev/null; then
  echo "gate: a fit-cpu-v3 training process is still running"; exit 1
fi
[ -s $CK ] || { echo "gate: missing checkpoint $CK"; exit 1; }
echo "gate: fit-cpu-v3 complete, checkpoint present"
run() {  # species seqid seqidsfile gff fna
  SP=$1; SEQ=$2; IDS=$3; GFF=$4; FNA=$5; N=${SEQ}_seg19_ov4096_w12288_float64
  taskset -c 2 /usr/bin/time -v -o $S/${N}_time.txt \
    /tmp/lenin-venv/bin/python -m model.a.train measure --config /tmp/lenin-scratch/fit/fit-cpu-v3/config.json \
    --species $SP --seqid $SEQ --checkpoint $CK \
    --profile chromosome --window 12288 --overlap 4096 --segments 19 --dtype float64 \
    --gff-out $S/$N.gff3 --json-out $S/measure_$N.json > $S/$N.out 2>&1
  M=$?
  echo "measure $SEQ exit=$M"
  # stalin-0113 P2: a failed measurement must not be scored, and must not be
  # reported as success by the wrapper; a stale $N.gff3 from an earlier run
  # would otherwise be scored against the current checkpoint's declaration.
  [ $M -eq 0 ] || return $M
  taskset -c 2 /usr/bin/time -v -o $S/score_${N}_time.txt /tmp/lenin-venv/bin/python benchmark/score.py \
    --reference $GFF --genome $FNA --prediction $S/$N.gff3 --species $SP --declaration $S/declaration.yaml \
    --seqids $S/$IDS --out $S/score_$N.json > $S/score_$N.stdout 2>&1
  Q=$?
  echo "score $SEQ exit=$Q"
  return $Q
}
SC=/tmp/lenin-scratch/sources/Saccharomyces_cerevisiae; CE=/tmp/lenin-scratch/sources/Caenorhabditis_elegans
run Saccharomyces_cerevisiae NC_001133.9 seqids_chrI.txt $SC/GCF_000146045.2_R64_genomic.gff.gz $SC/GCF_000146045.2_R64_genomic.fna.gz || exit $?
run Caenorhabditis_elegans NC_003283.11 seqids_chrV.txt $(ls $CE/*.gff.gz) $(ls $CE/*.fna.gz) || exit $?
echo "all workloads exit=0"
