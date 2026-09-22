#!/bin/bash
# Calibration sweep (section 3.7): the same fit-cpu-v3 best.pt, the same tiling and
# the same pooled section 3.5 matrix, decoded with the PWM multiplied by a scalar.
# $1 = the scale (e.g. 0.25); no hard mask, i.e. the accepted finite-score support,
# which is the default and the regime score-pwm-only/ measured at scale 1. Core 2.
set -u
cd /tmp/wt014
SCALE=${1:?scale}
TAG=$(printf '%03d' "$(python3 -c "print(round(float('$SCALE')*100))")")
S=/tmp/lenin-scratch/fit/score-scale-$TAG-v3
F=/tmp/lenin-scratch/fit/fit-cpu-v3
CK=$F/best.pt
PWM=/tmp/lenin-scratch/fit/pwm-v1/splicepwm.json
mkdir -p $S
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# Completion gate, inherited from run_pwm.sh (engels-0108 P3): taskset only sets
# affinity; require the fit's own recorded status before scoring.
grep -q '^train exit=0 ' $F/run.log || { echo "gate: no 'train exit=0' in $F/run.log"; exit 1; }
grep -q 'Exit status: 0' $F/train_time.txt || { echo "gate: /usr/bin/time did not record Exit status: 0"; exit 1; }
if pgrep -f 'model[.]a[.]train train --config .*fit-cpu-v3' > /dev/null; then
  echo "gate: a fit-cpu-v3 training process is still running"; exit 1
fi
[ -s $CK ] || { echo "gate: missing checkpoint $CK"; exit 1; }
# The matrix is an input of this sweep, not an output of it: every point must
# score the same declared pooled matrix, only the scalar differs.
[ "$(sha256sum $PWM | cut -d' ' -f1)" = a964e3e64c1ee3386fa421d7a9fd9630b29dbc3bb94fb08f41288e6bcb5a4f3a ] \
  || { echo "gate: splicepwm.json is not the declared pooled matrix"; exit 1; }
echo "gate: fit-cpu-v3 complete, checkpoint and declared PWM present (scale $SCALE -> $S)"
cp /tmp/lenin-scratch/fit/score-pwm-only-v3/seqids_chrI.txt /tmp/lenin-scratch/fit/score-pwm-only-v3/seqids_chrV.txt $S/
# Declaration: the scale-1 one with the scalar stated, so each scored point
# records which calibration produced it.
sed "s/67,325 junctions of the train split)/67,325 junctions of the train split, scaled by $SCALE)/" \
  /tmp/lenin-scratch/fit/score-pwm-only-v3/declaration.yaml > $S/declaration.yaml
grep -q "scaled by $SCALE" $S/declaration.yaml || { echo 'gate: declaration did not record the scale'; exit 1; }
run() {  # species seqid seqidsfile gff fna
  SP=$1; SEQ=$2; IDS=$3; GFF=$4; FNA=$5; N=${SEQ}_seg19_ov4096_w12288_float64
  taskset -c 2 /usr/bin/time -v -o $S/${N}_time.txt \
    /tmp/lenin-venv/bin/python -m model.a.train measure --config $F/config.json \
    --species $SP --seqid $SEQ --checkpoint $CK \
    --profile chromosome --window 12288 --overlap 4096 --segments 19 --dtype float64 \
    --splice-pwm $PWM --splice-pwm-scale $SCALE \
    --gff-out $S/$N.gff3 --json-out $S/measure_$N.json > $S/$N.out 2>&1
  M=$?
  echo "measure $SEQ exit=$M"
  # stalin-0113 P2: a failed measurement must not be scored.
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
echo "all workloads exit=0 (scale $SCALE)"
