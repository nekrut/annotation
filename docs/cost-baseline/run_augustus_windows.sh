#!/bin/bash
# Time AUGUSTUS on one long sequence split into K prediction windows, P at a time.
# For chromosomes of gigabase genomes, where one sequence would hold a core
# for longer than a bounded run allows. Each window is a separate AUGUSTUS
# process over the same FASTA with --predictionStart/--predictionEnd, so every
# process reads the whole sequence (a few seconds each) and genes spanning a
# window boundary are cut; both are disclosed in the row that uses this.
# usage: run_augustus_windows.sh <one-sequence FASTA> <augustus --species value> <K> <P> <work dir> [extra augustus args]
# Needs: augustus, /usr/bin/time (GNU time), bc.
set -u
FA=$1; PARAM=$2; K=$3; P=$4; W=$5; shift 5; EXTRA="$*"
rm -rf "$W"; mkdir -p "$W/parts"; cd "$W"
SEQ=$(head -1 "$FA" | sed 's/^>//; s/ .*//')
LEN=$(grep -v '^>' "$FA" | tr -d '\n' | wc -c)
NONN=$(grep -v '^>' "$FA" | tr -d '\n' | tr -d 'Nn' | wc -c)
STEP=$(( (LEN + K - 1) / K ))
: > windows.txt
for i in $(seq 0 $((K-1))); do
  s=$((i*STEP+1)); e=$(((i+1)*STEP)); [ $e -gt $LEN ] && e=$LEN
  echo "$i $s $e" >> windows.txt
done
echo "fasta=$FA seq=$SEQ param=$PARAM windows=$K nproc=$P extra=$EXTRA cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') start=$(date -u +%FT%TZ)" > run.log
T0=$(date +%s.%N)
xargs -P "$P" -L1 sh -c \
  "/usr/bin/time -v -o parts/w\$0.time augustus --species=$PARAM --gff3=on --genemodel=complete --UTR=off --predictionStart=\$1 --predictionEnd=\$2 $EXTRA $FA > parts/w\$0.gff3 2> parts/w\$0.err" < windows.txt
T1=$(date +%s.%N)
echo "wall_s=$(echo "$T1 - $T0" | bc)" >> run.log
echo "end=$(date -u +%FT%TZ)" >> run.log
# AUGUSTUS restarts gene numbering per invocation; prefix ids with the window.
(echo '##gff-version 3'; for i in $(seq 0 $((K-1))); do grep -v '^#' "parts/w$i.gff3" | sed "s/ID=/ID=w${i}_/; s/Parent=/Parent=w${i}_/"; done) > augustus.gff3
awk -F': ' '/User time/{u+=$2} /System time/{s+=$2} /Maximum resident/{if($2>m)m=$2} END{printf "sum_user_s=%.1f\nsum_sys_s=%.1f\nmax_part_rss_kb=%d\n", u, s, m}' parts/*.time >> run.log
echo "genes=$(grep -c $'\tgene\t' augustus.gff3)" >> run.log
echo "seq_bp=$LEN" >> run.log
echo "non_n_bp=$NONN" >> run.log
cat run.log
