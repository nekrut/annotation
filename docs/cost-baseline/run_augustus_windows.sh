#!/bin/bash
# Time AUGUSTUS on one long sequence split into K prediction windows, P at a time.
# For chromosomes of gigabase genomes, where one sequence would hold a core
# for longer than a bounded run allows. Each window is a separate AUGUSTUS
# process over the same FASTA with --predictionStart/--predictionEnd, so every
# process reads the whole sequence (a few seconds each) and genes spanning a
# window boundary are cut; both are disclosed in the row that uses this.
# usage: run_augustus_windows.sh <one-sequence FASTA> <augustus --species value> <K> <P> <work dir> [extra augustus args]
# Needs: augustus, /usr/bin/time (GNU time), bc.
#
# Failure handling: the wrapper exits 1, and writes no augustus.gff3 and no
# success summary, unless every one of the K windows ran to completion with
# exit status 0 as recorded by GNU time. A window that predicts no gene is a
# success. The per-window record (windows.tsv: exit status, CPU, memory,
# genes) is what makes a run auditable without the genome.
set -euo pipefail
die() { echo "run_augustus_windows.sh: $*" >&2; exit 1; }
[ $# -ge 5 ] || die "usage: $0 <FASTA> <species> <K> <P> <work dir> [extra augustus args]"
FA=$1; PARAM=$2; K=$3; P=$4; W=$5; shift 5; EXTRA="$*"
for tool in augustus /usr/bin/time bc; do command -v "$tool" >/dev/null || die "$tool not found"; done
[ -s "$FA" ] || die "FASTA $FA is missing or empty"
[ "$K" -ge 1 ] && [ "$P" -ge 1 ] || die "K and P must be positive integers"
FA=$(readlink -f "$FA")
rm -rf "$W"; mkdir -p "$W/parts"; cd "$W"
SEQ=$(head -1 "$FA" | sed 's/^>//; s/ .*//')
[ -n "$SEQ" ] || die "no FASTA header in $FA"
[ "$(grep -c '^>' "$FA")" -eq 1 ] || die "$FA must hold exactly one sequence"
LEN=$(grep -v '^>' "$FA" | tr -d '\n' | wc -c)
NONN=$(grep -v '^>' "$FA" | tr -d '\n' | tr -d 'Nn' | wc -c)
[ "$LEN" -ge "$K" ] || die "sequence of $LEN bp cannot be split into $K windows"
STEP=$(( (LEN + K - 1) / K ))
: > windows.txt
for i in $(seq 0 $((K-1))); do
  s=$((i*STEP+1)); e=$(((i+1)*STEP)); [ $e -gt $LEN ] && e=$LEN
  echo "$i $s $e" >> windows.txt
done
CPU=$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 | sed 's/^ //' || true)
echo "fasta=$FA fasta_md5=$(md5sum "$FA" | cut -d' ' -f1) seq=$SEQ param=$PARAM windows=$K nproc=$P extra=$EXTRA cpu=${CPU:-unknown} start=$(date -u +%FT%TZ)" > run.log
T0=$(date +%s.%N)
set +e
xargs -P "$P" -L1 sh -c \
  "/usr/bin/time -v -o parts/w\$0.time augustus --species=$PARAM --gff3=on --genemodel=complete --UTR=off --predictionStart=\$1 --predictionEnd=\$2 $EXTRA $FA > parts/w\$0.gff3 2> parts/w\$0.err" < windows.txt
XRC=$?
set -e
T1=$(date +%s.%N)
# Audit every scheduled window before anything is combined: GNU time records
# "Exit status: N", or "Command terminated by signal N" when the child died.
printf 'window\tstart\tend\texit_status\tuser_s\tsys_s\tmax_rss_kb\tgenes\n' > windows.tsv
FAILED=0
while read -r i s e; do
  t="parts/w$i.time"; g="parts/w$i.gff3"
  status=missing; user=NA; sys=NA; rss=NA; genes=NA
  if [ -s "$t" ]; then
    if grep -q 'Command terminated by signal' "$t"; then
      status="signal:$(sed -n 's/.*terminated by signal //p' "$t" | head -1)"
    else
      status=$(awk -F': ' '/Exit status/{print $2}' "$t"); status=${status:-missing}
    fi
    user=$(awk -F': ' '/User time/{print $2}' "$t"); sys=$(awk -F': ' '/System time/{print $2}' "$t")
    rss=$(awk -F': ' '/Maximum resident/{print $2}' "$t")
  fi
  # grep -c prints 0 and exits 1 on a window with no gene: keep the 0.
  [ -f "$g" ] && genes=$(grep -c $'\tgene\t' "$g" || true)
  [ "$status" = 0 ] || FAILED=$((FAILED+1))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$i" "$s" "$e" "$status" "${user:-NA}" "${sys:-NA}" "${rss:-NA}" "$genes" >> windows.tsv
done < windows.txt
if [ "$XRC" -ne 0 ] || [ "$FAILED" -ne 0 ]; then
  echo "FAILED: xargs exit $XRC; $FAILED of $K windows did not complete with exit status 0 (see windows.tsv, parts/*.err)" | tee -a run.log >&2
  exit 1
fi
echo "wall_s=$(echo "$T1 - $T0" | bc)" >> run.log
echo "end=$(date -u +%FT%TZ)" >> run.log
# AUGUSTUS restarts gene numbering per invocation; prefix ids with the window.
# A window with no gene has only comment lines, so grep -v may match nothing.
(echo '##gff-version 3'; for i in $(seq 0 $((K-1))); do { grep -v '^#' "parts/w$i.gff3" || true; } | sed "s/ID=/ID=w${i}_/; s/Parent=/Parent=w${i}_/"; done) > augustus.gff3
awk -F': ' '/User time/{u+=$2} /System time/{s+=$2} /Maximum resident/{if($2>m)m=$2} END{printf "sum_user_s=%.1f\nsum_sys_s=%.1f\nmax_part_rss_kb=%d\n", u, s, m}' parts/*.time >> run.log
echo "genes=$(grep -c $'\tgene\t' augustus.gff3 || true)" >> run.log
echo "seq_bp=$LEN" >> run.log
echo "non_n_bp=$NONN" >> run.log
echo "windows_ok=$K/$K" >> run.log
cat run.log
