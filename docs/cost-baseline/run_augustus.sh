#!/bin/bash
# Time AUGUSTUS on one panel genome, one process per sequence, P at a time.
# usage: run_augustus.sh <Species_name> <augustus --species value> <P> <panel dir> <work dir> [extra augustus args]
# Needs: augustus, /usr/bin/time (GNU time), bc. Panel files come from benchmark/fetch.py.
set -u
SP=$1; PARAM=$2; P=$3; PANEL=$4; W=$5; shift 5; EXTRA="$*"
rm -rf "$W"; mkdir -p "$W/parts"; cd "$W"
zcat "$PANEL/$SP"/*_genomic.fna.gz > genome.fna
awk '/^>/{n=substr($1,2); f="parts/"n".fa"} {print > f}' genome.fna
ls parts/*.fa | sed 's#parts/##; s#\.fa$##' > seqs.txt
echo "species=$SP param=$PARAM nproc=$P extra=$EXTRA cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //') start=$(date -u +%FT%TZ)" > run.log
T0=$(date +%s.%N)
xargs -P "$P" -I{} sh -c \
  "/usr/bin/time -v -o parts/{}.time augustus --species=$PARAM --gff3=on --genemodel=complete --UTR=off $EXTRA parts/{}.fa > parts/{}.gff3 2> parts/{}.err" < seqs.txt
T1=$(date +%s.%N)
echo "wall_s=$(echo "$T1 - $T0" | bc)" >> run.log
echo "end=$(date -u +%FT%TZ)" >> run.log
# AUGUSTUS restarts gene numbering per invocation; prefix ids with the sequence name.
(echo '##gff-version 3'; for n in $(cat seqs.txt); do grep -v '^#' "parts/$n.gff3" | sed "s/ID=/ID=${n}_/; s/Parent=/Parent=${n}_/"; done) > augustus.gff3
awk -F': ' '/User time/{u+=$2} /System time/{s+=$2} /Maximum resident/{if($2>m)m=$2} END{printf "sum_user_s=%.1f\nsum_sys_s=%.1f\nmax_part_rss_kb=%d\n", u, s, m}' parts/*.time >> run.log
echo "genes=$(grep -c $'\tgene\t' augustus.gff3)" >> run.log
echo "genome_bp=$(grep -v '^>' genome.fna | tr -d '\n' | wc -c)" >> run.log
echo "nseq=$(wc -l < seqs.txt)" >> run.log
cat run.log
