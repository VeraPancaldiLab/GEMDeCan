#!/bin/bash

outdir=$1
indir=$2
for sample_number in $(find "$indir"/*.fastq.gz | awk -F'/' '{print $NF}'|  sed 's/_.*//'|tail -n+1 | sort -u ); do
  echo $sample_number
  for subgroup in R1 R2; do
    cat "$indir/$sample_number"_L*_"$subgroup"_001.fastq.gz > "$outdir/$sample_number"_"$subgroup.fastq.gz"
  done
done