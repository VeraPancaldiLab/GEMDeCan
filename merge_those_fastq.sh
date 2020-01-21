#!/bin/bash

outdir=$1

function join_by { local IFS="$1"; shift; echo "$*"; }
SAMPLE_SUFFIX="_001"
EXTENSION_SUFFIX=".fastq.gz"


for sample_number in $(find $outdir/bcl2raw/*.fastq.gz | awk -F'/' '{print $NF}'| cut -d'-' -f3 | cut -d'_' -f1 | sort  -u ); do
  echo $sample_number
  BASE_FORMAT="$(find $outdir/bcl2raw/*.fastq.gz | awk -F'/' '{print $NF}'| cut -d'-' -f 1-2 | sort -u)-${sample_number}" 
  LINE_FORMAT="L00\d"
  for r_number in {1..2}; do
    R_FORMAT="${BASE_FORMAT}_${LINE_FORMAT}_R${r_number}"
    FILE_TEMPLATE="${R_FORMAT}${SAMPLE_SUFFIX}${EXTENSION_SUFFIX}"
    echo $FILE_TEMPLATE
    FILES="$(find $1 | grep -P $FILE_TEMPLATE | sort -u)"
    echo $FILES
    JOINED_FILES=$(join_by ' ' $FILES)
   cat ${JOINED_FILES} > "$outdir/raw2merge/${BASE_FORMAT}_R${r_number}${EXTENSION_SUFFIX}"
	done
done