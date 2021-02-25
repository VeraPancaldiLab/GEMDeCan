#!/bin/bash
# $1 genome, $2 CDNA $3 threads $4 output
# Index genome for SALMON
grep "^>" <(gunzip -c "$1") | cut -d " " -f 1 >decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat "$2" "$1" >gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p "$3" --gencode -i "$4"
