#!/bin/bash
# $1 genome, $2 CDNA $3 threads $4 output
# Index genome for SALMON
mkdir -p ${4}
grep "^>" <(gunzip -c "$1") | cut -d " " -f 1 >${4}/decoys.txt
sed -i.bak -e 's/>//g' ${4}/decoys.txt
cat ${2} ${1} >${4}/gentrome.fa.gz
salmon index -t ${4}/gentrome.fa.gz -d ${4}/decoys.txt -p ${3} --gencode  -i ${4}
