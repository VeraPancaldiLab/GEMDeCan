#!/bin/bash

# Index genome for Kallisto
kallisto index -i data/genome/transcript.idx data/genome/homo_cdna.fa.gz

# Index genome for STAR
STAR --runThreadN 8 \ # change the number of threads you can use
    --runMode genomeGenerate \
    --genomeDir data/star \
    --genomeFastaFiles data/genome/homo_dna.fasta


# Index genome for SALMON
grep "^>" <(gunzip -c homo_dna.fa.g) | cut -d " " -f 1 > decoys.txt
grep "^>" <(gunzip -c homo_dna.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat homo_cdna.fa.gz homo_dna.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 --index salmon_index --gencode
# -p 12 refers to the number of threads to be used (here 12), feel free to change it