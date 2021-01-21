#!/bin/bash
# $2 input, $3 output

if [ "$1" == "kallisto" ]; then
    # Index genome for Kallisto
    kallisto index -i "$3" "$2"
elif [ "$1" == "star" ]; then
    # Index genome for STAR
    STAR --runThreadN "$4" \ # change the number of threads you can use
    --runMode genomeGenerate \
        --genomeDir "$3" \
        --genomeFastaFiles "$2"
elif [ "$1" == "salmon" ]; then
    # Index genome for SALMON
    grep "^>" <(gunzip -c "$2") | cut -d " " -f 1 >decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat "$3" "$2" >gentrome.fa.gz
    salmon index -t gentrome.fa.gz -d decoys.txt -p "$4" --index salmon_index --gencode
elif [ "$1" == "-h" ]; then
    echo "Usage : $(basename "$0") [method] [input] [output]"
else
    echo "Please chose a method to compute the index."
    echo "Usage : $(basename "$0") [method] [input] [output] [threads]"
    echo "Usage for salmon: $(basename "$0") salmon [input_genome] [input_CDNA] [threads]"
fi
