You can get exemple files in the following links :
* Fasta : `wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.22.fa.gz`
* GTF : `wget ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz`
* CDNA : `wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz`
* Reads : http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

You will need to rename the reads files, we suggest you use the `rename` command in bash like so :  
`rename 's/.read1/_R1/g' *.fastq.gz`  
`rename 's/.read2/_R2/g' *.fastq.gz`

Once this is done, copy the `samples.txt` and `config.yaml` files of this directory.
* `samples.txt` goes in `Data/`
* `config.yaml` shall be placed at the same level of the `Snakefile`

Make sure all data (GTF, Genome fasta, CDNA, samples.txt  and reads) are in the `data/` folder.
