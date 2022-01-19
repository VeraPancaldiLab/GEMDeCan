You can get exemple files in the following links :
* Genome Fasta : `wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.22.fa.gz`
* GTF : `wget ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz`
* CDNA : `wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz`
* Reads : 
`wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar`
`tar xvf HBR_UHR_ERCC_ds_5pc.tar`

Here, you will need to rename the reads files - as they don't end by 'R1' and 'R2', we suggest to use the `rename` command in Bash.
You may need to install it : `sudo apt-get install rename` or 'conda install rename'.
And on this example dataset :
`rename 's/.read1/_R1/g' *.fastq.gz`
`rename 's/.read2/_R2/g' *.fastq.gz`

Once this is done, copy the `samples.txt` and `config.yaml` files of this directory.
* `samples.txt` should be located in the `data/` folder located at the level of the `Snakefile`
* `config.yaml` should be located at the level of the `Snakefile`

To reproduce this example based on the config.yaml, make sure all other source data (GTF, Genome Fasta, CDNA, samples.txt and reads) are also in the same `data/` folder.

In the end you should end up with a `data/` folder looking like this :
```
data/
├── HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22_R1.fastq.gz
├── HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22_R2.fastq.gz
├── HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22_R1.fastq.gz
├── HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22_R2.fastq.gz
├── HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22_R1.fastq.gz
├── HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22_R2.fastq.gz
├── Homo_sapiens.GRCh38.103.chr.gtf
├── Homo_sapiens.GRCh38.cdna.all.fa.gz
├── Homo_sapiens.GRCh38.dna_rm.chromosome.22.fa
├── samples.txt
├── signatures
│   ├── BPprom_314.txt
│   ├── BPRNACan_32.txt
│   ├── BPRNACan3DMet.txt
│   ├── BPRNACan_52.txt
│   ├── BPRNACanProMet3D.txt
│   ├── BPRNACan.txt
│   └── MCPcounter
│       └── MCPcounter-genes.txt
├── TruSeqAdapt.fa
├── UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22_R1.fastq.gz
├── UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22_R2.fastq.gz
├── UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22_R1.fastq.gz
├── UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22_R2.fastq.gz
├── UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22_R1.fastq.gz
└── UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22_R2.fastq.gz
```
