You can get exemple files in the following links :
* Fasta : `wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.chromosome.22.fa.gz`
* GTF : `wget ftp://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr.gtf.gz`
* CDNA : `wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz`
* Reads : http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

You will need to rename the reads files, we suggest you use the `rename` command (you may need to install it : `sudo apt-get install rename`) in bash like so :  
`rename 's/.read1/_R1/g' *.fastq.gz`  
`rename 's/.read2/_R2/g' *.fastq.gz`

Once this is done, copy the `samples.txt` and `config.yaml` files of this directory.
* `samples.txt` goes in `Data/`
* `config.yaml` shall be placed at the same level of the `Snakefile`

Make sure all data (GTF, Genome fasta, CDNA, samples.txt  and reads) are in the `data/` folder.

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
