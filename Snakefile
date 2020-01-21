#######################################
### Common part for pipeline ###
#######################################
# Julien Pernet 2020 for Pancaldi lab - CRCT Team 21 - INSERM
# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################

from snakemake.utils import validate

configfile: "config.yaml"
validate(config, "schema.yaml")

# Define paths
OUTDIR = config["Output_Directory"]
INDIR = config["Input_Directory"]
SAMPLESHEET = config["Sample_Sheet"]
THREADS = config["THREADS"]
INDEXK = config["Index_kallisto"]
INDEXS = config["Index_salmon"]
INDEXSTAR = config["Index_STAR"]
GENOME = config["Genome"]
GTF = config["GTF"]
ADAPTER = config["Adapter"]

OUTbcl = OUTDIR+"/bcl2raw"
OUTmerge = OUTDIR+"/raw2merge"
OUTfastqc = OUTDIR+"/fastqc_raw"
OUTfastqc2 = OUTDIR+"/fastqc_cutadapter"
OUTmultiqc = OUTDIR+"/multiqc_raw"
OUTmultiqc2 = OUTDIR+"/multiqc_after_cutadapter"
OUTcut = OUTDIR+"/data_after_cutadapter"
QUANTIF = OUTDIR+"Quantification"

## Outputs
rule all:
    input:
        OUTmultiqc+"/multiqc_report.html",
        OUTmultiqc2+"/multiqc_report.html",
        expand(QUANTIF+"/{samples}/quantif.txt", samples=config["Samples"])

## Converts base call (.BCL) files into FASTQ
if config["Convert_bcl2fastq"] == "yes" :
    rule bcl2fastq:
        input:
            INDIR = INDIR,
            SHEET = SAMPLESHEET
        output:
            touch("bcl2fastq.done")
        params:
            OUTbcl
        message:
            "Converting BCL2 Illumina files to fastQ"
        conda:
            "Tools/bcl2fastq.yaml"
        threads: THREADS
        shell:
            """
            bcl2fastq -w {threads} -R {input.INDIR} -o {params} --sample-sheet {input.SHEET}
            rm {params}/Undetermined*
            rename "s/_S[0-9]+//" {params}/*.fastq.gz
            """

##  Merge into one fastq file from different lanes
rule merging:
    input:
        "bcl2fastq.done"
    output: 
        OUTmerge+"/{samples}_R1.fastq.gz",
        OUTmerge+"/{samples}_R2.fastq.gz"
    params:
        OUT = OUTDIR
    message:
        "Merging files from lanes"
    shell:
        "bash merge_those_fastq.sh {params.OUT}"

## Quality control for raw fastq data
rule fastqc1:
    input:
        OUTmerge+"/{samples}.fastq.gz"
    output:
        html = OUTfastqc+"/{samples}_fastqc.html",
        zip = OUTfastqc+"/{samples}_fastqc.zip"
    threads: THREADS
    message:
        "Quality control before trimming"
    params:
        "-t 8"
    wrapper:
        "0.47.0/bio/fastqc"

rule multiqc1:
    input:
        expand(OUTfastqc+"/{sample}_R1_fastqc.html", sample=config["Samples"]),
        expand(OUTfastqc+"/{sample}_R2_fastqc.html", sample=config["Samples"])
    output:
        MainOut = OUTmultiqc+"/multiqc_report.html"
    wrapper:
        "0.47.0/bio/multiqc"


if config["Trim_with"] == "Trimmomatic" :
    ## Read trimming by Trimmomatic (Paired-End)
    rule Trimmomatic:
        input:
            r1 = OUTmerge+"/{samples}_R1.fastq.gz",
            r2 = OUTmerge+"/{samples}_R2.fastq.gz"
        output:
            r1 = OUTcut+"/{samples}_R1.fastq.gz",
            r1_unpaired = OUTcut+"/{samples}_R1.unpaired.fastq.gz",
            r2 = OUTcut+"/{samples}_R2.fastq.gz",
            r2_unpaired = OUTcut+"/{samples}_R2.unpaired.fastq.gz"
        threads: THREADS
        message:
            "Trimming using Trimmomatic"
        params:
            trimmer = ["TRAILING:20", "LEADING:20", "MINLEN:36", "CROP:10000", "ILLUMINACLIP:"+ADAPTER+":2:30:10"],
            compression_level="-9",
            extra = "-phred33"
        wrapper:
            "0.47.0/bio/trimmomatic/pe"

elif config["Trim_with"] == "Trimgalore" :
    ## Read trimming by Trim-galore (Paired-end)
    rule trimgalore:
        input:
            [expand(OUTmerge+"/{sample}_R1.fastq.gz", sample=config["Samples"]),
            expand(OUTmerge+"/{sample}_R2.fastq.gz", sample=config["Samples"])]
        output:
            OUTcut+"/{sample}_R1_val_1.fq.gz",
            OUTcut+"/{sample}_R1.fastq.gz_trimming_report.txt",
            OUTcut+"/{sample}_R2_val_2.fq.gz",
            OUTcut+"/{sample}_R2.fastq.gz_trimming_report.txt"
        params:
            extra= '--phred33 --illumina --paired --quality 20 --length 36'
        threads: THREADS
        message:
            "Trimming using Trim-Galore"
        wrapper:
            "0.47.0/bio/trim_galore/pe"

    rule rename:
        input:
            R1 = expand(OUTcut+"/{sample}_R1_val_1.fq.gz", sample=config["Samples"]),
            R2 = expand(OUTcut+"/{sample}_R2_val_2.fq.gz", sample=config["Samples"])
        output:
            R1out = OUTcut+"/{sample}_R1.fastq.gz",
            R2out = OUTcut+"/{sample}_R2.fastq.gz"
        shell:
            """
            mv -f {input.R1} {output.R1out}
            mv -f {input.R2} {output.R2out}
            """

## Quality control after trimming
rule fastqc2:
    input:
        OUTcut+"/{samples}.fastq.gz"
    output:
        html = OUTfastqc2+"/{samples}_fastqc.html",
        zip = OUTfastqc2+"/{samples}_fastqc.zip"
    threads: THREADS
    message:
        "Quality control after trimming"
    params:
        "-t 8"
    wrapper:
        "0.47.0/bio/fastqc"
rule multiqc2:
    input:
        expand(OUTfastqc2+"/{sample}_R1_fastqc.html", sample=config["Samples"]),
        expand(OUTfastqc2+"/{sample}_R2_fastqc.html", sample=config["Samples"])
    output:
        OUTmultiqc2+"/multiqc_report.html"
    wrapper:
        "0.47.0/bio/multiqc"

# Quantification
if config["Quantification_with"] == "kallisto" :
    rule kallisto_quant:
        input:
            R1 = expand(OUTcut+"/{sample}_R1.fastq.gz", sample=config["Samples"]),
            R2 = expand(OUTcut+"/{sample}_R2.fastq.gz", sample=config["Samples"]),
            INDEXK = INDEXK,
        threads: THREADS
        output:
            QUANTIF+"/{sample}/abundance.tsv",
            QUANTIF+"/{sample}/abundance.h5"
        params:
            OUTDIRE = QUANTIF+"/{sample}"
        message:
            "Quantification with Kallisto"
        conda:
            "Tools/kallisto.yaml"
        shell:
            "kallisto quant -t {threads} -i {input.INDEXK} -b 30 "
            "-o {params.OUTDIRE} "
            "{input.R1} {input.R2}"

    rule quant_to_gene:
        input:
            QUANTIF+"/{sample}/abundance.h5"
        output:
            QUANTIF+"/{sample}/quantif.txt"
        params:
            QUANTIF+"/{sample}"
        script:
            "quant_for_kallisto.R"
    
elif config["Quantification_with"] == "salmon" :
    rule salmon:
        input:
            r1 = expand(OUTcut+"/{sample}_R1.fastq.gz", sample=config["Samples"]),
            r2 = expand(OUTcut+"/{sample}_R2.fastq.gz", sample=config["Samples"]),
            index = INDEXS
        output:
            quant = QUANTIF+"/{sample}/quant.sf",
            lib = QUANTIF+"/{sample}/lib_format_counts.json"
        params:
            DIR = QUANTIF+"/{sample}",
            libtype ="A",
            extra=" --validateMappings"
        threads: THREADS
        message:
            "Quantification with Salmon"
        conda:
            "Tools/salmon.yaml"
        shell:
            "salmon quant -i {input.index} -l {params.libtype} "
            "-1 {input.r1} -2 {input.r2} "
            "-o {params.DIR} "
            "-p {threads} --validateMappings"
    
    rule salmon_quant:
        input:
            QUANTIF+"/{sample}/quant.sf"
        output:
            QUANTIF+"/{sample}/quantif.txt"
        params:
            QUANTIF+"/{sample}"
        script:
            "quant_for_salmon.R"

elif config["Quantification_with"] == "STAR":
    rule star_map_reads:
        input:
            fq1 = expand(OUTcut+"/{sample}_R1.fastq.gz", sample=config["Samples"]),
            fq2 = expand(OUTcut+"/{sample}_R2.fastq.gz", sample=config["Samples"])
        output:
            "star/sam/{sample}/Aligned.out.sam"
        params:
            index= INDEXSTAR,
            extra=""
        threads: THREADS
        message:
            "Quantification with STAR"
        wrapper:
            "0.47.0/bio/star/align"

    rule samtools_faidx:
        input:
            GENOME
        output:
            GENOME+".fai"
        params: ""
        message:
            "Samtools faidx ..."
        wrapper:
            "0.47.0/bio/samtools/faidx"

    rule samtools_view:
        input:
            "star/sam/{sample}/Aligned.out.sam",
            BT = GENOME+".fai"
        output:
            "star/bam/{sample}/Aligned.out.bam"
        threads: THREADS
        message:
            "Samtools view ..."
        params:
            "-bt "+GENOME+".fai -@ 12"
        wrapper:
            "0.47.0/bio/samtools/view"

    rule samtools_sort:
        input:
            "star/bam/{sample}/Aligned.out.bam"
        output:
            "star/bam/{sample}/Aligned.out.sorted.bam"
        params:
            ""
        threads:  THREADS
        message:
            "Samtools sort ..."
        wrapper:
            "0.47.0/bio/samtools/sort"

    rule samtools_index:
        input:
            "star/bam/{sample}/Aligned.out.sorted.bam"
        output:
            "star/bam/{sample}/Aligned.out.sorted.bam.bai"
        threads: THREADS
        message:
            "Samtools index ..."
        params:
            "-@ 12"
        wrapper:
            "0.47.0/bio/samtools/index"

    rule htseqcount:
        input:
            BAM = "star/bam/{sample}/Aligned.out.sorted.bam",
            GTF = GTF,
            BAI ="star/bam/{sample}/Aligned.out.sorted.bam.bai"
        output:
            QUANTIF+"/{sample}/quantif.txt"
        conda:
            "Tools/htseq.yaml"
        message:
            "Running HTseq-count ..."
        shell:
            "htseq-count -f bam "
            "-s reverse -r pos -i gene_name "
            "{input.BAM} "
            "{input.GTF} "
            "> {output}"
