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
from os.path import basename

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
GENES = config["Genes_signature"]
sampledir = config["Samples"]
SIGNATURE = config["Signature"]
GENELENGTH = config["Gene_length_file"]
MEANLENGTH = config["Mean_fragment_length"]
QUANTIFTOOL = config["Deconvolution_method"]

OUTbcl = OUTDIR+"/bcl2raw"
OUTmerge = OUTDIR+"/raw2merge"
OUTfastqc = OUTDIR+"/fastqc_raw"
OUTfastqc2 = OUTDIR+"/fastqc_cutadapter"
OUTmultiqc = OUTDIR+"/multiqc_raw"
OUTmultiqc2 = OUTDIR+"/multiqc_after_cutadapter"
OUTcut = OUTDIR+"/data_after_cutadapter"
QUANTIF = OUTDIR+"/Quantification"

SAMPLES = list(open(sampledir).read().splitlines())

SIG_name = basename(SIGNATURE)

## Outputs
if config["Do_deconv"] == "yes" and config["Do_rnaseq"] == "yes":
    rule all:
        input:
            expand(OUTmultiqc+"/{sample}_multiqc_report.html", sample=SAMPLES),
            expand(OUTmultiqc2+"/{sample}_multiqc_report.html", sample=SAMPLES),
            OUTDIR+"/deconvolution_"+QUANTIFTOOL+"_"+SIG_name
elif config["Do_deconv"] == "yes" and config["Do_rnaseq"] == "no":
    rule all:
        input:
            OUTDIR+"/deconvolution_"+QUANTIFTOOL+"_"+SIG_name
elif config["Do_deconv"] == "no" and config["Do_rnaseq"] == "yes":
    rule all:
        input:
            expand(OUTmultiqc+"/{sample}_multiqc_report.html", sample=SAMPLES),
            expand(OUTmultiqc2+"/{sample}_multiqc_report.html", sample=SAMPLES),
            QUANTIF+"/all_sample_quantified.txt"  

if config["Do_rnaseq"] == "yes" :
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
            container: "docker://continuumio/miniconda3:4.8.2"
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
                "bash Tools/merge_those_fastq.sh {params.OUT} {params.OUT}/bcl2raw"

    if config["Convert_bcl2fastq"] == "no" and config["Need_merging"] == "yes":
        rule merge_fastq:
            input:
                INDIR
            output:
                OUTmerge+"/{samples}_R1.fastq.gz",
                OUTmerge+"/{samples}_R2.fastq.gz"
            params:
                    OUT = OUTDIR
            message:
                "Merging files from lanes"
            shell:
                "bash Tools/merge_those_fastq.sh {params.OUT} {input}"

    if config["Need_merging"] == "yes":
        QCINPUT = OUTmerge+"/{samples}.fastq.gz"
    else:
        QCINPUT = INDIR

    ## Quality control for raw fastq data
    rule fastqc1:
        input:
            QCINPUT+"/{samples}.fastq.gz"
        output:
            html = OUTfastqc+"/{samples}_fastqc.html",
            zip = OUTfastqc+"/{samples}_fastqc.zip"
        threads: THREADS
        benchmark:
            "benchmarks/benchmark.fastqc1_{samples}.txt"
        message:
            "Quality control before trimming"
        params:
            "-t 8"
        wrapper:
            "0.47.0/bio/fastqc"

    rule multiqc1:
        input:
            OUTfastqc+"/{samples}_R1_fastqc.html",
            OUTfastqc+"/{samples}_R2_fastqc.html"
        benchmark:
            "benchmarks/benchmark.multiqc1_{samples}.txt"
        output:
            MainOut = OUTmultiqc+"/{samples}_multiqc_report.html"
        wrapper:
            "0.47.0/bio/multiqc"


    if config["Trim_with"] == "Trimmomatic" :
        ## Read trimming by Trimmomatic (Paired-End)
        rule Trimmomatic:
            input:
                r1 = QCINPUT+"/{samples}_R1.fastq.gz",
                r2 = QCINPUT+"/{samples}_R2.fastq.gz"
            output:
                r1 = OUTcut+"/{samples}_R1.fastq.gz",
                r1_unpaired = OUTcut+"/{samples}_R1.unpaired.fastq.gz",
                r2 = OUTcut+"/{samples}_R2.fastq.gz",
                r2_unpaired = OUTcut+"/{samples}_R2.unpaired.fastq.gz"
            threads: THREADS
            message:
                "Trimming using Trimmomatic"
            benchmark:
                "benchmarks/benchmark.trimmomatic_{samples}.txt"
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
                [OUTmerge+"/{sample}_R1.fastq.gz", OUTmerge+"/{sample}_R2.fastq.gz"]
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
            benchmark:
                "benchmarks/benchmark.trimgalore_{samples}.txt"
            wrapper:
                "0.47.0/bio/trim_galore/pe"

        rule rename:
            input:
                R1 = OUTcut+"/{sample}_R1_val_1.fq.gz",
                R2 = OUTcut+"/{sample}_R2_val_2.fq.gz"
            output:
                R1out = OUTcut+"/{sample}_R1.fastq.gz",
                R2out = OUTcut+"/{sample}_R2.fastq.gz"
            benchmark:
                "benchmarks/benchmark.rename_{samples}.txt"
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
        benchmark:
            "benchmarks/benchmark.fastqc2_{samples}.txt"
        params:
            "-t 8"
        wrapper:
            "0.47.0/bio/fastqc"
    rule multiqc2:
        input:
            OUTfastqc2+"/{samples}_R1_fastqc.html",
            OUTfastqc2+"/{samples}_R2_fastqc.html"
        output:
            OUTmultiqc2+"/{samples}_multiqc_report.html"
        benchmark:
            "benchmarks/benchmark.multiqc2_{samples}.txt"
        wrapper:
            "0.47.0/bio/multiqc"

    # Quantification
    if config["Quantification_with"] == "kallisto" :
        rule kallisto_quant:
            input:
                R1 = OUTcut+"/{samples}_R1.fastq.gz",
                R2 = OUTcut+"/{samples}_R2.fastq.gz",
                INDEXK = INDEXK,
            threads: THREADS
            output:
                QUANTIF+"/{samples}/abundance.tsv",
                QUANTIF+"/{samples}/abundance.h5"
            params:
                OUTDIRE = QUANTIF+"/{samples}"
            message:
                "Quantification with Kallisto"
            benchmark:
                "benchmarks/benchmark.kallisto_{samples}.txt"
            conda:
                "Tools/kallisto.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            shell:
                "kallisto quant -t {threads} -i {input.INDEXK} -b 30 "
                "-o {params.OUTDIRE} "
                "{input.R1} {input.R2}"

        rule quant_to_gene:
            input:
                QUANTIF+"/{samples}/abundance.h5"
            output:
                QUANTIF+"/{samples}/quantif.txt"
            params:
                QUANTIF+"/{samples}"
            benchmark:
                "benchmarks/benchmark.quant_to_gene_{samples}.txt"
            conda:
                "Tools/quantif.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/quant_for_kallisto.R"
        
    elif config["Quantification_with"] == "salmon" :
        rule salmon:
            input:
                r1 = OUTcut+"/{samples}_R1.fastq.gz",
                r2 = OUTcut+"/{samples}_R2.fastq.gz",
                index = INDEXS
            output:
                quant = QUANTIF+"/{samples}/quant.sf",
                lib = QUANTIF+"/{samples}/lib_format_counts.json"
            params:
                DIR = QUANTIF+"/{samples}",
                libtype ="A",
                extra=" --validateMappings"
            threads: THREADS
            message:
                "Quantification with Salmon"
            benchmark:
                "benchmarks/benchmark.salmon_{samples}.txt"
            conda:
                "Tools/salmon.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            shell:
                "salmon quant -i {input.index} -l {params.libtype} "
                "-1 {input.r1} -2 {input.r2} "
                "-o {params.DIR} "
                "-p {threads} --validateMappings"
        
        rule salmon_quant:
            input:
                QUANTIF+"/{samples}/quant.sf"
            output:
                QUANTIF+"/{samples}/quantif.txt"
            params:
                QUANTIF+"/{samples}"
            benchmark:
                "benchmarks/benchmark.quant_to_gene_{samples}.txt"
            conda:
                "Tools/quantif.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/quant_for_salmon.R"

    elif config["Quantification_with"] == "STAR":
        rule star_map_reads:
            input:
                fq1 = OUTcut+"/{samples}_R1.fastq.gz",
                fq2 = OUTcut+"/{samples}_R2.fastq.gz"
            output:
                "star/sam/{samples}/Aligned.out.sam"
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
                "star/sam/{samples}/Aligned.out.sam",
                BT = GENOME+".fai"
            output:
                "star/bam/{samples}/Aligned.out.bam"
            threads: THREADS
            message:
                "Samtools view ..."
            params:
                "-bt "+GENOME+".fai -@ 12"
            wrapper:
                "0.47.0/bio/samtools/view"

        rule samtools_sort:
            input:
                "star/bam/{samples}/Aligned.out.bam"
            output:
                "star/bam/{samples}/Aligned.out.sorted.bam"
            params:
                ""
            threads:  THREADS
            message:
                "Samtools sort ..."
            wrapper:
                "0.47.0/bio/samtools/sort"

        rule samtools_index:
            input:
                "star/bam/{samples}/Aligned.out.sorted.bam"
            output:
                "star/bam/{samples}/Aligned.out.sorted.bam.bai"
            threads: THREADS
            message:
                "Samtools index ..."
            params:
                "-@ 12"
            wrapper:
                "0.47.0/bio/samtools/index"

        rule htseqcount:
            input:
                BAM = "star/bam/{samples}/Aligned.out.sorted.bam",
                GTF = GTF,
                BAI ="star/bam/{samples}/Aligned.out.sorted.bam.bai"
            output:
                QUANTIF+"/{samples}/count_quantif.txt"
            conda:
                "Tools/htseq.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            message:
                "Running HTseq-count ..."
            shell:
                "htseq-count -f bam "
                "-s reverse -r pos "
                "-i gene_name "
                "{input.BAM} "
                "{input.GTF} "
                "> {output}"
        
        rule count_to_tpm:
            input:
                QUANTIF+"/{samples}/count_quantif.txt"
            output:
                QUANTIF+"/{samples}/quantif.txt"
            params:
                QUANTIF+"/{samples}",
                GENELENGTH,
                MEANLENGTH
            script:
                "Tools/count_to_tpm.R"

    rule merge_quantif:
        input:
            expand(QUANTIF+"/{samples}/quantif.txt", samples= SAMPLES)
        output:
            QUANTIF+"/all_sample_quantified.txt"
        benchmark:
                "benchmarks/benchmark.merge.txt"
        script:
            "Tools/merge_quantif.R"

if config["Do_deconv"] == "yes":
    if config["Do_rnaseq"] == "yes" :
        DECONV_INPUT = QUANTIF+"/all_sample_quantified.txt"
    else:
        DECONV_INPUT = INDIR
    
    if config["Deconvolution_method"] == "quantiseq":
        rule quantiseq:
            input:
                DECONV_INPUT
            output:
                OUTDIR+"/deconvolution_"+QUANTIFTOOL+"_"+SIG_name
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.quantiseq.txt"
            conda:
                "Tools/immunedeconv.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_quantiseq.R"

    elif config["Deconvolution_method"] == "mcpcounter":
        rule mcpcounter:
            input:
                DECONV_INPUT
            output:
                OUTDIR+"/deconvolution_"+QUANTIFTOOL+"_"+SIG_name
            params:
                GENES
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.mcp.txt"
            script:
                "Tools/deconvolution_mcpcounter.R"

    elif config["Deconvolution_method"] == "deconRNAseq":
        rule deconRNAseq:
            input:
                DECONV_INPUT
            output:
                OUTDIR+"/deconvolution_"+QUANTIFTOOL+"_"+SIG_name
            params:
                SIGNATURE
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.deconRNA.txt"
            conda:
                "Tools/RNAdeconv.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_deconrnaseq.R"

    elif config["Deconvolution_method"] == "epidish":
        rule epidish:
            input:
                DECONV_INPUT
            output:
                OUTDIR+"/deconvolution_"+QUANTIFTOOL+"_"+SIG_name
            params:
                SIGNATURE
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.epidish.txt"
            conda:
                "Tools/epidish.yaml"
            container: "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_epidish.R"