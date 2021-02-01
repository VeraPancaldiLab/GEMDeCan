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

from os.path import basename
from os.path import abspath
from re import sub
import csv

configfile: "config.yaml"

##########################
####     PARAMETERS     ####
#########################

# Define paths
OUTDIR = config["Output_Directory"]
INDIR = config["Input_Directory"]
SAMPLESHEET = config["Sample_Sheet"]
THREADS = config["THREADS"]
INDEX = config["Index_rnaseq"]
GENOME = config["Genome"]
GTF = config["GTF"]
ADAPTER = config["Adapter"]
sampledir = config["Samples"]
SIGNATURE = config["Signature"]
QUANTIFTOOL = config["Deconvolution_method"]

OUTbcl = OUTDIR + "/bcl2raw"
OUTfastqc = OUTDIR + "/fastqc_raw"
OUTfastqc2 = OUTDIR + "/fastqc_cutadapter"
OUTmultiqc = OUTDIR + "/multiqc_raw"
OUTmultiqc2 = OUTDIR + "/multiqc_after_cutadapter"
OUTcut = OUTDIR + "/data_after_cutadapter"
QUANTIF = OUTDIR + "/Quantification"

deconv_with_sign = ["epidish", "deconRNAseq"]
deconv_without_sign = ["mcpcounter", "quantiseq"]

##########################
#### Exceptions handling ####
#########################

def exit_error(message):
    exit("ERROR: Exiting Snakemake procedure due to missing \"{}\" parameter in the config.yaml file.".format(message))

if config["Do_deconv"] == "yes":

    if  config["Deconvolution_method"] is None:
        exit_error("Deconvolution_method")

    if config["Signature"] is None and config["Deconvolution_method"] in deconv_with_sign:
        exit_error("Signature")

    if config["Deconvolution_method"] in deconv_with_sign:
        SIG_name = basename(SIGNATURE)
        SIG_name = sub(".txt", "", SIG_name, 1)

if config["Do_rnaseq"] == "yes":

    if config["Trim_with"] is None:
        exit_error("Trim_with")

    if config["Quantification_with"] is None:
        exit_error("Quantification_with")

    if config["Quantification_with"] == "STAR":
        if config["GTF"] is None:
            exit_error("GTF")
        if config["Genome"] is None:
            exit_error("Genome")

    if config["Index_rnaseq"] is None:
        exit_error("Index_rnaseq")

    if config["Convert_bcl2fastq"] is None:
        exit_error("Convert_bcl2fastq")

    if config["Convert_bcl2fastq"] == "yes" and config["Sample_Sheet"] is None:
        exit_error("Sample_Sheet")

    if config["Convert_bcl2fastq"] == "no" and config["Samples"] is None:
        exit_error("Samples")

    if config["Convert_bcl2fastq"] == "yes":
        # extract samples names from Illumina Sample Sheet.
        file_out = basename(SAMPLESHEET)
        file_out = sub(".csv", "", file_out, 1)
        file_out = "Samples_" + file_out + ".txt"
        print(file_out)
        with open(SAMPLESHEET, "r") as f_in, open(file_out, "w") as f_out:
            for skip in range(18):
                next(f_in)
            reader = csv.reader(f_in, delimiter=',')
            SAMPLES = list(set([row[2] for row in reader]))
            f_out.write(SAMPLES[0])
            for item in SAMPLES[1:]:
                    f_out.write('\n{}'.format(item))

    ## else if Do_rnaseq = yes and Convert_bcl2fastq = no, Path to file of Samples' names is required
    else:
        SAMPLES = list(open(sampledir).read().splitlines())


##########################
####       OUTPUTS          ####
#########################
if config["Do_deconv"] == "yes" and config["Do_rnaseq"] == "yes" and config["Deconvolution_method"] in deconv_without_sign:
    rule all:
        input:
            expand(OUTmultiqc + "/{sample}_multiqc_report.html", sample=SAMPLES),
            expand(OUTmultiqc2 + "/{sample}_multiqc_report.html", sample=SAMPLES),
            OUTDIR + "/deconvolution_" + QUANTIFTOOL + ".txt",
            directory(OUTDIR + "/HTML_REPORT_" + QUANTIFTOOL)
elif config["Do_deconv"] == "yes" and config["Do_rnaseq"] == "no" and config["Deconvolution_method"] in deconv_without_sign:
    rule all:
        input:
            OUTDIR + "/deconvolution_" + QUANTIFTOOL + ".txt",
            directory(OUTDIR + "/HTML_REPORT_" + QUANTIFTOOL)

elif config["Do_deconv"] == "yes" and config["Do_rnaseq"] == "yes" and config["Deconvolution_method"] in deconv_with_sign:
    rule all:
        input:
            expand(OUTmultiqc + "/{sample}_multiqc_report.html", sample=SAMPLES),
            expand(OUTmultiqc2 + "/{sample}_multiqc_report.html", sample=SAMPLES),
            OUTDIR + "/deconvolution_" + QUANTIFTOOL + "_" + SIG_name + ".txt",
            directory(OUTDIR + "/HTML_REPORT_" + QUANTIFTOOL + "_" + SIG_name)
elif config["Do_deconv"] == "yes" and config["Do_rnaseq"] == "no" and config["Deconvolution_method"] in deconv_with_sign:
    rule all:
        input:
            OUTDIR + "/deconvolution_" + QUANTIFTOOL + "_" + SIG_name + ".txt",
            directory(OUTDIR + "/HTML_REPORT_" + QUANTIFTOOL + "_" + SIG_name)

elif config["Do_deconv"] == "no" and config["Do_rnaseq"] == "yes":
    rule all:
        input:
            expand(OUTmultiqc + "/{sample}_multiqc_report.html", sample=SAMPLES),
            expand(OUTmultiqc2 + "/{sample}_multiqc_report.html", sample=SAMPLES),
            OUTDIR + "/all_sample_quantified.txt"

########################
####### RNA-SEQ #######
#######################

if config["Do_rnaseq"] == "yes":

    ## Converts base call (.BCL) files into FASTQ
    if config["Convert_bcl2fastq"] == "yes":
        rule bcl2fastq:
            input:
                INDIR = INDIR,
                SHEET = SAMPLESHEET
            output:
                OUTbcl + "/Reports/html/index.html"
            params:
                OUTbcl
            message:
                "Converting BCL2 Illumina files to fastQ"
            conda:
                "Tools/bcl2fastq.yaml"
            threads: THREADS
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            shell:
                """
                bcl2fastq -w {threads} -R {input.INDIR} -o {params} --sample-sheet {input.SHEET} --no-lane-splitting
                rm {params}/Undetermined*
                """

        rule rename_raw:
            input:
                OUTbcl + "/Reports/html/index.html"
            output:
                expand(OUTbcl + "/{samples}_R1.fastq.gz", samples= SAMPLES),
                expand(OUTbcl + "/{samples}_R2.fastq.gz", samples= SAMPLES)
            params:
                OUTbcl
            conda:
                "Tools/rename.yaml"
            shell:
                """
                rename "s/_S[0-9]+_R1_001/_R1/g" {params}/*.fastq.gz
                rename "s/_S[0-9]+_R2_001/_R2/g" {params}/*.fastq.gz
                """

    # declare QCINPUT for next processing
    if  config["Convert_bcl2fastq"] == "yes":
        QCINPUT = OUTbcl
    else:
        QCINPUT = INDIR

    ## Quality control for raw fastq data
    rule fastqc1:
        input:
            QCINPUT + "/{samples}.fastq.gz"
        output:
            html = OUTfastqc + "/{samples}_fastqc.html",
            zip = OUTfastqc + "/{samples}_fastqc.zip"
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
            OUTfastqc + "/{samples}_R1_fastqc.html",
            OUTfastqc + "/{samples}_R2_fastqc.html"
        benchmark:
            "benchmarks/benchmark.multiqc1_{samples}.txt"
        output:
            MainOut = OUTmultiqc + "/{samples}_multiqc_report.html"
        wrapper:
            "0.47.0/bio/multiqc"

    if config["Trim_with"] == "Trimmomatic":
        ## Read trimming by Trimmomatic (Paired-End)
        rule Trimmomatic:
            input:
                r1 = QCINPUT + "/{samples}_R1.fastq.gz",
                r2 = QCINPUT + "/{samples}_R2.fastq.gz"
            output:
                r1 = OUTcut + "/{samples}_R1.fastq.gz",
                r1_unpaired = OUTcut + "/{samples}_R1.unpaired.fastq.gz",
                r2 = OUTcut + "/{samples}_R2.fastq.gz",
                r2_unpaired = OUTcut + "/{samples}_R2.unpaired.fastq.gz"
            threads: THREADS
            message:
                "Trimming using Trimmomatic"
            benchmark:
                "benchmarks/benchmark.trimmomatic_{samples}.txt"
            params:
                trimmer = ["TRAILING:20", "LEADING:20", "MINLEN:36", "CROP:10000", "ILLUMINACLIP:" + ADAPTER + ":2:30:10"],
                extra = "-phred33"
            wrapper:
                "0.47.0/bio/trimmomatic/pe"

    elif config["Trim_with"] == "Trimgalore":
        ## Read trimming by Trim-galore (Paired-end)
        rule trimgalore:
            input:
                QCINPUT + "/{samples}_R1.fastq.gz",
                QCINPUT + "/{samples}_R2.fastq.gz"
            output:
                OUTcut + "/{samples}_R1_val_1.fq.gz",
                OUTcut + "/{samples}_R1.fastq.gz_trimming_report.txt",
                OUTcut + "/{samples}_R2_val_2.fq.gz",
                OUTcut + "/{samples}_R2.fastq.gz_trimming_report.txt"
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
                R1 = OUTcut + "/{samples}_R1_val_1.fq.gz",
                R2 = OUTcut + "/{samples}_R2_val_2.fq.gz"
            output:
                R1out = OUTcut + "/{samples}_R1.fastq.gz",
                R2out = OUTcut + "/{samples}_R2.fastq.gz"
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
            OUTcut + "/{samples}.fastq.gz"
        output:
            html = OUTfastqc2 + "/{samples}_fastqc.html",
            zip = OUTfastqc2 + "/{samples}_fastqc.zip"
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
            OUTfastqc2 + "/{samples}_R1_fastqc.html",
            OUTfastqc2 + "/{samples}_R2_fastqc.html"
        output:
            OUTmultiqc2 + "/{samples}_multiqc_report.html"
        benchmark:
            "benchmarks/benchmark.multiqc2_{samples}.txt"
        wrapper:
            "0.47.0/bio/multiqc"

    # Quantification
    if config["Quantification_with"] == "kallisto":
        rule kallisto:
            input:
                R1 = OUTcut + "/{samples}_R1.fastq.gz",
                R2 = OUTcut + "/{samples}_R2.fastq.gz",
                INDEX = INDEX,
            threads: THREADS
            output:
                QUANTIF + "/{samples}/abundance.tsv",
                QUANTIF + "/{samples}/abundance.h5"
            params:
                OUTDIRE = QUANTIF + "/{samples}"
            message:
                "Quantification with Kallisto"
            benchmark:
                "benchmarks/benchmark.kallisto_{samples}.txt"
            conda:
                "Tools/kallisto.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            shell:
                "kallisto quant -t {threads} -i {input.INDEX} -b 30 "
                "-o {params.OUTDIRE} "
                "{input.R1} {input.R2}"

        rule kallisto_quant:
            input:
                expand(QUANTIF + "/{samples}/abundance.h5", samples= SAMPLES)
            output:
                OUTDIR + "/all_sample_quantified.txt"
            params:
                QUANTIF,
                SAMPLES
            benchmark:
                "benchmarks/benchmark.quant_to_gene.txt"
            conda:
                "Tools/quantif.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/quant_for_kallisto.R"

    elif config["Quantification_with"] == "salmon":
        rule salmon:
            input:
                r1 = OUTcut + "/{samples}_R1.fastq.gz",
                r2 = OUTcut + "/{samples}_R2.fastq.gz",
                index = INDEX
            output:
                quant = QUANTIF + "/{samples}/quant.sf",
                lib = QUANTIF + "/{samples}/lib_format_counts.json"
            params:
                DIR = QUANTIF + "/{samples}",
                libtype ="A",
                extra=" --validateMappings"
            threads: THREADS
            message:
                "Quantification with Salmon"
            benchmark:
                "benchmarks/benchmark.salmon_{samples}.txt"
            conda:
                "Tools/salmon.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            shell:
                "salmon quant -i {input.index} -l {params.libtype} "
                "-1 {input.r1} -2 {input.r2} "
                "-o {params.DIR} "
                "-p {threads} --validateMappings"

        rule salmon_quant:
            input:
                expand(QUANTIF + "/{samples}/quant.sf", samples= SAMPLES)
            output:
                OUTDIR + "/all_sample_quantified.txt"
            params:
                QUANTIF,
                SAMPLES
            benchmark:
                "benchmarks/benchmark.quant_to_gene_{samples}.txt"
            conda:
                "Tools/quantif.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/quant_for_salmon.R"

    elif config["Quantification_with"] == "STAR":
        rule star:
            input:
                fq1 = OUTcut + "/{samples}_R1.fastq.gz",
                fq2 = OUTcut + "/{samples}_R2.fastq.gz"
            output:
                OUTDIR + "/star/{samples}/Aligned.toTranscriptome.out.bam"
            params:
                OUT = OUTDIR + "/star/{samples}/",
                GENOMEdir = INDEX,
                GTF = GTF
            threads:
                THREADS
            benchmark:
                "benchmarks/star_{samples}.txt"
            conda:
                "Tools/star.yaml"
            shell:
                "STAR "
                "--runThreadN {threads} "
                "--runMode alignReads "
                "--readFilesCommand zcat "
                "--outSAMtype BAM SortedByCoordinate "
                "--quantMode TranscriptomeSAM "
                "--quantTranscriptomeBan IndelSoftclipSingleend "
                "--outFileNamePrefix {params.OUT} "
                "--genomeDir {params.GENOMEdir} "
                "--sjdbGTFfile {params.GTF} "
                "--readFilesIn {input.fq1} {input.fq2}"

        rule RSEM_ref:
            input:
                GTF
            output:
                "data/rsem/gen.seq"
            params:
                GEN = GENOME
            threads:
                THREADS
            conda:
                "Tools/rsem.yaml"
            shell:
                "rsem-prepare-reference "
                "-p {threads} "
                "--gtf {input} "
                "{params} "
                "{output}"

        rule RSEM:
            input:
                BAM = OUTDIR + "/star/{samples}/Aligned.toTranscriptome.out.bam",
                REF = "data/rsem/gen.seq"
            output:
                OUTDIR + "/rsem/{samples}.genes.results"
            params:
                OUT = OUTDIR + "/rsem/{samples}",
                REF = "data/rsem/gen"
            threads:
                THREADS
            benchmark:
                "benchmarks/rsem_{samples}.txt"
            conda:
                "Tools/rsem.yaml"
            shell:
                "rsem-calculate-expression "
                "-p {threads} "
                "--paired-end --alignments  "
                "--estimate-rspd "
                "--no-bam-output "
                "--strandedness reverse "
                "{input.BAM} "
                "{params.REF} "
                "{params.OUT}"

        rule star_quant:
            input:
                expand(OUTDIR + "/rsem/{samples}.genes.results", samples= SAMPLES)
            output:
                OUTDIR + "/all_sample_quantified.txt"
            params:
                OUTDIR + "/rsem",
                SAMPLES
            conda:
                "Tools/quantif.yaml"
            script:
                "Tools/quant_for_star.R"

##########################
#### DECONVOLUTION ####
#########################

if config["Do_deconv"] == "yes":
    if config["Do_rnaseq"] == "yes":
        DECONV_INPUT = OUTDIR + "/all_sample_quantified.txt"
    else:
        DECONV_INPUT = INDIR

    if config["Deconvolution_method"] == "quantiseq":
        rule quantiseq:
            input:
                DECONV_INPUT
            output:
                OUTDIR + "/deconvolution_" + QUANTIFTOOL + ".txt"
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.quantiseq.txt"
            conda:
                "Tools/immunedeconv.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_quantiseq.R"

    elif config["Deconvolution_method"] == "mcpcounter":
        rule mcpcounter:
            input:
                DECONV_INPUT
            output:
                OUTDIR + "/deconvolution_" + QUANTIFTOOL + ".txt"
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.mcp.txt"
            conda:
                "Tools/mcpcounter.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_mcpcounter.R"

    elif config["Deconvolution_method"] == "deconRNAseq":
        rule deconRNAseq:
            input:
                DECONV_INPUT
            output:
                OUTDIR + "/deconvolution_" + QUANTIFTOOL + "_" + SIG_name + ".txt"
            params:
                SIGNATURE
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.deconRNA.txt"
            conda:
                "Tools/RNAdeconv.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_deconrnaseq.R"

    elif config["Deconvolution_method"] == "epidish":
        rule epidish:
            input:
                DECONV_INPUT
            output:
                OUTDIR + "/deconvolution_" + QUANTIFTOOL + "_" + SIG_name + ".txt"
            params:
                SIGNATURE
            message:
                "Running deconvolution"
            benchmark:
                "benchmarks/benchmark.epidish.txt"
            conda:
                "Tools/epidish.yaml"
            singularity:
                "docker://continuumio/miniconda3:4.8.2"
            script:
                "Tools/deconvolution_epidish.R"

    if config["Deconvolution_method"] in deconv_without_sign:
        report_input = OUTDIR + "/deconvolution_" + QUANTIFTOOL + ".txt"
        report_output = directory(OUTDIR + "/HTML_REPORT_" + QUANTIFTOOL)
        path_to_deconv = abspath(OUTDIR + "/deconvolution_" + QUANTIFTOOL + ".txt")
    elif config["Deconvolution_method"] in deconv_with_sign:
       report_input = OUTDIR + "/deconvolution_" + QUANTIFTOOL + "_" + SIG_name + ".txt"
       report_output = directory(OUTDIR + "/HTML_REPORT_" + QUANTIFTOOL + "_" + SIG_name)
       path_to_deconv = abspath(OUTDIR + "/deconvolution_" + QUANTIFTOOL + "_" + SIG_name + ".txt")
    else:
        exit("Please check the config.yaml 'Deconcolution_method' parameter.")

    rule report:
        input:
            report_input
        output:
            report_output
        params:
            path_to_deconv
        message:
            "Generating report"
        conda:
            "Tools/analyse.yaml"
        shell:
            "Rscript -e \"rmarkdown::render('Tools/analyses.R', output_dir='{output}')\" {params}"
