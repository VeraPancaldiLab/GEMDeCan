# quick pipeline for bam files quantif with salmon (and deconvolution because why not)

configfile: "config.yaml"


sampledir = config["Samples"]
SAMPLES = list(open(sampledir).read().splitlines())

#output
rule all:
    input:
        "output/deconvolution.txt"

rule salmon:
    input:
        r1 = "/{samples}.bam",
        transcript = INDEXS
    output:
        quant = "quantif/{samples}/quant.sf",
        lib = "quantif/{samples}/lib_format_counts.json"
    params:
        DIR = "quantif/{samples}",
        libtype ="A",
        extra=" --validateMappings"
    threads: THREADS
    message:
        "Quantification with Salmon"
    conda:
        "Tools/salmon.yaml"
    shell:
        "salmon quant -t {input.transcript} -l {params.libtype} "
        "-a {input.r1}"
        "-o {params.DIR} "
        "-p {threads} --validateMappings"

rule salmon_quant:
    input:
        "quantif/{samples}/quant.sf"
    output:
        "quantif/{samples}/quantif.txt"
    params:
        "quantif/{samples}"
    conda:
        "Tools/quantif.yaml"
    script:
        "Tools/quant_for_salmon.R"

rule merge_quantif:
    input:
        expand("quantif/{samples}/quantif.txt", samples= SAMPLES)
    output:
        "quantif/all_sample_quantified.txt"
    script:
        "Tools/merge_quantif.R"

rule quantiseq:
    input:
        "quantif/all_sample_quantified.txt"
    output:
        "output/deconvolution.txt"
    message:
        "Running deconvolution"
    conda:
        "Tools/immunedeconv.yaml"
    script:
        "Tools/deconvolution_quantiseq.R"