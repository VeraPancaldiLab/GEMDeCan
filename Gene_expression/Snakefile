import yaml
import os
configfile: 'snakemake.yml'

CBSx_SIGNATURES = config["CBSx_SIGNATURES"]

MIXTURES = config["MIXTURES"]

INPUT_FOLDER = config["INPUT_FOLDER"]
OUTPUT_FOLDER = config["OUTPUT_FOLDER"]
CIBERSORTx_OUTPUT_FOLDER = OUTPUT_FOLDER  + "/cibersortx"

credentials_yml = None
with open(config["CREDENTIALS_FILE"]) as f:
    credentials_yml = yaml.safe_load(f)

MAIL = credentials_yml["MAIL"]
TOKEN = credentials_yml["TOKEN"]

rule cibersortx:
  input:
    signature = INPUT_FOLDER + "/{signature}.txt",
    mixture = INPUT_FOLDER + "/{mixture}.txt"
  output:
    CIBERSORTx_OUTPUT_FOLDER + "/{signature}___CIBERSORTx_fractions___{mixture}/CIBERSORTx_Results.txt"
  shell:
    "Rscript scripts/run_cibersort_local_container.R --signature_file {input.signature} --expression_file {input.mixture} --mail " + MAIL + " --token " + TOKEN + " --output_folder " + CIBERSORTx_OUTPUT_FOLDER

rule curate_cibersortx_results:
  input:
    CIBERSORTx_OUTPUT_FOLDER + "/{signature}___CIBERSORTx_fractions___{mixture}/CIBERSORTx_Results.txt",
  output:
    CIBERSORTx_OUTPUT_FOLDER + "/{signature}___CIBERSORTx_fractions___{mixture}/CIBERSORTx_Results_curated.txt",
  shell:
    "Rscript scripts/curate_cibersortx_results.R --cibersortx_results_file {input} --signature {wildcards.signature} --output_file {output}"


rule deconvolutions:
  input:
    mixture = INPUT_FOLDER + "/{mixture}.txt"
  output:
    OUTPUT_FOLDER + "/deconvolutions/deconvolutions_{mixture}.txt"
  shell:
    "Rscript scripts/deconvolution/deconvolution.R --expression_file {input.mixture} --output_file {output}"


rule all_deconvolutions_one_table:
  input:
    OUTPUT_FOLDER + "/deconvolutions/deconvolutions_{mixture}.txt",
    expand(CIBERSORTx_OUTPUT_FOLDER + "/{signature}___CIBERSORTx_fractions___{{mixture}}/CIBERSORTx_Results_curated.txt", signature=CBSx_SIGNATURES)
  output:
    OUTPUT_FOLDER + "/all_deconvolutions/all_deconvolutions_{mixture}.txt"
  shell:
    "Rscript scripts/merge_all_deconvolutions.R --deconvolution_files {input} --output_file {output}"

real_user = os.environ["SUDO_USER"]
rule fix_root_permissions:
  input:
    expand(OUTPUT_FOLDER + "/all_deconvolutions/all_deconvolutions_{mixture}.txt", mixture = MIXTURES),
  output:
    OUTPUT_FOLDER + "/.fixed_permissions"
  shell:
    "touch {output} && chown -R " + real_user + ":" + real_user + " " + OUTPUT_FOLDER

rule all:
  input:
    expand(expand(CIBERSORTx_OUTPUT_FOLDER + "/{signature}___CIBERSORTx_fractions___{{mixture}}/CIBERSORTx_Results_curated.txt", signature=CBSx_SIGNATURES), mixture = MIXTURES),
    OUTPUT_FOLDER + "/.fixed_permissions",
    expand(OUTPUT_FOLDER + "/deconvolutions/deconvolutions_{mixture}.txt", mixture = MIXTURES),
    expand(OUTPUT_FOLDER + "/all_deconvolutions/all_deconvolutions_{mixture}.txt", mixture = MIXTURES),


