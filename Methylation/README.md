# Pipeline for Methylation And Deconvolution analysis

## Introduction

This computational pipeline takes as input `.bcl` files and process them into a beta vlues matrix before running deconvolution with EpiDish.
<p align="center">
  <img src="/diagram.png?raw=true" />
</p>

## Installation 

Snakemake allows for a very efficient and user friendly way of using pipelines. It is designed so all you need to install is Conda (required to install Snakemake) and Snakemake itself.
You can refer to this link for installation, using the "Installation via Conda" chapter : [Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Note that if you don't use Conda, this pipeline won't work. In case you need it, here's the installation process :
* [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

Once the snakemake environment is activated, you need to add the following channels :
* `conda config --add channels bioconda`
* `conda config --add channels conda-forge`

Note that officialy only Linux is supported for this pipeline. This also requires an internet connection in order to use conda auto-generated environements for all necessary softwares and packages.
All R scripts within this pipeline were written with R version 4.0.3. You may experience packages conflicts while running previous version of R.

## Dependencies
All dependencies are automatically installed by snakemake using the .yaml configuration files. The following versions of libraries are used : 
* bioconductor-epidish= 2.2.0
* bioconductor-champ =2.16.0
* bioconductor-minfi =1.32.0
* r-zoo =1.8_8

## Configure your workspace
The snakefile shouldn't be modified. A provided `config.yaml` file takes as input all needed directories and files.

### General informations
 * **Output** : directory for all outputs. If you have multiple dataset, make one for each dataset.
 * **Data** : directory where your data are located (`.idat`). If using only deconvolution, the path to your quantification matrix file (tab-separated beta values).
 * **Threads** : number of threads allowed for each job.
 * **Array_type** : the type of Illumina methylation array used for sequencing : 450K or EPIC
 * **Phenotype_file** : `.csv` file with annotations for your data. An exemple is provided to make sure you use the correct formatting. 
 * **Annotation_file** : Illumina methylation array annotation file (.csv format) which can be found [here on Illumina website.](https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html)
 * **Do_beta_values** : wether you need to compute beta values from raw `.idat`data or not
 * **Do_deconv** : f you want to perform deconvolution using EpiDish
 * **Signature** : additional signatures to perform the deconvolution with. Please note that as this is for methylation analysis, your signatures need to feature CPGs.
 
 ## Usage
Once everything is configured and installed, open a terminal on the `snakefile` location.
Activate your conda environement and launch the pipeline using a single bash line :
`snakemake -j <number_of_threads> --use-conda`
The `<number_of_threads>` parameter in the command line can be different from the parameter in the config file. For exemple, if you give 4 threads in the config file and 8 in the snakemake command line, 2 jobs can run in parallel. 

## Deconvolution
Last part of the pipeline runs the EpiDish deconvolution algorithm on the quantified samples. 
EpiDish is ran with its signature (see publication) and also run other signatures you might provide in a separate folder. We provide our BPMetCan signature in this folder.

Please do take note that all methods require a quantification matrix as input, in a tabulated format with beta values. The signature is also in a tab separated file, check the `/data/signature`folder for exemple. 

## Output files
This pipeline will generate 2 output files : 
* `deconvolution.txt` : results of the deconvolution step
* `results/beta_values.txt` : matrix of beta values from the methylation data
