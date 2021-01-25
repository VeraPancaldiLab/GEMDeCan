README

# Pipeline for lungProject RNA-seq analysis

## What am I looking at ?
This pipeline goes from converting raw Illumina data to quantification and deconvolution with user-chosen software between :
* Kallisto *
* STAR + RSEM *
* Salmon *
* QuantiSeq
* MCP Counter
* deconRNAseq
* EpiDish

User also gets to chose his favorite trimmer between:
* Trim-galore
* Trimmomatic

\* Please do take note that genome indexing isn't performed in this pipeline. If you don't have a genome index for the software you chose, it won't perform. An additional script `compute_indexes.sh` is provided to help you into building one if you need it.
This script is for command line usage like so : `bash compute_indexes.sh [method] [input] [output] [number_of_threads]`

## Installation

Snakemake allows for a very efficient and user friendly way of using pipelines. It is designed so all you need to install is Conda (required to install Snakemake) and Snakemake itself.
You can refer to this link for installation, using the "Installation via Conda" chapter : [Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Note that if you don't use Conda, this pipeline won't work. In case you need it, here's the installation process :
* [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

Once the snakemake environment is activated, you need to add the following channels :
* `conda config --add channels bioconda`
* `conda config --add channels conda-forge`

Note that officialy only Linux is supported for this pipeline. This also requires an internet connection in order to use conda auto-generated environements for all necessary softwares and packages.
This pipeline was tested under R version 4. You may experience packages conflicts while running previous version of R.

## Configure your workspace
The snakefile shouldn't be modified. A provided `config.yaml` file takes as input all needed directories and files.

### General informations
 * **Output directory** : directory for all outputs. If you have multiple dataset, make one for each dataset.
 * **Input** : directory where your RNASeq data are located (`.fastq` or `.bcl`). If using only deconvolution, the path to your quantification matrix file (tab-separated TPM values).
 * **Threads** : number of threads allowed for each job.

### Options 
 * **Do_deconv** : should the pipeline run a deconvolution method ?
 * **Do_rnaseq** : should the pipeline run the RNAseq analysis part ?

### RNASeq
 * **Quantification_with** : STAR, Kallisto or Salmon to be used for quantification analysis
 * **Index** : index location for Kallisto or Salmon
 * **Sample_sheet** results from illumina sequencing. It is needed for Illumina `.bcl` to `.fastq` conversion.
 * **Adapter** : adapter used for illumina sequencing that is to be trimmed. Required for Trimmomatic, but not for Trim-galore.
 * **Samples** : the list of all samples to be analysed. It should be a path to a `.txt` file with the list of samples in it. If you are using bcl2fastq, simply copy/paste the column sample_name from your sample sheet.
  * **Convert_bcl2fastq** : do you need to convert `.bcl` to `.fastq` ? "yes" or "no"
  * **Trim_with** : chose between one of the two trimmer 
 
 ### STAR specific files
 * **Genome** file in fasta format
 * **GTF** file in absolute path (mendatory !). Can be compressed or not.
 
### Deconvolution 
  * **Deconvolution_method** : run the deconvolution with QuantiSeq, deconRNAseq or MCPCounter. 
 * **Signature** : a regular signature to use with deconRNAseq or EpiDish
 


## Usage
Once everything is configured and installed, open a terminal on the `snakefile` location.
Activate your conda environement and launch the pipeline using a single bash line :
`snakemake -j <number_of_threads> --use-conda`
The `<number_of_threads>` parameter in the command line can be different from the parameter in the config file. For exemple, if you give 4 threads in the config file and 8 in the snakemake command line, 2 jobs can run in parallel. 

Note that Kallisto and Salmon and faster than STAR, as they are pseudo-aligner.

### Containers
The pipeline is ready to be executed using **[Singularity](https://sylabs.io/singularity/)** containers\
Simply install the last version of _Singularity_ following the [official documentation](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps) and add the **`--use-singularity`** flag.\
`snakemake -j <number_of_threads> --use-conda --use-singularity`

## Deconvolution
Last part of the pipeline runs a deconvolution algorithm on the quantified samples. 
We here chose to run QuantiSeq through the R `immunedeconv` package which wraps several other algorithms.\
DeconRNASeq and EpiDish both allow the user to chose his own signature. 

Please do take note that all methods require a quantification matrix as input, in a tabulated format with a preference for **TPM** normalization (we also provide a quick R function of convert from FPKM). The signature is also in a tab separated file. You can find our signature from the publication [link to publication] in the `Signature/` directory of this repository.

## Exemple
You can run an exemple of this pipeline using ressources provided in the `Exemple/` directory of this repository.
