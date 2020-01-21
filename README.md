# Pipline for lungProject RNA-seq analysis

## What am I looking at ?
This pipeline goes from converting raw Illumina data to quantification with user-chosen software between :
* Kallisto
* STAR + HTseq-count
* Salmon

User also gets to chose his favorite trimmer between:
* Trim-galore
* Trimmomatic

Please do take note that genome indexing isn't performed in this pipeline. If you don't have a genome index for the software you chose, it won't perform. An additional script `compute_indexes.sh` is provided if you need to build it.

## Installation

Snakemake allows for a very efficient and user friendly way of using pipelines. It is designed so all you need to install is Conda (required to install Snakemake) and Snakemake itself.
You can refer to this link for installation, using the "Installation via Conda" chapter : [Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Note that if you don't use Conda, this pipeline won't work. In case you need it, here's the installation process :
* [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
* [Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)
* [MacOS](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)

Note that officialy only Linux is supported for this pipeline.

You may also need to install bash function `rename` to use the BCL to fastq conversion part :
`sudo apt install rename`
## Configure your workspace
The snakefile shouldn't be modified. A provided `config.yaml` file takes as input all needed directories and files.
 * **Samples** is the list of all samples to be analysed
 * **Input directory** : if you change this, please also change all other paths starting with `data` to the corresponding new path
 * **Sample_sheet** results from illumina sequecing. It is needed for Illumina `.bcl` to `.fastq` conversion
 * **Adapter** : adapter used for illumina sequecing that is to be trimmed. Required for Trimmomatic, but not for Trim-galore.
 * **Genome** file in fasta format
 * **GTF** file in absolute path (mendatory !). Can be compressed or not.
 * **Threads** : number of threads allowed for the rules. 
 * **Trim_with** : chose between one of the two trimmer
 * **Quantification_with** : STAR, Kallisto or Salmon to be used for quantification analysis.
 * **Index** : index location for named software
 * **Convert_bcl2fastq** : do you need to convert `.bcl` to `.fastq` ? "yes" or "no"

## Usage
Once everything is configured and installed, open a terminal on the `snakefile` location and launch the pipeline using a single bash line :
`snakemake -j <number_of_threads> --use-conda`


