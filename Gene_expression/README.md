# Pipeline for RNA-seq and deconvolution analysis

## Introduction

This computational pipeline takes as input BCL or FASTQ files of RNA-seq reads, performs trimming, quantification and deconvolution with the following softwares :
<p align="center">
  <img src="diagram2.png?raw=true" />
</p>

**Trimming** :
* Trim-galore
* Trimmomatic

**RNASeq quantification** :
* Kallisto
* STAR + RSEM
* Salmon

**Deconvolution** (all methods are performed):
* QuantiSeq
* MCP Counter
* deconRNAseq
* EpiDish



## Installation
Snakemake allows for a very efficient and user friendly way of using pipelines. It is designed so all you need to install is Conda (required to install Snakemake and mamba) and Snakemake itself.

Note that officially only Linux is supported for this pipeline. This also requires an Internet connection in order to use conda auto-generated environments for all necessary software and packages.
All R scripts within this pipeline were written with R version 4.1.2. You may experience packages conflicts while running previous version of R.

Conda is required for this pipeline, here's the installation process :
* [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

Mamba
* [Linux](https://github.com/mamba-org/mamba)
`conda install mamba -n base -c conda-forge`

Snakemake
* [Linux](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
`conda activate base`
`mamba install snakemake -n base -c conda-forge -c bioconda`

## Dependencies
All dependencies are automatically installed by snakemake using the .yaml configuration files in Tools. So the user does not need to configure the conda environments as it is done within the pipeline. For information, here are the versions of libraries that are used :
* r-rmarkdown=2.7
* r-hmisc=4.4.1
* r-viridis=0.5.1
* bcl2fastq=2.19.0
* r-immunedeconv=2.0.1
* bioconductor-deconrnaseq=1.28.0
* bioconductor-epidish=2.2.0
* r-mcpcounter=1.1.0
* fastqc=0.11.9
* bioconductor-tximport=1.18.0
* bioconductor-rhdf5=2.34.0
* bioconductor-organism.dplyr=1.18.0
* bioconductor-org.hs.eg.db=3.12.0
* bioconductor-txdb.hsapiens.ucsc.hg38.knowngene=3.10.0
* rename=1.601
* rsem=1.3.3
* salmon=1.1.0
* star=2.7.6a

## Configure your workspace
The snakefile shouldn't be modified. A provided `config.yaml` file takes as input all needed directories and files.

### General informations
 * **Output directory** : directory for all outputs. If you have multiple dataset, make one for each dataset.
 * **Input** : directory where your RNASeq data are located (`.fastq` or `.bcl`). If using only deconvolution, the path to your quantification matrix file (tab-separated TPM values). Note that your .fastq files should end with "R1" and "R2" to be properly processed by the pipeline.
 * **Threads** : number of threads allowed for each job.

### Options
 * **Do_deconv** : should the pipeline run a deconvolution method ? "yes" or "no"
 * **Do_rnaseq** : should the pipeline run the RNAseq analysis part ? "yes" or "no"
 * **Convert_bcl2fastq** : do you need to convert `.bcl` to `.fastq` ? "yes" or "no"
 * **Compute_index** : do you need to perform genome indexing ? "yes" or "no". If "yes", provide the required files to compute the indexes according to your quantification method. If "no", i.e. if you already have a genome index ready to use, simply the path(s) to the index or folder in the configuration file.

### RNASeq
 * **Quantification_with** : STAR, Kallisto or Salmon to be used for quantification analysis
 * **Index_rnaseq** : index file location for Kallisto or Salmon, index folder location for STAR method. Make sure the index you're using was build with the right version of the method selected.
 * **Sample_sheet** results from illumina sequencing. It is needed for Illumina `.bcl` to `.fastq` conversion (Convert_bcl2fastq = "yes").
 * **Adapter** : Path to the adapter used for illumina sequencing that is to be trimmed. Required for Trimmomatic, but not for Trim-galore.
 * **Samples** : the list of all samples to be analysed. It should be a path to a `.txt` file with the list of samples in it. Only required if you don't run `bcl2fastq`.
  * **Trim_with** : chose between one of the two trimmer

 ### Index building
 * **CDNA** : Required for Salmon and Kallisto. Path to the CDNA file. E.g : Homo_sapiens.GRCh38_cdna.fa.gz
 * **Genome** : Path to genome file in fasta format. Note that this is also required to run STAR + RSEM, even if you don't build the index. E.g. : Homo_sapiens.GRCh38.dna.primary_assembly.fa

 ### STAR specific files
 * **GTF** : Path to file (absolute path is mandatory). E.g. : Homo_sapiens.GRCh38.103.gtf

### Deconvolution
 * **Signatures** : Path to folder containing signatures to use for the deconvolution



## Usage
Once everything is configured and installed, open a terminal on the `snakefile` location.
Activate your conda environment and launch the pipeline using a single bash line :
`snakemake -j <number_of_threads> --use-conda`
The `<number_of_threads>` parameter in the command line can be different from the parameter in the config file. For example, if you give 4 threads in the config file and 8 in the snakemake command line, 2 jobs can run in parallel.

Note that Kallisto and Salmon are faster than STAR+RSEM, as they are pseudo-aligner.

### Containers
The pipeline is ready to be executed using **[Singularity](https://sylabs.io/singularity/)** containers\
Simply install the last version of _Singularity_ following the [official documentation](https://sylabs.io/guides/3.6/user-guide/quick_start.html#quick-installation-steps) and add the **`--use-singularity`** flag.\
`snakemake -j <number_of_threads> --use-conda --use-singularity`

### Running in an HPC/Cluster/Supercomputer
To facilitate the process, a *Singularity Definition File* is provided so **[Singularity](https://sylabs.io/singularity/)** has to be installed in the local machine to generate the container and in the HPC to be able to run it.

To generate the Singularity container simply execute the next command in your local machine (sudo/root permissions are required)\
`sudo singularity build --force singularity_container.sif singularity.def`

To run the pipeline in the HPC (without sudo/root permissions):\
`singularity run singularity_container.sif snakemake --cores all`

## Deconvolution
Last part of the pipeline runs a deconvolution algorithm on the quantified samples.
We here chose to run QuantiSeq through the R `immunedeconv` package which wraps several other algorithms.\
DeconRNASeq and EpiDish both allow the user to choose his own signature.

Please do take note that all methods require a quantification matrix as input, in a tabulated format with a preference for **TPM** normalization (we also provide a quick R function of convert from FPKM). The signature is also in a tab separated file. You can find our signature from the publication [link to publication] in the `Signature/` directory of this repository.

## Output files
### RNASeq quantification
This part outputs 2 files :
* **TPM.txt** : matrix of quantification normalized in TPM, this is the one used to run dconvolution.
* **gene_counts.txt** : matrix of quantification without normalization. Not used in the rest of the pipeline but might be useful for the user.
### Deconvolution
* **deconvolution.txt** : results of deconvolution : a table of cell type per sample.
* An HTML summary report of the deconvolution results.

## Example
You can run an example of this pipeline using ressources provided in the `Example/` directory of this repository.
