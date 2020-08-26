README

# Pipeline for lungProject RNA-seq analysis

## What am I looking at ?
This pipeline goes from converting raw Illumina data to quantification and deconvolution with user-chosen software between :
* Kallisto
* STAR + HTseq-count
* Salmon
* QuantiSeq
* MCP Counter
* deconRNAseq
* EpiDish

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

Once the snakemake environment is activated, you need to add the following channels :
* `conda config --add channels bioconda`
* `conda config --add channels conda-forge`

Note that officialy only Linux is supported for this pipeline.

You may also need to install bash function `rename` to use the BCL to fastq conversion part :
`sudo apt install rename`

## Configure your workspace
The snakefile shouldn't be modified. A provided `config.yaml` file takes as input all needed directories and files.

### General informations

 * **Input directory** : if you change this, please also change all other paths starting with `data` to the corresponding new path
 * **Threads** : number of threads allowed for the rules 
 * **Sample_sheet** results from illumina sequencing. It is needed for Illumina `.bcl` to `.fastq` conversion
 * **Adapter** : adapter used for illumina sequencing that is to be trimmed. Required for Trimmomatic, but not for Trim-galore.
 * **Samples** is the list of all samples to be analysed. It should be a path to a `.txt` file with a list. A dedicated script in the `Tools` directory, `sample_parser.py` can generate it from Illumina files.
### Options
 * **Trim_with** : chose between one of the two trimmer 
 * **Do_deconv** : should the pipeline run a deconvolution method ?
 * **Do_rnaseq** : should the pipeline run the RNAseq analysis part ?
 * **Convert_bcl2fastq** : do you need to convert `.bcl` to `.fastq` ? "yes" or "no"
 * **Need_merging** : do you need to merge the different sequencing lines for each sample ?
### Quantification
 * **Quantification_with** : STAR, Kallisto or Salmon to be used for quantification analysis
 * **Index_kallisto** : index location for Kallisto
 * **Index_salmon** : index location for Salmon
 
 ### Deconvolution 
  * **Deconvolution_method** : run the deconvolution with QuantiSeq, deconRNAseq or MCPCounter. In the last case, you need to provide signature files. See the [MCPCounter doc](https://github.com/ebecht/MCPcounter) for more information
 * **Genes_signature** : a marker based signature file to use with MCPCounter
 * **Signature** : a regular signature to use with deconRNAseq
 
### STAR spcific files
 * **Genome** file in fasta format
 * **GTF** file in absolute path (mendatory !). Can be compressed or not.
 * **Gene length file** : if using STAR to compute deconvolution, you need to provide a file with the length of every gene in order to convert from gene counts to TPM
 * **Mean fragment length** : like above, for STAR with deconvolution you need to provide the mean length of the sequenced reads (in bp)
 * **Index_STAR**: index location for STAR

## Usage
Once everything is configured and installed, open a terminal on the `snakefile` location and launch the pipeline using a single bash line :
`snakemake -j <number_of_threads> --use-conda`

For the sample parser tool, you simply have to use it like a regular python script, with the path to your illumina files as an argument. You might need to modify it to suit your naming convention.
`python3 sample_parser.py <path_to_files> `

If you really need to use STAR before running the deconvolution, you need to provide aditionnal informations described above. The gene length can be obtained by parsing a GTF annotation. An exemple script to make such a file is provided in `Tools/compute_geneLength.R`.

## Deconvolution
Last part of the pipeline runs a deconvolution algorithm on the quantified samples. 
We here chose to run QuantiSeq through the R `immunedeconv` package which wraps several other algorithms.
MCPCounter allows you to use your own signature files, for usage see the link to the official package Git above.
DeconRNASeq and EpiDish both allow the user to chose his own signature. 

Please do take note that all methods require a quantification matrix as input, in a tabulated format. The signature is also in a tab separated file