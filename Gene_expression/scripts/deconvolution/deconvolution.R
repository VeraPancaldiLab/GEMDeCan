library(argparse)
library(readr)
library(tibble)
library(dplyr)
source("scripts/deconvolution/deconvolution_algorithms.R")

parser <- ArgumentParser(description = "Run all deconvolutions algorithms and signatures")
parser$add_argument("--expression_file", type = "character", required = "True", help = "Expression file (required)")
parser$add_argument("--output_file", type = "character", required = "True", help = "Output file (required)")
# args <- parser$parse_args(c("--input_ref", "inputs_paper/CIBERSORTx_ref.txt", "--input_expr", "inputs_paper/GSE107011.txt", "--ref_type", "rna-seq", "--mail", "ting.xie@inserm.fr", "--token","4e9bc3069acda02c6e51980cb3bd5625"))
# "--input_ref inputs_paper/CIBERSORTx_ref.txt --input_expr inputs_paper/GSE107011.txt --ref_type rna-seq --mail ting.xie@inserm.fr --token 4e9bc3069acda02c6e51980cb3bd5625"
args <- parser$parse_args(commandArgs(trailingOnly = T))

rnaseq_tpm <- read_tsv(args$expression_file)

if (length(rnaseq_tpm %>% pull(1)) != length(unique(rnaseq_tpm %>% pull(1)))) {
  rnaseq_tpm_matrix <- rnaseq_tpm %>%
    group_by_at(colnames(.)[1]) %>%
    summarise_all(mean) %>%
    column_to_rownames(colnames(.)[1]) %>%
    as.matrix()
} else {
  rnaseq_tpm_matrix <- rnaseq_tpm %>%
    column_to_rownames(colnames(.)[1]) %>%
    as.matrix()
}

rnaseq_tpm <- NULL

rnaseq_tpm_deconvolutions <- deconvolutions(rnaseq_tpm_matrix)

rnaseq_tpm_deconvolutions %>%
  as_tibble(rownames = "Samples") %>%
  write_tsv(args$output_file)
