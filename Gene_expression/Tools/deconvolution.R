#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages({
  library(dplyr) # rename
  library(tidyr) # pivot_longer
  library(readr) # read_tsv read_csv write_tsv
  library(stringr) # str_c
  library(magrittr) # %<>%
  library(EpiDISH) # epidish
  library(MCPcounter) # MCPcounter.estimate
  library(immunedeconv) # deconvolute
  library(DeconRNASeq) # DeconRNASeq
  library(purrr) # reduce
  library(parallel) # mclapply
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  print("./deconvolution.R TPM_file output_file signatures_path threads_number")
  quit(status = 1)
}

TPM_path <- args[1]
output_path <- args[2]
signatures_path <- args[3]
threads <- args[4]

signature_files <- list.files(signatures_path, full.names = T)
signature_files <- setdiff(signature_files, signature_files[dir.exists(signature_files)])

TPM <- read_tsv(TPM_path)
TPM_matrix <- as.matrix(TPM[, -1])
rownames(TPM_matrix) <- TPM %>% pull(1)

computeQuantiseq <- function(TPM_matrix) {
  quantiseq <- as_tibble(deconvolute(TPM_matrix, "quantiseq")) %>%
    pivot_longer(-cell_type) %>%
    pivot_wider(names_from = cell_type, values_from = value) %>%
    rename(sample = name)

  colnames(quantiseq)[-1] <- str_c("Quantiseq_", colnames(quantiseq)[-1])
  colnames(quantiseq) <- sapply(colnames(quantiseq), . %>% {
    str_replace_all(., " ", "_")
  })
  quantiseq
}


computeMCP <- function(TPM_matrix, signatures_path) {
  genes <- read.table(file.path(signatures_path, "MCPcounter", "MCPcounter-genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  mcp <- as_tibble(MCPcounter.estimate(TPM_matrix, genes = genes, featuresType = "HUGO_symbols", probesets = NULL), rownames = "cell_type") %>%
    pivot_longer(-cell_type) %>%
    pivot_wider(names_from = cell_type, values_from = value) %>%
    rename(sample = name)

  colnames(mcp)[-1] <- str_c("MCP_", colnames(mcp)[-1])
  colnames(mcp) <- sapply(colnames(mcp), . %>% {
    str_replace_all(., " ", "_")
  })
  mcp
}


methods_with_variable_signatures <- function(TPM_matrix, signature_files, threads) {
  TPM_df <- as.data.frame(TPM_matrix)
  all_methods_and_signatures <- mclapply(signature_files, function(signature_file) {
    signature <- as.matrix(read.table(signature_file, header = TRUE, row.names = 1, sep = "\t"))

    signature_name <- basename(signature_file)
    signature_name <- str_split(signature_name, "\\.")[[1]][1]

    epi <- epidish(TPM_matrix, signature, method = "RPC", maxit = 200)
    epi <- as_tibble(epi$estF, rownames = "sample")
    colnames(epi)[-1] <- str_c("Epidish", signature_name, colnames(epi)[-1], sep = "_")
    colnames(epi) <- sapply(colnames(epi), . %>% {
      str_replace_all(., " ", "_")
    })

    decon <- DeconRNASeq(TPM_df, as.data.frame(signature))
    decon <- bind_cols(colnames(TPM_df), as_tibble(decon$out.all)) %>% rename(sample = ...1)
    colnames(decon)[-1] <- str_c("DeconRNASeq", signature_name, colnames(decon)[-1], sep = "_")
    colnames(decon) <- sapply(colnames(decon), . %>% {
      str_replace_all(., " ", "_")
    })

    inner_join(epi, decon, by = "sample")
  }, mc.cores = threads)
  reduce(all_methods_and_signatures, inner_join, "sample")
}

all_deconvolutions_table <- mclapply(c("Quantiseq", "MCP", "rest"), function(method) {
  if (method == "Quantiseq") {
    computeQuantiseq(TPM_matrix)
  } else if (method == "MCP") {
    computeMCP(TPM_matrix, signatures_path)
  } else if (method == "rest") {
    methods_with_variable_signatures(TPM_matrix, signature_files, threads)
  }
}, mc.cores = threads)


all_deconvolutions_table %<>% reduce(inner_join, "sample")

all_deconvolutions_table %>% write_tsv(output_path)
