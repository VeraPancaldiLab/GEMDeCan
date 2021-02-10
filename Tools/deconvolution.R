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
})

TPM_path <- snakemake@input[[1]]
signatures_path <- snakemake@params[[1]]
output_path <- snakemake@output[[1]]

signature_files <- list.files(signatures_path, full.names = T)


TPM <- readr::read_tsv(TPM_path)
TPM_matrix <- as.matrix(TPM[, -1])
rownames(TPM_matrix) <- TPM %>% pull(1)

quantiseq <- as_tibble(deconvolute(TPM_matrix, "quantiseq")) %>%
  pivot_longer(-cell_type) %>%
  pivot_wider(names_from = cell_type, values_from = value) %>%
  rename(sample = name)

colnames(quantiseq)[-1] <- str_c("Quantiseq_", colnames(quantiseq)[-1])
colnames(quantiseq) <- sapply(colnames(quantiseq), . %>% {
  str_replace_all(., " ", "_")
})

mcp <- as_tibble(MCPcounter.estimate(TPM_matrix, featuresType = "HUGO_symbols"), rownames = "cell_type") %>%
  pivot_longer(-cell_type) %>%
  pivot_wider(names_from = cell_type, values_from = value) %>%
  rename(sample = name)

colnames(mcp)[-1] <- str_c("MCP_", colnames(mcp)[-1])
colnames(mcp) <- sapply(colnames(mcp), . %>% {
  str_replace_all(., " ", "_")
})

all_deconvolutions_table <- inner_join(quantiseq, mcp, by = "sample")

TPM_df <- as.data.frame(TPM_matrix)
for (signature_file in signature_files) {
  signature <- as.matrix(read.table(signature_file, header = TRUE, row.names = 1, sep = "\t"))

  signature_name <- basename(signature_file)
  signature_name <- str_split(signature_name, "-")[[1]][2]
  signature_name <- str_sub(signature_name, end = -5)

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

  all_deconvolutions_table <- inner_join(all_deconvolutions_table, epi, by = "sample") %>% inner_join(decon, by = "sample")
}



all_deconvolutions_table %>% write_tsv(output_path)
