library(readr)
library(argparse)
library(dplyr)
library(stringr)
library(magrittr)


parser <- ArgumentParser(description = "Curate CIBERSORTx_Results file")
parser$add_argument("--cibersortx_results_file", type = "character", required = "True", help = "cibersortx_results.txt file")
parser$add_argument("--signature", type = "character", required = "True", help = "Reference")
parser$add_argument("--output_file", type = "character", required = "True", help = "Output file")

args <- parser$parse_args(commandArgs(trailingOnly = T))

filename_sig <- args$signature
# print("filename_sig")
# print(filename_sig)

results <- read_tsv(args$cibersortx_results_file)
colnames(results) <- sapply(colnames(results), . %>%
  {
    str_replace_all(str_c("CBSX__", filename_sig, "__", .), " ", "_")
  })

colnames(results)[1] <- "Samples"
results %<>% select(!contains("RMSE") & !contains("P-value") & !contains("Correlation"))
results %>% write_tsv(args$output_file)
