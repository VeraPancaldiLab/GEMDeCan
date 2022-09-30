library(readr)
library(argparse)
library(dplyr)
library(magrittr)
library(purrr)


parser <- ArgumentParser(description = "Join all deconvolution tables for a dataset in one")
parser$add_argument("--deconvolution_files", type = "character", nargs = "+", required = "True", help = "Deconvolution tables")
parser$add_argument("--output_file", type = "character", required = "True", help = "Output file")
args <- parser$parse_args(commandArgs(trailingOnly = T))
deconvolution_files <- args$deconvolution_files


deconvolution_tibbles <- lapply(deconvolution_files, read_tsv)
deconvolutions_one_table <- deconvolution_tibbles %>% purrr::reduce(inner_join, by = "Samples")
deconvolutions_one_table %>% write_tsv(args$output_file)
