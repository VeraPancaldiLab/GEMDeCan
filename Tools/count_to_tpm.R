if (!require(tidyverse)) {
  install.packages("tidyverse")
}
library(tidyverse)

# HTseq-count outputs a file with a gene count. This is not suitable for deconvolution for many reasons,
# so we convert it to TPM, like Kallisto and Salmon does

setwd(snakemake@params[[1]])
sample <- basename(getwd())

counts.df <- read_tsv("count_quantif.txt", col_names = F)
names(counts.df) <- c("gene_name", "HTseq")
gene.length <- read_tsv(snakemake@params[[2]])
meanFragmentLength <- snakemake@params[[3]]

counts.df <- inner_join(counts.df, gene.length, by = c("gene_name" = "Gene"))
counts.df[, 3] <- floor(counts.df[, 3])

# in order to convert values, we need the length of each gene and  the mean size of sequenced fragments
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {

  # Ensure valid arguments.
  stopifnot(nrow(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))

  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1, function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx, ]
  effLen <- effLen[idx, ]
  featureLength <- featureLength[idx, ]

  # Process one column at a time.
  counts <- as.data.frame(counts)
  effLen <- as.data.frame(effLen)
  tpm <- do.call(cbind, lapply(1, function(i) {
    rate <- log(counts[, i]) - log(effLen[, i])
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the column names from the original matrix.
  tpm <- cbind(Gene = rownames(counts), tpm)
  return(tpm)
}

featureLength <- counts.df[, 3]
counts <- as.matrix(counts.df[, 2])
rownames(counts) <- counts.df$gene_name

TPM.from.HTSeq <- counts_to_tpm(counts, featureLength, meanFragmentLength)
colnames(TPM.from.HTSeq) <- c("Gene", sample)
TPM.from.HTSeq %>%
  as.data.frame() %>%
  write_tsv("quantif.txt")
