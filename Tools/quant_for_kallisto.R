suppressMessages({
  library(tximport)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Organism.dplyr)
})
setwd(snakemake@params[[1]])
samples <- snakemake@params[[2]]
files <- sort(paste0(samples, "/abundance.h5"))
names(files) <- sort(samples)

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

tx2gene <- dplyr::inner_join(tbl(src, "id"), tbl(src, "ranges_tx"), by = "entrez") %>%
  dplyr::select(tx_name, symbol) %>%
  as_tibble()

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, txIn = TRUE)
# Get TPM
txi_TPM <- as.data.frame(cbind(Gene = rownames(txi$abundance), txi$abundance))
write.table(txi_TPM, "../TPM.txt", sep = "\t", quote = F, row.names = F)
# Output Gene counts for good measure
txi_count <- as.data.frame(cbind(Gene = rownames(txi$counts), txi$counts))
write.table(txi_count, "../gene_counts.txt", sep = "\t", quote = F, row.names = F)
