suppressMessages({
  require(tximport)
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Organism.dplyr)
})
setwd(snakemake@params[[1]])
samples <- snakemake@params[[2]]
files <- sort(paste0(samples, "/quant.sf"))
names(files) <- sort(samples)

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene", overwrite=T)
k <- keys(src, keytype = "tx_id")
tx2gene <- select(src, keys = k, columns = c("tx_name", "symbol"), keytype = "entrez")
tx2gene <- tx2gene[, -1]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# Get TPM
txi_TPM <- as.data.frame(cbind(Gene = rownames(txi$abundance), txi$abundance))
write.table(txi_TPM, "../TPM.txt", sep = "\t", quote = F, row.names = F)
# Output Gene counts for good measure
txi_count <- as.data.frame(cbind(Gene = rownames(txi$counts), txi$counts))
write.table(txi_count, "../gene_counts.txt", sep = "\t", quote = F, row.names = F)
