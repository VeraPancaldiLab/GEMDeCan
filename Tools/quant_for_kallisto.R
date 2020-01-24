suppressMessages({
  require(tximport)
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Organism.dplyr)
})
setwd(snakemake@params[[1]])
sample = basename(getwd())

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
src <- src_ucsc("Homo sapiens")
k <- keys(src, keytype = "tx_id")
tx2gene<- select(src, keys = k, columns = c("tx_name","symbol"), keytype = "entrez")
tx2gene<-tx2gene[,-1]
txi <- tximport("abundance.h5", type = "kallisto", tx2gene = tx2gene, txIn=TRUE)
txi_data<-as.data.frame(cbind(Gene=rownames(txi[[1]]), txi[[1]]))
colnames(txi_data) = c("Gene",sample)
write.table(txi_data, "quantif.txt", sep="\t", quote=F, row.names = F)
