suppressMessages({
  require(tximport)
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Organism.dplyr)
  require(tidyverse)
})

# Convert ENST to HUGO symbols, quantif.sf is the output from Salmon
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
src <- src_ucsc("Homo sapiens")
k <- keys(src, keytype = "tx_id")
tx2gene<- select(src, keys = k, columns = c("tx_name","symbol"), keytype = "entrez")
tx2gene<-tx2gene[,-1]
txi <- tximport("quant.sf", type = "salmon", tx2gene = tx2gene)
write_tsv(as.data.frame(cbind(Gene=rownames(txi[[3]]), floor(txi[[3]]) )), "gene_length.txt")

