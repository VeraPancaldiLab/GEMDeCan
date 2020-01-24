suppressMessages({
  if(!require(BiocManager)){
    install.packages("BiocManager")
  }
  if(!require(tximport)){
    BiocManager::install("tximport")
    require(tximport)
  }
  if(!require(org.Hs.eg.db)){
    BiocManager::install("org.Hs.eg.db")
    require(org.Hs.eg.db)
  }
  if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  }
  if(!require(Organism.dplyr)){
    BiocManager::install("Organism.dplyr")
    require(Organism.dplyr)
  }
})

setwd(snakemake@params[[1]])
sample = basename(getwd())

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
src <- src_ucsc("Homo sapiens")
k <- keys(src, keytype = "tx_id")
tx2gene<- select(src, keys = k, columns = c("tx_name","symbol"), keytype = "entrez")
tx2gene<-tx2gene[,-1]
txi <- tximport("quant.sf", type = "salmon", tx2gene = tx2gene)
txi_data<-as.data.frame(cbind(Gene=rownames(txi[[1]]), txi[[1]]))
colnames(txi_data) = c("Gene",sample)
write.table(txi_data, "quantif.txt", sep="\t", quote=F, row.names = F)


