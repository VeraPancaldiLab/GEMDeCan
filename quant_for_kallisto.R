suppressMessages({
  if(!require(BiocManager)){
    install.packages("BiocManager")
  }
  if(!require(tximport)){
    BiocManager::install("tximport")
    require(tximport)
  }
  if(!require(tidyverse)){
    BiocManager::install("tidyverse")
    require(tidyverse)
  }
  if(!require(Organism.dplyr)){
    BiocManager::install("Organism.dplyr")
    require(Organism.dplyr)
  }
})
snakemake@params[[1]] %>% setwd()

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
src <- src_ucsc("Homo sapiens")
k <- keys(src, keytype = "tx_id")
tx2gene<- select(src, keys = k, columns = c("tx_name","symbol"), keytype = "entrez")
tx2gene<-tx2gene[,-1]
txi <- tximport("abundance.h5", type = "kallisto", tx2gene = tx2gene, txIn=TRUE)
txi_data<-as.data.frame(cbind(Gene=rownames(txi[[1]]), txi[[1]]))
write.table(txi_data, "quantif.txt", sep="\t", quote=F, row.names = F)
