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


txdb = GenomicFeatures::makeTxDbFromGFF(snakemake@params[[2]])

k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")
txi = tximport::tximport("abundance.h5", type = "kallisto", txIn = TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE)
table <- txi$counts %>% round(4) %>% as.data.frame() %>% tibble::rownames_to_column()
names(table)[1]="Gene" 
table  %>% write_tsv("quantif.txt", quote=F)
