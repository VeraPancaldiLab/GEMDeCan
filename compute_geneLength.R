suppressMessages({
  require(tximport)
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Organism.dplyr)
  require(tidyverse)
  require(GenomicRanges)
})

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
src <- src_ucsc("Homo sapiens")
k <- keys(src, keytype = "tx_id")
tx2gene <- select(src, keys = k, columns = c("symbol"), keytype = "entrez")
tx_by_gene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns = "gene_id")
gene_length <- data.frame("entrez" = tx_by_gene$gene_id, "length" = width(tx_by_gene))
gene_length <- inner_join(tx2gene, gene_length, by = "entrez")[, -1]
write_tsv(gene_length, "gene_length.txt")