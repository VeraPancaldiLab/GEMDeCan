sample <- snakemake@input

merge <- read.table(sample[[1]], sep = "\t", header = T)
for (i in 2:length(sample)) {
  x <- read.table(sample[[i]], sep = "\t", header = T)
  merge <- merge(merge, x, by="Gene")
}

write.table(merge, snakemake@output[[1]], sep = "\t", quote = F, row.names = F)
