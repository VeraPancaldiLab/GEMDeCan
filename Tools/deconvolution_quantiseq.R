require(immunedeconv)

t1 = as.matrix(read.table(snakemake@input[[1]], header = TRUE, row.names = 1))
sample = basename(snakemake@params[[1]])
res = deconvolute(t1, "quantiseq")
colnames(res) = c("cell_type", sample)
res[,2] = round(res[,2], 2)
write.table(res, snakemake@output[[1]], sep="\t", quote=F, row.names = F)

