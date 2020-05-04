require(immunedeconv)

t1 <- as.matrix(read.table(snakemake@input[[1]], header = TRUE, row.names = 1))
res <- deconvolute(t1, "quantiseq")
res[, -1] <- round(res[, -1], 3)
res <- t(res)
write.table(res, snakemake@output[[1]], sep = "\t", quote = F, row.names = T, col.names = F)
