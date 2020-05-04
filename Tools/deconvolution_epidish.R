require(EpiDISH)

signature <- as.matrix(read.table(snakemake@params[[1]], header = TRUE, row.names = 1))
t1 <- read.table(snakemake@input[[1]], header = TRUE, row.names = 1)

res <- epidish(t1, signature, method = "RPC")
Fres <- as.data.frame(res$estF)
Fres <- cbind(Sample = rownames(Fres), Fres)
Fres[, -1] <- round(Fres[, -1], 3)

write.table(Fres, snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)
