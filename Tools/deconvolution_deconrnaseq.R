require(DeconRNASeq)

signature <- read.table(snakemake@params[[1]], header = TRUE, row.names = 1)
t1 <- read.table(snakemake@input[[1]], header = TRUE, row.names = 1)
# signature need to be formated to the correct input format for deconrna, as it is picky
res <- DeconRNASeq(t1, signature)
Fres <- as.data.frame(res$out.all)
Fres <- cbind(Sample = colnames(t1), Fres)
Fres[, -1] <- round(Fres[, -1], 3)
write.table(Fres, snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)
