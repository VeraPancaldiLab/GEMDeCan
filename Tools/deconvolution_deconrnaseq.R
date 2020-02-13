require(DeconRNASeq)

signature = read.table(snakemake@params[[2]], header = TRUE)
t1 = read.table(snakemake@input[[1]], header = TRUE, row.names = 1)
# signature need to be formated to the correct input format for deconrna, as it is picky
rownames(signature)= signature[,1]
res = DeconRNASeq(t1, signature[,-1])
Fres = as.data.frame(res$out.all)
Fres= cbind(Cell_type=colnames(t1), Fres)
Fres[,-1] = round(Fres[,-1], 3)
write.table(Fres, snakemake@output[[1]], sep="\t", quote=F, row.names = F, col.names = T)