require(immunedeconv)

t1 = as.matrix(read.table(snakemake@input[[1]], header = TRUE, row.names = 1))
sample = basename(snakemake@params[[1]])
# replace mcp counter with snakemake@params
res = deconvolute(t1, "mcp_counter", feature_types="HUGO_symbols", 
                  probesets=read.table(snakemake@params[[2]],sep="\t",stringsAsFactors=FALSE,colClasses="character"), 
                  genes=read.table(snakemake@params[[3]], sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))
colnames(res) = c("cell_type", sample)
res[,2] = round(res[,2], 2)
write.table(res, snakemake@output[[1]], sep="\t", quote=F, row.names = F)