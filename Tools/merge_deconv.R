sample = read.table(snakemake@input[[1]], sep="\t", header=TRUE)

if (!file.exists(snakemake@params[[1]])){
merge = data.frame(t(sample))
  } else{
  merge = t(read.table(snakemake@params[[1]], sep="\t", header=TRUE))
  merge = rbind(merge, sample[2,, drop=FALSE])
}

write.table(merge, snakemake@params[[1]], sep= "\t", quote=F, row.names = F)
