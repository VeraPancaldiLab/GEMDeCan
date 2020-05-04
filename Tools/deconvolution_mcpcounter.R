install.packages(c("devtools", "curl"))
require(devtools)
install_github("ebecht/MCPcounter", ref = "master", subdir = "Source")
require(MCPcounter)

t1 <- as.matrix(read.table(snakemake@input[[1]], header = TRUE, row.names = 1))
res_mcp <- MCPcounter.estimate(t1,
  featuresType = "HUGO_symbols",
  genes = read.table(snakemake@params[[1]], sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
)
res_mcp[] <- round(res_mcp[], 3)
res_mcp <- t(res_mcp)
res_mcp <- rbind(Sample = colnames(res_mcp), res_mcp)
write.table(res_mcp, snakemake@output[[1]], sep = "\t", quote = F, row.names = T, col.names = F)
