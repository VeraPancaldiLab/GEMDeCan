suppressPackageStartupMessages({
  library(EpiDISH)
})

message(">>> Reading data")
t1 <- read.table(snakemake@input[[1]], header = TRUE, sep="\t", row.names = 1)
# t1 <- read.table("/home/julien/Documents/lungMethylation/pipeline/output/results/beta_values.txt", header = TRUE, sep="\t", row.names = 1)

signatures_path = snakemake@params[[1]]
signature_files <- list.files(signatures_path, full.names = T)
signature_files <- setdiff(signature_files, signature_files[dir.exists(signature_files)])
data(centDHSbloodDMC.m)
#########################
###### Deconvolution ----
#########################

message(">>> Running deconvolution on EpiDish.")

res <- epidish(t1, centDHSbloodDMC.m, method = "RPC", maxit = 100)
res_epi <- as.data.frame(res$estF)
res_epi <- cbind(Sample = rownames(res_epi), res_epi)
res_epi[, -1] <- round(res_epi[, -1], 5)
colnames(res_epi) = paste0("Epidish_", colnames(res_epi))
colnames(res_epi)[1] = "Row.names"
final_deconv = res_epi

for (i in 1:length(signature_files)){
  # get signature info
  signature <- as.matrix(read.table(signature_files[i], header = TRUE, row.names = 1, sep="\t"))
  sign_name = basename(signature_files[i]) 
  sign_name = gsub(".txt", "", sign_name)
  message(">>> Running ", sign_name, " signature.")
  
  # Do deconv epidish
  res <- epidish(t1, signature, method = "RPC", maxit = 100)
  res_epi <- as.data.frame(res$estF)
  res_epi <- cbind(Sample = rownames(res_epi), res_epi)
  res_epi[, -1] <- round(res_epi[, -1], 5)
  colnames(res_epi) = paste0("Epidish_", sign_name, "_", colnames(res_epi))
  final_deconv = merge(final_deconv, res_epi[,-1], by.x = "Row.names", by.y = "row.names")

}

rownames(final_deconv) = final_deconv[,1]
final_deconv = final_deconv[,-1]

message(">>> Deconvolution sucessfull. Writing output.")
write.table(final_deconv, snakemake@output[[1]], sep = "\t", quote = F, row.names = F, col.names = T)
