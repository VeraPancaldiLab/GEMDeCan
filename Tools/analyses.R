#' ---
#' title: Analysis on deconvolution results
#' author: Julien Pernet
#' date:
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
knitr::opts_chunk$set(echo = F, warning = F, message = F)
#'
#' # Data visualization
#'
# heatmap(cells.mat, main = "Estimation of cell type abundance", col = viridis::inferno(50), scale = "column")


args <- commandArgs(trailingOnly = TRUE)
cells <- read.table(args[1], header = T, sep = "\t", row.names = 1)
cells.mat <- as.matrix(cells)

cells.mat.no.zero <- cells.mat[, -which(colSums(cells.mat) == 0)]

pheatmap::pheatmap(t(cells.mat.no.zero), scale = "row", main = "Estimation of cell type abundance", fontsize_row = 10 - nrow(cells.mat.no.zero) / 5, fontsize_col = 10 - ncol(cells.mat.no.zero) / 20)

cor_samples <- Hmisc::rcorr(t(cells.mat))
# cor_samples <- cor_samples$r[cor_samples$P < 0.05, cor_samples$P < 0.05]
# cor$r[is.nan(cor$r)] <- 0
# corrplot::corrplot(cor$r, tl.cex = 25 / ncol(cor$r), title = "Correlation between cell types", mar = c(0, 0, 1, 0))

corrplot::corrplot(cor_samples$r, tl.cex = 25 / ncol(cor_samples$r), title = "Correlation between cell types", mar = c(0, 0, 1, 0))


pheatmap::pheatmap(cor_samples$r, main = "Correlation between samples", fontsize_row = 10 - nrow(cells.mat.no.zero) / 15, fontsize_col = 10 - ncol(cor_samples$r) / 15)

deconvolution_methods <- unique(sapply(colnames(cells.mat.no.zero), function(colname) {
  stringr::str_split(colname, "_")[[1]][[1]]
}))

for (method in deconvolution_methods) {
  mat <- cells.mat.no.zero[, grepl(method, colnames(cells.mat.no.zero))]
  cor_cell_types <- Hmisc::rcorr(mat)
  pheatmap::pheatmap(cor_cell_types$r, main = stringr::str_c(method, " correlation between cell types"), fontsize_row = 10 - nrow(cor_cell_types$r) / 5, fontsize_col = 10 - ncol(cor_cell_types$r) / 10)
}

# cor <- Hmisc::rcorr(cells.mat)
# stop(cor$r)
# cor <- cor[!is.nan(cor$P) & cor$P < 0.05]
# cor$r[is.nan(cor$r)] <- 0
# corrplot::corrplot(cor$r, tl.cex = 25 / ncol(cor$r), title = "Correlation between cell types", mar = c(0, 0, 1, 0))
# pheatmap::pheatmap(cor$r, main = "Correlation between cell types", fontsize_row = 3, fontsize_col = 3)
# heatmap(cor$r, main = "Correlation between cell types", col = RColorBrewer::brewer.pal(n = 50, "RdBu"))
#'
#' # Full raw data
#'
knitr::kable(cells)
