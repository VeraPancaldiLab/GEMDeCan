#' ---
#' title: Analysis on deconvolution results
#' author: Julien Pernet
#' date: 
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

args = commandArgs(trailingOnly = TRUE)
cells <- read.table(args[1], header = T, sep = "\t", row.names = 1)
cells.mat <- as.matrix(cells)
#'
#' # Data visualization 
#'
heatmap(cells.mat, main = "Estimation of cell type abundance", col = viridis::inferno(50), scale = "column")

cor = Hmisc::rcorr(cells.mat)
cor$r[is.nan(cor$r) ] = 0
heatmap(cor$r, main = "Correlation between cell types", col = viridis::inferno(50))
#'
#' # Full raw data
#'
knitr::kable(cells)