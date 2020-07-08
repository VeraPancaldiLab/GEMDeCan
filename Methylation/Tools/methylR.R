suppressPackageStartupMessages({
  require(ChAMP)
  require(minfi)
  require(BEclear)
  require(factoextra)
  require(zoo)
  require(dplyr)
  require(DNAmArray)
  require(sva)
  require(Hmisc)
  require(pheatmap)
  require(limma)
})

# load data
wd <- snakemake@input[[1]] # "/home/julien/Documents/lungMethylation/Methylomes_Lucille_Vera/DATA/Rawdata"
epicAnnot <- snakemake@params[[1]] # "/home/julien/Documents/lungMethylation/Methylomes_Lucille_Vera/MethylationEPIC_v-1-0_B4.csv"
pheno <- snakemake@params[[2]] # "Phenotype.csv"
cores = snakemake@params[[3]]
dir.create("Outputs/results/plots", showWarnings = FALSE, recursive = TRUE)

Phenotype <- read.csv(pheno, header = T, sep = ",")
anno_epic <- read.csv(epicAnnot, as.is = TRUE, skip = 7)
anno_epic <- anno_epic[, c("CHR", "MAPINFO", "Name", "UCSC_RefGene_Name")]

RGSet <- read.metharray.exp(wd, recursive = T, force = TRUE)

# if input file is IlluminaHumanMethylation450k, change  #
#         outType = "IlluminaHumanMethylation450k"       #
#         outType = "IlluminaHumanMethylationEPIC"       #

RGSet <- convertArray(RGSet, outType = "IlluminaHumanMethylationEPIC")
colnames(RGSet) <- Phenotype$Sample_Name


# QC removes bad samples ----------
# Samples with median Meth and Umneth signals below this cutoff will be labelled ‘bad’
MSet <- preprocessRaw(RGSet)
qc <- getQC(MSet)
plotQC(qc)
meds <- (qc$uMed + qc$mMed) / 2
keepIndex <- which(meds > 10.5)
Good_samples <- colnames(RGSet)[which(index(colnames(RGSet)) %in% keepIndex)]

RawBeta <- champ.load(wd, arraytype = "EPIC")
champ.QC(beta = RawBeta$beta, pheno = RawBeta$pd$Sample_Group, resultsDir = "Outputs/Rawbeta_QC")

#screeplot(RGSet)

# Normalization by different methods ------------

BMIQ <- champ.norm(beta = RawBeta$beta, arraytype = "EPIC", cores = cores, resultsDir = "Outputs/BMIQ_Normalization/", method = "BMIQ")
gc()

BetaNoob <- preprocessNoob(RGSet, offset = 15, dyeCorr = TRUE, verbose = FALSE, dyeMethod = "single") %>%
  getBeta(.) %>%
  champ.filter(beta = ., pd = Phenotype, arraytype = "EPIC")
BetaNoob$beta <- BetaNoob$beta[(rownames(BetaNoob$beta) %in% rownames(RawBeta$beta)), ]

BetaFN <- preprocessFunnorm(RGSet, nPCs = 4) %>%
  getBeta(.) %>%
  champ.filter(beta = ., pd = Phenotype, arraytype = "EPIC")
BetaFN$beta <- BetaFN$beta[(rownames(BetaFN$beta) %in% rownames(RawBeta$beta)), ]

BetaFNm <- preprocessFunnorm.DNAmArray(RGSet, nPCs = 4, keepCN = FALSE) %>%
  getBeta(.) %>%
  champ.filter(beta = ., pd = Phenotype, arraytype = "EPIC")
BetaFNm$beta <- BetaFNm$beta[(rownames(BetaFNm$beta) %in% rownames(RawBeta$beta)), ]

# function for plot beta densities
plot.beta.densities <- function(beta, title, name) {
  if (!is.null(dim(beta))) {
    densities <- apply(beta, 2, function(x) {
      density(x, na.rm = TRUE)
    })
    xmax <- max(sapply(densities, function(d) {
      max(d$x)
    }))
    xmin <- min(sapply(densities, function(d) {
      min(d$x)
    }))
    ymax <- max(sapply(densities, function(d) {
      max(d$y)
    }))
    png(name)
    plot(NA, xlim = c(xmin, xmax), ylim = c(0, ymax), main = title, ylab = "")
    colors <- rainbow(10)
    for (i in 1:ncol(beta)) {
      lines(densities[[i]], col = colors[i %% 10 + 1])
    }
    dev.off()
  } else if (length(beta) > 1) {
    pdf(name)
    plot(density(beta, na.rm = TRUE), main = title)
    dev.off()
  }
}
plot.beta.densities(RawBeta$beta, "Densities of raw values per sample", "Outputs/results/plots/rawValues.png")
plot.beta.densities(BMIQ, "Densities of BMIQ normalization per sample", "Outputs/results/plots/BMIQ_norm.png")
plot.beta.densities(BetaNoob$beta, "Densities of Noob normalization per sample", "Outputs/results/plots/noob_norm.png")
plot.beta.densities(BetaFN$beta, "Densities of Function normalization per sample", "Outputs/results/plots/function_norm.png")
plot.beta.densities(BetaFNm$beta, "Densities of Function DNAmArray normalization per sample", "Outputs/results/plots/DNAmArray_norm.png")

# Heatmap *** -----

cor_BMIQ <- rcorr(BMIQ, type = "pearson")$r

cor_Noob <- rcorr(BetaNoob$beta, type = "pearson")$r

cor_FN <- rcorr(BetaFN$beta, type = "pearson")$r

cor_FNm <- rcorr(BetaFNm$beta, type = "pearson")$r

pheatmap(cor_BMIQ, main = "Heatmap BMIQ", filename = "Outputs/results/plots/heatmap_bmiq.png")
pheatmap(cor_Noob, main = "Heatmap Noob", filename = "Outputs/results/plots/heatmap_noob.png")
pheatmap(cor_FN, main = "Heatmap Function", filename = "Outputs/results/plots/heatmap_function.png")
pheatmap(cor_FNm, main = "Heatmap Function DNAmArray", filename = "Outputs/results/plots/heatmap_DNAmArray.png")

# DMPs analysis ****** -------

ct <- factor(Phenotype$Cell_line)
design <- model.matrix(~ 0 + ct)
colnames(design) <- levels(ct)
fit <- lmFit(BMIQ, design)
contrasts <- makeContrasts(HL60.basal - HL60WT.basal, HL60.Ag120 - HL60WT.basal, MOLM.basal - MOLMWT.basal, MOLM.AGI120 - MOLMWT.basal, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend = TRUE)
summary(decideTests(fit2))

move <- function(data, cols, ref, side = c("before", "after")) {
  if (!requireNamespace("dplyr")) {
    stop("Make sure package 'dplyr' is installed to use function 'move'")
  }
  side <- match.arg(side)
  cols <- rlang::enquo(cols)
  ref <- rlang::enquo(ref)
  if (side == "before") {
    dplyr::select(data, 1:!!ref, -!!ref, -!!cols, !!cols, dplyr::everything())
  } else {
    dplyr::select(data, 1:!!ref, -!!cols, !!cols, dplyr::everything())
  }
}

FitList <- list()
for (i in 1:ncol(contrasts)) {
  FitList[[i]] <- topTable(fit2, coef = i, number = nrow(BMIQ)) %>%
    mutate(ID = rownames(.)) %>%
    move(., ID, logFC, side = "before") %>%
    filter(adj.P.Val < 0.05)
}

### Different cell lines ### ----
Transformed <- data.frame(t(RawBeta$beta))
Split <- split(Transformed, Phenotype$Cell_line)
Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
Split <- do.call(cbind, Split)
rownames(Split) <- rownames(RawBeta$beta)

# Then I need to actually annotate each one of these comparison topTables to get beta values

dbList <- list()
for (i in 1:ncol(contrasts)) {
  dB <- with(data.frame(Split), eval(parse(text = colnames(contrasts)[[i]])))
  dB <- data.frame(dB = dB, ID = rownames(Split))
  dbList[[i]] <- dB
}

# Filter by thresholds ----

dbList <- lapply(dbList, function(x) filter(x, abs(dB) > 0.2))
FitListFinal <- list()
for (i in 1:length(FitList)) {
  A1 <- FitList[[i]]
  A1 <- filter(A1, ID %in% dbList[[i]]$ID)
  A1 <- A1 %>% .[rev(order(.$t)), ]
  FitListFinal[[i]] <- A1
}

AnnoEpic <- move(anno_epic, Name, CHR, "before") %>% move(., UCSC_RefGene_Name, CHR, "before")

for (i in 1:ncol(contrasts)) {
  FitListFinal[[i]] <- merge(AnnoEpic, FitListFinal[[i]], by.x = "Name", by.y = "ID") %>%
    merge(., BMIQ, by.x = "Name", by.y = "row.names") %>%
    merge(., dbList[[i]], by.x = "Name", by.y = "ID")
  rownames(FitListFinal[[i]]) <- FitListFinal[[i]]$Name
  FitListFinal[[i]] <- FitListFinal[[i]][, -1]
}

for (i in 1:length(FitListFinal)) {
  names(FitListFinal) <- colnames(summary(decideTests(fit2)))
  write.table(data.frame("Outputs/CpGs" = rownames(FitListFinal[[i]]), FitListFinal[[i]]), file = paste(names(FitListFinal)[i], ".txt"), sep = "\t", quote = F, row.names = F)
}

# final output for snakemake is written at the end
write.table(BMIQ, snakemake@output[[1]], sep="\t", quote=F, row.names = T, col.names = T)
