suppressPackageStartupMessages({
  require(ChAMP)
  require(minfi)
   require(zoo)
})


# load data
message(">>>  Loading data")
wd <- snakemake@input[[1]]
epicAnnot <- snakemake@params[[1]]
pheno <- snakemake@params[[2]] 
cores = snakemake@params[[3]]
array_type = snakemake@params[[4]]

Phenotype <- read.csv(pheno, header = T, sep = ",")
anno_epic <- read.csv(epicAnnot, as.is = TRUE, skip = 7)
anno_epic <- anno_epic[, c("CHR", "MAPINFO", "Name", "UCSC_RefGene_Name")]

RGSet <- read.metharray.exp(wd, recursive = T, force = TRUE)

# if input file is IlluminaHumanMethylation450k, change  #
#         outType = "IlluminaHumanMethylation450k"       #
#         outType = "IlluminaHumanMethylationEPIC"       #

if (array_type == "450K"){
  RGSet <- convertArray(RGSet, outType = "IlluminaHumanMethylation450k") # because OF COURSE this one has no capital letter
}else {
  RGSet <- convertArray(RGSet, outType = paste0("IlluminaHumanMethylation", array_type))
}

colnames(RGSet) <- Phenotype$Sample_Name


# QC removes bad samples ----------
# Samples with median Meth and Umneth signals below this cutoff will be labelled ‘bad’
message(">>> Processing data")
MSet <- preprocessRaw(RGSet)
qc <- getQC(MSet)
meds <- (qc$uMed + qc$mMed) / 2
keepIndex <- which(meds > 10.5)
Good_samples <- colnames(RGSet)[which(index(colnames(RGSet)) %in% keepIndex)]

RawBeta <- champ.load(wd, arraytype = array_type)

# Normalization ------------
message(">>> Computing beta values")
BMIQ <- champ.norm(beta = RawBeta$beta, arraytype = array_type, cores = cores, resultsDir = paste0(snakemake@output[[1]],"/BMIQ_Normalization/"), method = "BMIQ")

# final output for snakemake is written at the end
message(">>> Writing final output")
write.table(BMIQ, snakemake@output[[1]], sep="\t", quote=F, row.names = T, col.names = T)
