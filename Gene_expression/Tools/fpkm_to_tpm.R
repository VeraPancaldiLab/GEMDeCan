# Convert FPKM values to transcripts per million (TPM).

fpkm2tpm <- function(fpkm) {
  tpm <- exp(log(fpkm) - log(sum(fpkm, na.rm = T)) + log(1e6))
  tpm[which(is.na(tpm))] <- 0
  return(tpm)
}

# TPM normalization

TPM_normalization <- function(data) {
  
  # TPM normalization
  
  TPM_data <- t(t(data) * 1e6 / apply(data, 2, sum))
  return(TPM_data)
}
