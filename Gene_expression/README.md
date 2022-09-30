---
noteId: "e830fdf03fec11edb2b3a33f1817fe71"
tags: []

---

# How to run the pipeline (sudo/root permissions are mandatory)
`sudo snakemake --cores 8 all`

# From where take the generated CIBERSORTx signatures
`results/cibersortx_signatures`
And add them to the EpiDISH and DeconRNASeq explained in the next section

# Where to add more signatures to EpiDISH and DeconRNASeq
`scripts/deconvolution/signatures`

# Where to find the inputs folder
`/mnt/SERVER-CRCT-STORAGE/CRCT21/Private/Miguel/GEMDeCan/inputs`

# Where to find the credentials.txt file
`/mnt/SERVER-CRCT-STORAGE/CRCT21/Private/Miguel/GEMDeCan/credentials.txt`

# Where to find the final tables
`results/all_deconvolutions`
