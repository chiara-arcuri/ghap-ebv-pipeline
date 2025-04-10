#############################################################
# Title: Run single-step GBLUP model
# Authors: Chiara Arcuri and Marco Milanesi
# Date: 10/04/2025
# Description:
#   This script performs a single-step GBLUP analysis by:
#     - Fitting a linear mixed model using GHap
# Requirements:
#   - GHap package 
#   - GHap-formatted PLINK genotype data ("output")
#   - Phenotype file: BufRenumEBV.dat
#   - Pedigree file: BufRenumEBVs.ped
#############################################################

# Prepare the environment
options(stringsAsFactors = FALSE)
options(digits = 15)
options(scipen = 999)
rm(list = ls())
gc()

# Load required library
library(GHap)

# Load PLINK data in GHap format
plink <- ghap.loadplink("output")

# Load phenotype data
df <- read.table(file = "BufRenumEBV.dat", header = FALSE, 
                 col.names = c("id","TC","PHYS","LC","CLASSDO","SIRECODPROV",
                               "pheno","GP","PP","G","P","RESA"))

# Load pedigree data (id, sire, dam)
ped <- read.table(file = "BufRenumEBVs.ped", header = FALSE, 
                  col.names = c("id", "sire", "dam"))

# Compute the genomic relationship matrix (GRM)
K <- ghap.kinship(plink)

# Prepare list of all individuals for the H matrix
ids <- unique(c(ped$id, ped$sire, ped$dam, colnames(K)))

# Compute the inverse of the H matrix (blended genomic and pedigree)
Hinv <- ghap.getHinv(K = K, ped = ped, include = ids)

# Add replication column for random effects
df$rep <- df$id

# Fit single-step GBLUP model
model3 <- ghap.lmm(
  formula = pheno ~ 1 + (1|id) + (1|rep),
  data = df,
  covmat = list(id = Hinv, rep = NULL),
  invcov = TRUE
)

# Save R environment with fitted model
save.image("model_env.RData")
