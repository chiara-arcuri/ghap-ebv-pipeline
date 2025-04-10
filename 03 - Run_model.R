#############################################################
# Title: Run single-step GBLUP model with GHap
# Authors: Chiara Arcuri and Marco Milanesi
# Date: 10/04/2025
# Description:
#   This script performs a single-step GBLUP analysis by:
#     - Loading genotype, phenotype, and pedigree data
#     - Loading and preparing the inverse of the H matrix (genomic + pedigree matrix)
#     - Fitting a linear mixed model using GHap to predict breeding values (EBVs)
#     - Saving the model environment for future use
#
# Requirements:
#   - GHap package
#   - Matrix package
#   - Input files:
#     - Genotype data in PLINK format
#     - Phenotype file (e.g., renf90GENOBU.dat)
#     - Pedigree file (e.g., renadd01GENOBU.ped)
#     - H-inverse matrix file (e.g., Hinv_ebv.txt)
#############################################################

# Prepare the environment
options(stringsAsFactors = FALSE)
options(digits = 15)
options(scipen = 999)
rm(list = ls())  
gc()  

# Load required libraries
library(Matrix) 
library(GHap)  

# Load PLINK genotype data (in GHap format)
plink <- ghap.loadplink("output")

# Load phenotype data from a text file
df <- read.table(file = "renf90GENOBU.dat", header = FALSE,
                 col.names = c("id","TC","HYS","LC","DO","pheno"))

# Load pedigree data from a text file
ped <- read.table(file = "renadd01GENOBU.ped", header = FALSE, 
                  col.names = c("id","sire","dam"))

# Load the inverse of the H matrix (this includes genomic and pedigree information)
Hinv <- read.table(file = "Hinv_ebv.txt", header = FALSE)

# Combine and sort unique levels (from V1 and V2 columns of the Hinv matrix)
common_levels <- sort(unique(c(Hinv$V1, Hinv$V2)))

# Map levels to numerical indices
row_indices <- match(Hinv$V1, common_levels)
col_indices <- match(Hinv$V2, common_levels)

# Create a sparse matrix using the H-inverse data
sparse_mat <- sparseMatrix(
  i = row_indices,
  j = col_indices,
  x = Hinv$V3,  # Values from the matrix
  symmetric = TRUE,
  dimnames = list(common_levels, common_levels)  # Set row and column names
)

# Display the sparse matrix
print(sparse_mat)

# Run the single-step GBLUP model
df$rep <- df$id  # Add the replication variable to the phenotype data

# Fit the GBLUP model using the formula with random effects for id and rep
model3 <- ghap.lmm(formula = pheno ~ 1 + TC + HYS + LC + DO + (1|id) + (1|rep),
                   data = df,
                   covmat = list(id = sparse_mat, rep = NULL),
                   invcov = TRUE,  # Use the inverse covariance matrix
                   errors = FALSE,  # Suppress error output
                   vcp.initial = list(id = 84706, rep = 52892, Residual = 105940),  # Initial variance components
                   vcp.estimate = FALSE, 
                   em.reml = 3)  

# Save the model environment
save.image(file = 'Env_GHap.RData')
