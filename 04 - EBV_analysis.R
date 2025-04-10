#############################################################
# Title: GHap EBV Pipeline
# Authors: Chiara Arcuri and Marco Milanesi
# Date: 10/04/2025
# Description:
#   This script:
#     - Loads genotype data and EBVs
#     - Calculates kinship matrix
#     - Estimates SNP-wise variance explained
#     - Produces GEBVs
#     - Plots Manhattan and histograms
# Requirements:
#   - GHap
#   - Input PLINK file (binary format)
#   - EBV file (xlsx format)
#############################################################

# Preparation of the environment 
options(stringsAsFactors = FALSE)
options(digits = 15)
options(scipen = 999)
rm(list = ls())  
gc()  

# Load necessary libraries
library(GHap)
library(openxlsx)
library(ggplot2)
library(ggExtra)

# Load PLINK genotype data
plink <- ghap.loadplink("output")

# Load EBV data
ebv <- read.xlsx(xlsxFile = "ebv_all.xlsx", colNames = TRUE)
ebv_gen <- ebv[which(ebv$IDA %in% plink$id), c(8, 5)]

# Ensure that EBV data and PLINK IDs match
identical(as.numeric(ebv_gen$IDA), as.numeric(plink$id))
ebv_gen <- ebv_gen[order(match(ebv_gen$IDA, plink$id)), ]
identical(as.numeric(ebv_gen$IDA), as.numeric(plink$id))

# Extract reference BLUPs and create histogram
refblup <- ebv_gen$solBLUP
names(refblup) <- ebv_gen$IDA
hist(refblup)

# Compute the genomic relationship matrix (kinship)
K <- ghap.kinship(plink, ncores = 40)

# Convert BLUPs of individuals into BLUPs of variants
mkrblup <- ghap.varblup(object = plink, gebv = refblup, covmat = K, ncores = 40, verbose = TRUE)

# Summarize the pVAR values
summary(mkrblup$pVAR)

# Create and save the Manhattan plot for pVAR
pdf("pVAR.pdf", width = 10, height = 7)
ghap.manhattan(data = mkrblup, chr = "CHR", bp = "BP", y = "pVAR")
dev.off()

# Plot percentage of variance explained by SNPs
pdf("percVAR.pdf", width = 10, height = 7)
mkrblup$percVAR <- mkrblup$pVAR * 100
ghap.manhattan(data = mkrblup, chr = "CHR", bp = "BP", y = "percVAR", main = "Explained Variance\nHaplotypes", xlab = "Chromosome", ylab = "% VAR")
dev.off()

# Save variants explaining more than 0.5% of the variance
write.table(x = mkrblup[which(mkrblup$percVAR > 0.5), ], file = "percVAR_over05.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Build GEBVs from variant effects and compare predictions with reference BLUPs
gebv <- ghap.profile(object = plink, score = mkrblup, ncores = 20, verbose = TRUE)
hist(gebv$SCORE)

# Correlation between predicted GEBVs and reference BLUPs
cor.test(gebv$SCORE, refblup)
plot(gebv$SCORE, refblup, main = "Correlation by Haplotypes", xlab = "EBVs computed by GHap", ylab = "EBVs computed by ANASB")
abline(0, 1)
