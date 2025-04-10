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

### 1. Setup environment
options(stringsAsFactors = FALSE)
options(digits = 15)
options(scipen = 999)
rm(list = ls())
gc()

# Load required packages
library(GHap)
library(openxlsx)
library(ggplot2)
library(ggExtra)

### 2. Load genotype data
plink <- ghap.loadplink("genobu_in")

### 3. Load EBVs
ebv <- read.xlsx(xlsxFile = "ebv_all.xlsx", colNames = TRUE)
ebv_gen <- ebv[which(ebv$IDA %in% plink$id), c(8,5)]
ebv_gen <- ebv_gen[order(match(ebv_gen$IDA, plink$id)), ]
refblup <- ebv_gen$solBLUP
names(refblup) <- ebv_gen$IDA

### 4. Calculate kinship matrix
K <- ghap.kinship(plink, ncores = 30)

### 5. Estimate variance explained by each marker
mkrblup <- ghap.varblup(object = plink, gebv = refblup, covmat = K, ncores = 30, verbose = TRUE)

mkrblup$percVAR <- mkrblup$pVAR * 100
write.table(mkrblup[mkrblup$percVAR > 0.5, ], "percVAR_over05.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Manhattan plots
pdf("pVAR.pdf", 10, 7)
ghap.manhattan(data = mkrblup, chr = "CHR", bp = "BP", y = "pVAR")
dev.off()

pdf("percVAR.pdf", 10, 7)
ghap.manhattan(data = mkrblup, chr = "CHR", bp = "BP", y = "percVAR",
               main = "Explained Variance\n1 marker window",
               xlab = "Chromosome", ylab = "% VAR")
dev.off()

### 6. Calculate GEBVs and compare to true EBVs
gebv <- ghap.profile(object = plink, score = mkrblup, ncores = 30, verbose = TRUE)
hist(gebv$SCORE)

# Correlation and visualization
cor.test(gebv$SCORE, refblup)
plot(gebv$SCORE, refblup); abline(0, 1)
