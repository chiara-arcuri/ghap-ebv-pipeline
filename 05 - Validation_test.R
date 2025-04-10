#############################################################
# Title: GHap Hold-Out GEBV Validation (400 animals)
# Authors: Chiara Arcuri and Marco Milanesi
# Date: 10/0472025
# Description:
#   This script performs:
#     - Training/test split (holdout of 400 animals)
#     - EBV variance partitioning on training set
#     - Prediction of GEBVs for all animals (including test set)
#     - Correlation analysis and visualization
#
# Requirements:
#   - GHap, openxlsx
#   - Input PLINK data (with and without 400 animals)
#   - EBV file (xlsx format)
#############################################################

### 1. Setup
options(stringsAsFactors = FALSE)
options(digits = 15)
options(scipen = 999)
rm(list = ls())
gc()

library(GHap)
library(openxlsx)

### 2. Load complete genotype dataset
plink_all <- ghap.loadplink("genobu_in")

### 3. [One-time step] Create test set of 400 animals
# remove <- sample(plink_all$id, 400, replace = FALSE)
# write(paste("IT", remove), file = "remove_400.txt")

# create filtered PLINK dataset
plink_cmd <- paste(
  "plink",
  "--bfile genobu_in",
  "--remove remove_400.txt",
  "--make-bed",
  "--out genobu_in_filtered"
)
system(plink_cmd)

### 4. Load filtered training dataset (4027 animals)
plink <- ghap.loadplink("genobu_in_filtered")

### 5. Load EBV and match to training set
ebv <- read.xlsx("ebv_all.xlsx", colNames = TRUE)
ebv_gen <- ebv[ebv$IDA %in% plink$id, c(8,5)]
ebv_gen <- ebv_gen[order(match(ebv_gen$IDA, plink$id)), ]
refblup <- ebv_gen$solBLUP
names(refblup) <- ebv_gen$IDA

### 6. Compute kinship and marker BLUPs on training set
K <- ghap.kinship(plink, ncores = 40)
mkrblup <- ghap.varblup(object = plink, gebv = refblup, covmat = K, ncores = 30, verbose = TRUE)

### 7. Predict GEBVs for training set
gebv <- ghap.profile(object = plink, score = mkrblup, ncores = 40, verbose = TRUE)

# Check accuracy on training
cor_train <- cor(gebv$SCORE, refblup)
print(paste("Training correlation:", cor_train)) 
plot(gebv$SCORE, refblup); abline(0, 1)

### 8. Predict GEBVs on full dataset (4427 animals)
test_gebv <- ghap.profile(object = plink_all, score = mkrblup, ncores = 50, verbose = TRUE)

### 9. Load EBV for full dataset and compute full correlation
ebv_gen_all <- ebv[ebv$IDA %in% plink_all$id, c(8,5)]
ebv_gen_all <- ebv_gen_all[order(match(ebv_gen_all$IDA, plink_all$id)), ]
refblup_all <- ebv_gen_all$solBLUP
names(refblup_all) <- ebv_gen_all$IDA

cor_all <- cor(test_gebv$SCORE, refblup_all)
print(paste("Overall correlation (4427):", cor_all)) 
plot(test_gebv$SCORE, refblup_all); abline(0,1)

### 10. Evaluate performance on test set (400 animals)
removed_ids <- read.table("remove_400.txt", header = FALSE, stringsAsFactors = FALSE)$V2
removed_ids <- as.character(removed_ids)

# Merge and filter to 400 test subjects
merged <- merge(test_gebv, ebv_gen_all, by.x = "ID", by.y = "IDA")
merged_400 <- merged[merged$ID %in% removed_ids, ]

cor_400 <- cor(merged_400$SCORE, merged_400$solBLUP)
print(paste("Test correlation (400 held-out):", cor_400)) # ≈ 0.928
plot(merged_400$SCORE, merged_400$solBLUP); abline(0, 1)

### 11. Save results
write.table(merged_400, "GEBV_testset_400.txt", sep = "\t", row.names = FALSE, quote = FALSE)
