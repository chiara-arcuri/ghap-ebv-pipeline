###########################################################
# Title: Phasing Pipeline with PLINK and Eagle
# Author: Chiara Arcuri and Marco Milanesi
# Date: 10/04/2025
# Description: This script performs genotype phasing 
#   chromosome-by-chromosome using PLINK and Eagle.
#   Designed for high-performance execution with multithreading.
#
# Requirements:
# - PLINK v1.9 
# - Eagle v2.4.1 
###########################################################

### 1. Environment setup
options(stringsAsFactors = FALSE)  
options(digits = 15)               
options(scipen = 999)              
rm(list = ls())  
gc()             

### 2. Parameters (customize these as needed)
ncores <- 48                     
pflie <- <Input PLINK file>         
outf <- "OUTFILE"                
nchr <- 24                      
xchr <- 25                       

### 3. Loop over each chromosome and perform phasing
for (i in 1:nchr) {
  
  cat("Processing chromosome", i, "...\n")  
  
  ## Step 1: Extract single chromosome using PLINK
  plink_cmd <- paste0(
    "plink --chr-set ", nchr,
    " --allow-no-sex --nonfounders",
    " --bfile ", pflie,
    " --chr ", i,
    " --make-bed --out tmp"
  )
  system(plink_cmd)
  
  ## Step 2: Run Eagle for phasing
  eagle_cmd <- paste0(
    "eagle --bfile tmp",
    " --geneticMapFile=USE_BIM",             
    " --chromX ", xchr,
    " --maxMissingPerSnp 1",
    " --maxMissingPerIndiv 1",
    " --numThreads ", ncores,
    " --outPrefix=", outf, "_", i, "_phased"
  )
  system(eagle_cmd)
  
  ## Step 3: Clean up temporary files
  system("rm tmp*")
  
  cat("Chromosome", i, "done.\n\n")
}

cat("Phasing pipeline completed for all chromosomes.\n")
