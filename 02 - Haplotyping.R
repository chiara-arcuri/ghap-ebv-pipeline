###########################################################
# Title: GHap haplotyping Processing 
# Authors: Chiara Arcuri and Marco Milanesi 
# Date: 10/04/2025
# Description:
#   This script processes phased haplotype data (Eagle output)
#   using the GHap package in R. It includes:
#     - Decompression of phased files
#     - Conversion to GHap binary format
#     - Haplotype block generation
#     - Haplotype genotype matrix creation
#     - Conversion to PLINK format
#
# Requirements:
#   - R packages: GHap, R.utils
#   - Input: .haps.gz and .sample files from Eagle
###########################################################

### 1. Environment setup
options(stringsAsFactors = FALSE)  
options(digits = 15)               
options(scipen = 999)              
rm(list = ls())  
gc()             

### 2. Load required packages
# install.packages("GHap")
library(GHap)

### 3. Set working directory and parameters                  
output_folder <- "analyses"          
nchr <- 24                             

### 4. Unzip phased haplotype files
for (i in 1:nchr) {
  file <- paste0("OUTFILE_", i, "_phased.haps.gz")
  gunzip(file) 
}

### 5. Convert Eagle output to GHap format
ghap.oxford2phase(
  haps.files   = paste0("OUTFILE_", 1:nchr, "_phased.haps"),
  sample.files = paste0("OUTFILE_", 1:nchr, "_phased.sample"),
  out.file     = file.path(output_folder, "output")
)

### 6. Compress GHap phase data
ghap.compress(
  input.file = "output",
  out.file   = file.path(output_folder, "output")
)

### 7. Load GHap phase data
phase <- ghap.loadphase(
  input.file = file.path(output_folder, "output")
)

### 8. Generate blocks of 4 markers sliding 1 markers at a time
blocks.mkr <- ghap.blockgen(phase, windowsize = 4, 
                            slide = 1, unit = "marker")

### 9. Create haplotype genotype matrix
ghap.haplotyping(
  object  = phase,
  blocks  = blocks.mkr,
  outfile = file.path(output_folder, "output"),
  binary  = TRUE
)

### 10. Load haplotype genotypes
haplo <- ghap.loadhaplo(
  file.path(output_folder, "output")
)

### 11. Convert haplotype genotypes to PLINK format
ghap.hap2plink(
  haplo,
  outfile = file.path(output_folder, "output")
)

cat("GHap processing completed.\n")
