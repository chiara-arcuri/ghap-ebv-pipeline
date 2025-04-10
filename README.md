This repository provides a complete pipeline to perform haplotyping-based analyses using the GHap R package, with PLINK-formatted genotype data as input. 
The pipeline includes data phasing, file conversion, kinship matrix construction, SNP effect estimation, GEBV prediction, cross-validation, and external validation.

## Requirements
- R (≥ 4.0)
- Packages: GHap 4.0.0.12, PLINK v1.9, Eagle v2.4.1

## Pipeline Summary
01 - Phasing_pipeline.R : Splits data by chromosome and performs phasing using Eagle.
02 - Haplotyping.R : Converts phased .haps and .sample files to GHap format. Compresses and loads data, defines 1-SNP sliding haplotype blocks.
03 - Run_model.R : Performs a single-step GBLUP analysis
04 - EBV_analysis.R : Produces GEBVs
05 - Validation_test.R : Calculates the correlation between the predicted GEBVs and reference BLUPs




