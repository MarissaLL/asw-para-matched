#!/usr/bin/env Rscript

library(SNPRelate)
library(tidyverse)


#############
# FUNCTIONS #
#############

###########
# GLOBALS #
###########

VCFfile <- snakemake@input[[1]]
GDSfile <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]
  
#dev
 # "output/050_stacks_pops/r0/populations.snps.vcf"  

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


snpgdsVCF2GDS(VCFfile, GDSfile, method = "biallelic.only" )