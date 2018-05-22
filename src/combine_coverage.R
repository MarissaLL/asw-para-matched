#!/usr/bin/env Rscript

library(tidyverse)


#############
# FUNCTIONS #
#############

read_coverage_stats <- function(x){
  my_tibble <- read_delim(coverage_files, 
                          delim = ",", 
                          skip = 0)
}

###########
# GLOBALS #
###########

coverage_files <- snakemake@input[["coverage_file"]]
output_file <- snakemake@output[["coverage_file"]]
log_file <- snakemake@log[[1]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")
print(coverage_files)

cov_stats_list <- lapply(coverage_files, read_coverage_stats)

cov_stats <- bind_rows(cov_stats_list) 

write_csv(cov_stats, path = output_file)

# write log
sessionInfo()
