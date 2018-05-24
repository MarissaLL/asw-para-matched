#!/usr/bin/env Rscript

library(tidyverse)


###########
# GLOBALS #
###########

coverage_files <- snakemake@input[["coverage_file"]]
full_popmap_file <- snakemake@input[["full_popmap"]]
output_file <- snakemake@output[["coverage_file"]]
output_popmap <- snakemake@output[["filtered_popmap"]]
log_file <- snakemake@log[[1]]

#dev
# coverage_files <-list.files('output/041_ustacks_coverage',pattern = '.csv',full.names = TRUE)
# full_popmap_file <-'output/010_config/full_popmap.txt'

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


cov_stats_list <- lapply(coverage_files, read_csv)

cov_stats <- bind_rows(cov_stats_list) 

full_popmap <- read_tsv(full_popmap_file, col_names = FALSE)

filtered_popmap <- filter(full_popmap, X1 %in% unique(cov_stats$individual))

write_csv(cov_stats, path = output_file)

write_tsv(filtered_popmap, path = output_popmap,col_names = FALSE)

# write log
sessionInfo()
