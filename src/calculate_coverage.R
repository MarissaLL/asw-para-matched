#!/usr/bin/env Rscript

library(tidyverse)

#############
# FUNCTIONS #
#############

calculate_coverage <- function(tags_file){
  my_tags <- read_tsv(tags_file, 
                      comment = "#", 
                      col_names = FALSE)
  individual <- gsub("^([^\\.]+).+", "\\1", basename(tags_file))
  
  kept_loci <- unique(filter(my_tags, X8 != 1 & X9 != 1)$X2)
  
  total_reads <- length(unique(filter(my_tags, X2 %in% kept_loci & X5 != "")$X5))
  
  final_coverage <- filter(my_tags, X2 %in% kept_loci & X5 != "") %>% 
    group_by(X2) %>% 
    summarise(reads = length(unique(X5)))
  
  tibble(individual = individual,
         final_coverage_mean = mean(final_coverage$reads),
         n_reads = total_reads)
}

###########
# GLOBALS #
###########

tags_file <- snakemake@input[["tags"]]
output_file <- snakemake@output[["csv"]]
log_file <- snakemake@log[[1]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

cov_stats <- calculate_coverage(tags_file)

write_csv(cov_stats, output_file)
