#!/usr/bin/env Rscript

library(tidyverse)


#############
# Functions #
#############

parse_bbduk_stats <- function(x){
my_tibble <- read_delim(x, 
                        delim = "\t", 
                        skip = 0, 
                        n_max = 3, 
                        col_names = c("variable", "value"))

spread(my_tibble, variable, value, convert = TRUE)
}



# Globals

stats_files <- snakemake@input[["stats_file"]]
output_file <- snakemake@output[["stats_file"]]
log_file <- snakemake@log[[1]]

# Main

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


bbduk_stats_list <- lapply(stats_files, parse_bbduk_stats)

bbduk_stats <- bind_rows(bbduk_stats_list) %>% 
  rename(`#Discarded` = `#Matched`) %>% 
  mutate(`#Kept` = `#Total`-`#Discarded`,
         `#Individual` = sub(".fq.gz","",basename(`#File`)))

write_csv(bbduk_stats, path = output_file)
  
  # write log
  sessionInfo()
