#!/usr/bin/env Rscript

library(tidyverse)

###########
# GLOBALS #
###########

popmap <- snakemake@input[["popmap"]]
reftable_indivs <- snakemake@output[[1]]
  
########
# MAIN #
########
# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Read in and re-format the data
indivs_file <- read_delim(popmap, delim = '\t',col_names = FALSE) %>% 
  mutate('sex' = 9) %>% 
  rename('indiv' = X1, 'pop' = X2) %>% 
  select(indiv, sex, pop)

# Write out the data
write_delim(indivs_file, path = reftable_indivs, delim = '\t')

# Log session info
sessionInfo()
