#!/usr/bin/env Rscript

library(tidyverse)

###########
# GLOBALS #
###########

pop_2pop <- snakemake@input[['popmap_2pops']] 
pop_r <- snakemake@input[['popmap_ruakura']]
pop_rpoa <-  snakemake@input[['popmap_ruakura_poa']]
pop_l <- snakemake@input[['popmap_lincoln']] 
pop_i <- snakemake@input[['popmap_invermay']]

i_2pop <- snakemake@output[['indivs_2pops']] 
i_r <- snakemake@output[['indivs_ruakura']]
i_rpoa <- snakemake@output[['indivs_ruakura_poa']]
i_l <- snakemake@output[['indivs_lincoln']]
i_i <- snakemake@output[['indivs_invermay']]

log_file <- snakemake@log[[1]]

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


 
indivs_2pops <- read_delim(pop_2pop, delim = '\t', col_names = FALSE) %>% 
  select(X1)

indivs_ruakura <- read_delim(pop_r, delim = '\t', col_names = FALSE) %>% 
  select(X1)

indivs_ruakura_poa <- read_delim(pop_rpoa, delim = '\t', col_names = FALSE) %>% 
  select(X1)

indivs_lincoln <- read_delim(pop_l, delim = '\t', col_names = FALSE) %>% 
  select(X1)

indivs_invermay <- read_delim(pop_i, delim = '\t', col_names = FALSE) %>% 
  select(X1)



write_delim(indivs_2pops, i_2pop, col_names = FALSE)
write_delim(indivs_ruakura, i_r, col_names = FALSE)
write_delim(indivs_ruakura_poa, i_rpoa, col_names = FALSE)
write_delim(indivs_lincoln, i_l, col_names = FALSE)
write_delim(indivs_invermay, i_i, col_names = FALSE)

# Log session info
sessionInfo()