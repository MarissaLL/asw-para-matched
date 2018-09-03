#!/usr/bin/env Rscript

library(adegenet)
library(tidyverse)


###########
# Globals #
###########

plink_file <- snakemake@input[["plink_file"]]
log_file <- snakemake@log[[1]]
whitelist <- snakemake@output[["whitelist"]]
popmap <-  snakemake@output[["popmap"]]

# plink_file <- 'output/060_pop_genet/plink.raw'
########
# Main #
########
# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


# Import the data 
SNP_data <- read.PLINK(plink_file)

# Extract individual names
ind_names <- indNames(SNP_data) %>% 
  as.tibble() %>% 
  rename(individual = value) %>% 
  mutate(pop = str_replace_all(individual, "[[:digit:]]", ""))

# Write the filtered popmap
write_delim(ind_names, 
            path = popmap, 
            delim = "\t", 
            col_names = FALSE)


# Extract locus ID from SNP names
locus_names <- locNames(SNP_data) %>% 
  as.tibble() %>% 
  mutate("locus_ID" = str_split(value, 
                                "_", 
                                simplify = TRUE)[,1]) %>% 
  mutate("SNP_ID" = str_split(value,
                              "_",
                              simplify = TRUE)[,2])

# Separate locus ID
whitelist_IDs <- select(locus_names, locus_ID, SNP_ID) 

 
# Write the locus whitelist
write_delim(whitelist_IDs, 
            path = whitelist, 
            delim = "\t", 
            col_names = FALSE)

sessionInfo()
