#!/usr/bin/env Rscript

library(tidyverse)

###########
# GLOBALS #
###########


ped_file <- snakemake@input[["ped_file"]]
para_data <- snakemake@input[["para_info"]]
output_ped <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]


########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Load the ped that needs population and parasitism data added to the first and sixth columns
ped_data <- read_delim(ped_file, delim = '\t', col_names = FALSE)

# Generate population column based on sample name
pop_para_ped <- ped_data %>% 
  mutate(X1 = str_replace_all(X2, "[[:digit:]]", ""))

# Read in info about which individuals are parasitised
para_info <-  read_delim(para_data, delim = '\t') %>% 
  select(Individual, Parasitism) %>% 
  mutate(State = if_else(Parasitism == "non-parasitized", 1 , 2))

# Add parasitism info to the ped file
indiv_status <-  para_info$State
names(indiv_status) <- para_info$Individual
pop_para_ped$X6 <- indiv_status[pop_para_ped$X2]

# Write out a new ped file
write_delim(pop_para_ped, output_ped , delim = '\t', col_names = FALSE )

# Log session info
sessionInfo()

