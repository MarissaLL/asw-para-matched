#!/usr/bin/env Rscript

library(tidyverse)

###########
# GLOBALS #
###########


para_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]


# dev
para_file <-  'sample_catalog.csv'

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Read in and reformat data
para <- read.csv(para_file) %>% 
  as_tibble() %>% 
  filter(Head =="TRUE") %>% 
  mutate_at("Parasitoid", replace_na, "FALSE") %>%       
  mutate("Parasitism" = if_else(Parasitoid=="TRUE", "parasitized","non-parasitized")) %>% 
  mutate("pop_para" = paste(Location,Parasitism, sep = "_")) %>% 
  mutate("pop_para_past" = paste(Location, Parasitism, Pasture, sep = "_"))  %>% 
  select(Individual, Location, Pasture, Parasitism, pop_para, pop_para_past)

# Export reformatted data
write_delim(para, output_file, delim = '\t')

# Log session info
sessionInfo()