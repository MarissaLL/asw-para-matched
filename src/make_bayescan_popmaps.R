#!/usr/bin/env Rscript

library(tidyverse)

###########
# GLOBALS #
###########


popmap_file <- snakemake@input[["popmap"]]
para_datafile <- snakemake@input[["para_data"]]
popmap_4pops <- snakemake@output[["popmap_4pops"]]
popmap_para <- snakemake@output[["popmap_para"]]
popmap_2pops <- snakemake@output[["popmap_2pops"]]
popmap_island <- snakemake@output[["popmap_island"]]
popmap_ruakura <- snakemake@output[["popmap_ruakura"]]
popmap_ruakura_poa <- snakemake@output[["popmap_ruakura_poa"]]
popmap_lincoln <- snakemake@output[["popmap_lincoln"]]
popmap_invermay <- snakemake@output[["popmap_invermay"]]
log_file <- snakemake@log[[1]]

#dev
# popmap_file <- 'output/060_pop_genet/r0.8_filtered_popmap.txt'
# para_datafile <- 'output/010_config/tidy_sample_info.tsv'


########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# Read in the input data
popmap_pop <- read_delim(popmap_file, delim = '\t', col_names = FALSE)
para_data <- read_delim(para_datafile, delim = '\t')

# Write out the same popmap with a different name to make running the downstream stuff together easier
write_delim(x = popmap_pop, path = popmap_4pops, delim = '\t', col_names = FALSE)

# Filter the popmap to only include Ruakura (merged) and Lincoln
RvL <- popmap_pop %>% 
  mutate(X2 = str_sub(X2, 1, 1)) %>% 
  filter(X2 == c('R', 'L'))

# Write new popmap for Ruakura vs Lincoln
write_delim(x = RvL, path = popmap_2pops, delim = '\t', col_names = FALSE)

# Replace population column of full filtered popmap with parasitism status info
indiv_status <-  para_data$Parasitism
names(indiv_status) <- para_data$Individual
para_map <- popmap_pop
para_map$X2 <- indiv_status[para_map$X1]

# Write new popmap with parasitism status info for all populations combined
write_delim(x = para_map, path = popmap_para, delim = '\t', col_names = FALSE)

# Split the popmap with parasitism info for all populations into separate populations
map_I <- para_map %>% 
  filter(str_detect(X1,"^I") == TRUE)

map_L <- para_map %>% 
  filter(str_detect(X1,"^L") == TRUE)

map_R <- para_map %>% 
  filter(str_detect(X1,"^R(?!p)") == TRUE)

map_Rpoa <- para_map %>% 
  filter(str_detect(X1,"^Rpoa") == TRUE)

# Write popmaps for all the separate populations
write_delim(x = map_I, path = popmap_invermay, delim = '\t', col_names = FALSE)
write_delim(x = map_L, path = popmap_lincoln, delim = '\t', col_names = FALSE)
write_delim(x = map_R, path = popmap_ruakura, delim = '\t', col_names = FALSE)
write_delim(x = map_Rpoa, path = popmap_ruakura_poa, delim = '\t', col_names = FALSE)

# Replace population information with which island samples were from
popmap_NS <- popmap_pop %>% 
  mutate(island = case_when(X2 == "I" | X2 == "L" ~ "S",
                            X2 == "R" | X2 == "Rpoa" ~ "N")) %>% 
  select(X1,island)

# Write out the popmap for comparison of north and south island
write_delim(x = popmap_NS, path = popmap_island, delim = '\t', col_names = FALSE)



# Log session info
sessionInfo()