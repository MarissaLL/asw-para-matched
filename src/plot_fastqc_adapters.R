#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

#############
# FUNCTIONS #
#############

read_adapter_section <- function(file_name){
  whole_file <- read_lines(file_name)
  start_line <- grep('^>>Adapter Content', whole_file)
  section_breaks <- grep('^>>END_MODULE', whole_file)
  stop_line <- sort(section_breaks[section_breaks > start_line])[[1]]
  adapter_section <- read_tsv(paste(whole_file[(start_line+1):(stop_line-1)], collapse = '\n')) %>% 
    mutate(position = as.numeric(substr(`#Position`, 1, 2)))
  return(adapter_section)
}

read_zip <- function(file_list){
  my_files <- unzip(file_list, list = TRUE)
  desc <- my_files[[grep('fastqc_data.txt', my_files$Name), "Name"]]
  zip <- unz(file_list, desc)
  read_adapter_section(zip)
}

###########
# Globals #
###########

data_dir <- 'output/022_fastqc'

########
# Main #
########

# List zipped directories produced by fastqc, containing the data
key_files <- list.files(data_dir,
                        recursive = FALSE,
                        pattern = ".zip",
                        full.names = TRUE) 

names(key_files) <-  sub(".zip", "", basename(key_files)) 

# Extract 
adapter_stats_list <- lapply(key_files, read_zip)

for_plotting <- melt(adapter_stats_list, id = 'position') %>% 
  filter(variable == 'Illumina Universal Adapter')

ggplot(for_plotting, aes(x = as.numeric(position), y = as.numeric(value), group = L1)) +
  geom_line() +
  labs(x = "Position along read", y = "What is this?")
