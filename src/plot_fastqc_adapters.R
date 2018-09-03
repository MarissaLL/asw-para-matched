#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

#############
# FUNCTIONS #
#############

# Read in the section of the fastqc data file that contains information about adapter
read_adapter_section <- function(file_name){
  whole_file <- read_lines(file_name)
  start_line <- grep('^>>Adapter Content', whole_file)
  section_breaks <- grep('^>>END_MODULE', whole_file)
  stop_line <- sort(section_breaks[section_breaks > start_line])[[1]]
  adapter_section <- read_tsv(paste(whole_file[(start_line+1):(stop_line-1)], collapse = '\n')) %>% 
    mutate(position = as.numeric(substr(`#Position`, 1, 2)))
  return(adapter_section)
}

# Access the fastqc data file without needing to unzip the folder it is in
read_zip <- function(file_list){
  my_files <- unzip(file_list, list = TRUE)
  desc <- my_files[[grep('fastqc_data.txt', my_files$Name), "Name"]]
  zip <- unz(file_list, desc)
  read_adapter_section(zip)
}

# Go through all the zipped folders in the directory and extract the adapter information into a format that can be plotted
extract_plot_data <- function(data_dir){
  key_files <- list.files(data_dir,
                          recursive = FALSE,
                          pattern = ".zip",
                          full.names = TRUE) 
  
  names(key_files) <-  sub(".zip", "", basename(key_files)) 
  
  adapter_stats_list <- lapply(key_files, read_zip)
  for_plotting <- melt(adapter_stats_list, id = 'position') %>% 
    filter(variable == 'Illumina Universal Adapter')
  return(for_plotting)
}

###########
# Globals #
###########

before_data_dir <- 'output/022_fastqc/before_filter'
after_data_dir <- 'output/022_fastqc/after_filter'

########
# Main #
########

# Extract data from the fastqc data files within all the zipped folders in the fastqc output directory and put it in a format that can be plotted
before_plot <- extract_plot_data(before_data_dir)
after_plot <- extract_plot_data(after_data_dir)


# Plot adapter content by position
ggplot(mapping = aes(x = as.numeric(position), y = as.numeric(value), group = L1)) +
  ylim(c(0,75)) +
  geom_line(data = after_plot, colour = "black") +
  geom_line(data = before_plot, colour = "#990000") +
  labs(x = "Position along read", y = "Adapter content for all reads within individual (%)") +
  theme_classic() 
  
