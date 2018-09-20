#!/usr/bin/env Rscript

library(data.table)
library(scales)
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

I11_gc_file <- 'fastqc_examples/I11_fastqc/out04'



########
# Main #
########

# Extract data from the fastqc data files within all the zipped folders in the fastqc output directory and put it in a format that can be plotted
before_plot <- extract_plot_data(before_data_dir) %>% 
  mutate(step = "Before filtering")

after_plot <- extract_plot_data(after_data_dir) %>% 
  mutate(step = "After filtering")

combined_fq_data <- rbind(before_plot, after_plot)


# Plot adapter content by position
ggplot(data = combined_fq_data, 
       mapping = aes(x = as.numeric(position), y = as.numeric(value), colour = step, group = interaction(L1, step))) +
  ylim(c(0,75)) +
  geom_line() +
  labs(x = "Position along read", y = "Adapter content for all reads within individual (%)") +
  theme_classic() +
  theme(axis.title = element_text(size = 14), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 12),
        legend.text=element_text(size=12)) +
  scale_colour_manual(name = "", values = c("black","#990000"), guide = guide_legend(reverse=TRUE))



# Read in an example file with an odd gc distribution
gc_content <- read_delim(I11_gc_file, delim = '\t', skip = 2)

# Parameters for the theoretical distribution. Currently just fudged, FIX THIS
mean_value <-  46.8 #gc_content$`#GC Content`[max(gc_content$Count)]
sd_value <- 12
n_value <-  6000000


# Plot the actual and theoretical gc distributions
ggplot(gc_content, aes(x = `#GC Content`, y = Count)) + 
  geom_line() +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  labs(x = "Mean GC content (%)",
       y = "GC Count") +
  stat_function(fun = function(x, mean, sd, n, bw){ 
    dnorm(x = x, mean = mean, sd = sd) * n * bw}, 
    args = c(mean = mean_value, 
             sd = sd_value, 
             n = n_value, 
             bw = 1), 
    colour = "#990000") +
  theme(axis.title = element_text(size = 14), 
             axis.title.y = element_text(margin = margin(r = 10)),
             axis.text = element_text(size = 12))
