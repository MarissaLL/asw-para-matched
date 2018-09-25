#!/usr/bin/env Rscript

library(scales)
library(tidyverse)

#############
# FUNCTIONS #
#############

extract_readlength <- function(my_file){
  my_data <- read_delim(my_file, delim = '\t,', col_names = TRUE)
  indiv_name <- my_file %>% 
    str_replace_all("output/021_filtered/readlength_hist/", "") %>% 
    str_replace_all("_readlength_hist.txt","")
  data_with_id <- my_data %>% 
    mutate(Individual = indiv_name)
  return(data_with_id)
}

###########
# Globals #
###########

data_dir <- 'output/021_filtered/readlength_hist'

#my_file <-  "output/021_filtered/readlength_hist/L33_readlength_hist.txt"

########
# Main #
########

# list all the readlength files (one for each individual)
key_files <- list.files(data_dir,
                        recursive = FALSE,
                        pattern = "readlength_hist.txt",
                        full.names = TRUE) 

# Extract the data from the readlength files
readlength_list <- lapply(key_files, extract_readlength)

# Convert to a tibble
readlengths <- bind_rows(readlength_list)

# Total read counts at each position (across all individuals)
all_combined <- readlengths %>% 
  group_by(`#Length` ) %>% 
  summarise(sum(Count)) %>% 
  rename(Length =`#Length`, Count = `sum(Count)`)

# Plot
ggplot(all_combined, aes(x = Length, y = Count)) +
  scale_y_continuous(limits=c(-1,25000000), labels = comma) + # to zoom in use limits=c(-1,30000000),
  scale_x_continuous(limits = c(-1,92)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Length of ASW sequence', y = 'Number of reads') 
