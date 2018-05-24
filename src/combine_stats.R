#!/usr/bin/env Rscript

library(tidyverse)


#############
# FUNCTIONS #
#############

get_mean_gc <-function(x) {
  my_mean <- read_delim(x, delim = "\t",  col_names = c("GC","count")) %>% 
    filter(GC == "#Mean") %>% 
    get("count",.) %>% 
    as.numeric()
  return(tibble(mean_gc=my_mean))
}

parse_bbduk_stats <- function(x){
my_tibble <- read_delim(x, 
                        delim = "\t", 
                        skip = 0, 
                        n_max = 3, 
                        col_names = c("variable", "value"))

spread(my_tibble, variable, value, convert = TRUE)
}

###########
# GLOBALS #
###########

stats_files <-snakemake@input[["stats_file"]]
  
gc_files<-snakemake@input[["gc_file"]]


output_file <- snakemake@output[["stats_file"]]
log_file <- snakemake@log[[1]]
 


########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

bbduk_stats_list <- lapply(stats_files, parse_bbduk_stats)

bbduk_stats <- bind_rows(bbduk_stats_list) %>% 
  rename(`#Discarded` = `#Matched`) %>% 
  mutate(`#Kept` = `#Total`-`#Discarded`,
         `#Individual` = sub(".fq.gz", "", basename(`#File`)))

names(gc_files) <- sub(".txt", "", basename(gc_files))

combined_gc_means <- lapply(gc_files, FUN = get_mean_gc) %>% 
  bind_rows(.id = "#Individual")

merged_stats <- left_join(bbduk_stats, combined_gc_means)

write_csv(merged_stats, path = output_file)
  
# write log
sessionInfo()
