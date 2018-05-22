library(tidyverse)

touch_flagfile <- function(sample_name, outdir){
  write.table(data.frame(), file=paste(outdir, sample_name, sep = "/"), col.names=FALSE)
}


combined_counts <- read_csv(snakemake@input[[1]], 
                            col_names = c("Sample","Lines"))

filtered_samples <- combined_counts %>% 
  filter(Lines/4 > 1e6) %>% 
  select(Sample)

sapply(unique(filtered_samples$Sample),
       FUN = touch_flagfile, 
       outdir =snakemake@params[["outdir"]])
