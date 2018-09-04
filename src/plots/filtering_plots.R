library(ggridges)
library(scales)
library(tidyverse)

# FIGURE OUT DIFFERENCE WITH GGRIDGES. RENAME POPULATIONS

########################## Reads and coverage ############################
cov_stuff <- read_csv("output/010_config/combined_coverage_ustacks.csv") %>% 
  mutate(pop = str_replace_all(individual, "[[:digit:]]", ""))

ggplot(data = cov_stuff, mapping = aes(x = n_reads, y = final_coverage_mean, colour = pop)) +
  geom_point() + 
  theme_classic()+
  labs(x = "Number of reads", y = "Final mean coverage") +
  geom_hline(yintercept = 8)

################## contamination #################################

contaminants <- read_csv("output/010_config/filtering_stats.csv", 
                         skip = 11, 
                         col_names = c("file","discarded","total","kept","individual")) %>% 
  mutate("pop" = str_replace_all(individual, "[[:digit:]]", ""), 
         discard=discarded/total) 


contam_ord <- contaminants[order(contaminants$pop, contaminants$discard),] %>% 
  mutate(sorted = row.names(contam_ord))


ggplot(data = contam_ord, mapping = aes(x = reorder(individual, as.numeric(sorted)), 
                                          y = discard, 
                                          fill = pop)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = percent_format(),
                     limits = c(0,1)) +
  theme_classic()+ 
  theme(axis.ticks.x.bottom = element_blank(), axis.text.x.bottom = element_blank(), axis.line.x.bottom = element_blank()) +
  labs(x = "Sample", y = 'Percentage of reads discarded due to contaminants')

# position = 'fill'
# labels = percent_format()

###### Contamination stats
# Total number of reads, min, max per sample
sum(contam_ord$total)
min(contam_ord$total)
max(contam_ord$total)

# Reads discarded
sum(contam_ord$discarded)

# Samples with fewer than 1 million reads kept
sum(contam_ord$kept < 1000000)

#

################################# GC content ############################


get_mean_gc <-function(x) {
  my_mean <- read_delim(x, delim = "\t",  col_names = c("GC","count")) %>% 
    filter(GC == "#Mean") %>% 
    get("count",.) %>% 
    as.numeric()
  return(tibble(mean_gc=my_mean))
}


file_list<-list.files("./output/021_filtered/gc_hist", full.names = TRUE)
names(file_list) <- sub(".txt", "", basename(file_list))

combined_gc_means <- lapply(file_list, FUN = get_mean_gc) %>% 
  bind_rows(.id = "individual")

q90 <- quantile(combined_gc_means$mean_gc, 0.9)

combined_gc_means %>% 
  mutate("keep" = mean_gc < q90)





individual_order <-  filter(combined_gc_means,!grepl("GBSNEG", individual)) %>% 
  arrange(mean_gc) %>% 
  select(individual) %>% 
  unlist()





combined_gc <- lapply(file_list, FUN = read_delim, delim = "\t", comment = "#", col_names = c("GC","count")) %>% 
  bind_rows(.id = "individual")

gc_plot_data <- filter(combined_gc,!grepl("GBSNEG", individual)) %>% 
  left_join(y = combined_gc_means) %>% 
  mutate("pop" = gsub("[[:digit:]]+", "", individual)) %>% 
  mutate("individual" = factor(individual, levels = individual_order)) %>% 
  mutate("keep" = mean_gc < q90)



ggplot(combined_gc_means,aes(x = individual, y = mean_gc)) +
  geom_point() +
  geom_hline(yintercept = q90)


########################## ggridges #################
#pal <- RColorBrewer::brewer.pal(n = 4,name = "OrRd")

ggplot(gc_plot_data, aes(x = GC, y = individual, height = count, fill = keep)) +
  facet_wrap(~pop, scales = "free") +
  scale_y_discrete(expand = expand_scale(mult = c(0.05,0.21))) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  geom_density_ridges(stat = "identity", scale = 12, alpha = 0.5 ) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "#bfbfbf", colour = "#bfbfbf"),
        axis.ticks.y.left = element_blank(), 
        axis.line.y.left = element_blank(), 
        axis.text.y = element_blank())

