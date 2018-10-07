#!/usr/bin/env Rscript

library(ggridges)
library(scales)
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

###########
# Globals #
###########

filter_stats_file <- "output/010_config/filtering_stats.csv"
coverage_file <- "output/010_config/combined_coverage_ustacks.csv"
gc_content_dir <- "output/021_filtered/gc_hist"

filtered_popmap <-  "output/010_config/filtered_popmap.txt"


########################## Reads and coverage ############################
cov_stuff <- read_csv(coverage_file) %>% 
  mutate(pop = str_replace_all(individual, "[[:digit:]]", ""))

ggplot(data = cov_stuff, mapping = aes(x = n_reads, y = final_coverage_mean, colour = pop)) +
  geom_point() + 
  theme_classic()+
  labs(x = "Number of reads", y = "Final mean coverage") +
  geom_hline(yintercept = 8)

################## contamination #################################

# Read in data on number of reads (total, kept and discarded)
contaminants <- read_csv(filter_stats_file, 
                         skip = 11, 
                         col_names = c("file","discarded","total","kept","individual")) %>% 
  mutate("pop" = str_replace_all(individual, "[[:digit:]]", ""), 
         discard=discarded/total) 

# Order bars by population and then by % reads discarded
contam_ord <- contaminants[order(contaminants$pop, contaminants$discard),] 
contam_ord$sorted <- row.names(contam_ord)

# Assign colours to the different populations
pop_colours <- RColorBrewer::brewer.pal(4, "PiYG")
#pop_colours <- c("pink", "pink3", "#B8E186", "#4DAC26")
names(pop_colours) <- c("I", "L", "R", "Rpoa")

# Plot 
ggplot(data = contam_ord, mapping = aes(x = reorder(individual, 
                                                    as.numeric(sorted)), 
                                          y = discard, 
                                          fill = pop)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = percent_format(),
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_fill_manual(values = pop_colours, 
                    labels = c("Invermay", "Lincoln", "Ruakura", "Ruakura (Poa)")) +
  theme_classic()+ 
  theme(axis.ticks.x.bottom = element_blank(), 
        axis.text.x.bottom = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) +
  labs(x = "Sample", 
       y = 'Percentage of reads discarded due to contaminants', 
       fill = "Population")

################# Number of reads before removing contamination ####################

# Order bars by population and then by total number of reads
contam_ord_tot <- contaminants[order(contaminants$pop, contaminants$total),] 
contam_ord_tot$sorted <- row.names(contam_ord_tot)


# Plot 
ggplot(data = contam_ord_tot, mapping = aes(x = reorder(individual, 
                                                    as.numeric(sorted)), 
                                        y = total, 
                                        fill = pop)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = comma,
                     expand = c(0,0)) +
  scale_fill_manual(values = pop_colours, 
                    labels = c("Invermay", "Lincoln", "Ruakura", "Ruakura (Poa)")) +
  theme_classic()+ 
  theme(axis.ticks.x.bottom = element_blank(), 
        axis.text.x.bottom = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) +
  labs(x = "Sample", 
       y = 'Number of reads in sample', 
       fill = "Population")




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

file_list<-list.files(gc_content_dir, full.names = TRUE)
names(file_list) <- sub(".txt", "", basename(file_list))

combined_gc_means <- lapply(file_list, FUN = get_mean_gc) %>% 
  bind_rows(.id = "individual")

q90 <- quantile(combined_gc_means$mean_gc, 0.9)

individual_order <-  filter(combined_gc_means,!grepl("GBSNEG", individual)) %>% 
  arrange(mean_gc) %>% 
  select(individual) %>% 
  unlist()


combined_gc <- lapply(file_list, FUN = read_delim, delim = "\t", comment = "#", col_names = c("GC","count")) %>% 
  bind_rows(.id = "individual")

gc_plot_data <- filter(combined_gc,!grepl("GBSNEG", individual)) %>% 
  left_join(y = combined_gc_means) %>% 
  mutate("pop" = gsub("[[:digit:]]+", "", individual)) %>% 
  mutate("Population" = case_when(pop == "I" ~ "Invermay", pop == "R" ~ "Ruakura", pop == "L" ~ "Lincoln", pop == "Rpoa" ~ "Ruakura (Poa)")) %>% 
  mutate("individual" = factor(individual, levels = individual_order)) %>% 
  mutate("keep" = mean_gc < q90)

final_kept <- read_delim(filtered_popmap, delim = '\t', col_names = FALSE) %>% 
  rename(individual = X1) %>% 
  select(individual) %>% 
  mutate(filtered = "kept")

all_filter_for_ggridges <- gc_plot_data %>% 
  left_join(y = final_kept) %>% 
  mutate_at("filtered", replace_na, "discarded") %>% 
  mutate("individual" = factor(individual, levels = individual_order)) 


ggplot(combined_gc_means,aes(x = individual, y = mean_gc)) +
  geom_point() +
  geom_hline(yintercept = q90)


########################## ggridges #################
#pal <- RColorBrewer::brewer.pal(n = 4,name = "OrRd")

# Ggridges indicating samples thrown out because of GC content
ggplot(gc_plot_data, aes(x = GC, y = individual, height = count, fill = keep)) +
  facet_wrap(~Population, scales = "free") +
  scale_y_discrete(expand = expand_scale(mult = c(0.05,0.21))) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  geom_density_ridges(stat = "identity", scale = 12, alpha = 0.5 ) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "#bfbfbf", colour = "#bfbfbf"),
        strip.text = element_text(size = 12),
        axis.ticks.y.left = element_blank(), 
        axis.line.y.left = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_text(size = 14)) + 
  labs(x = "GC content", y = "Sample")


# Ggridges indicating samples thrown out because of both remaining read count and GC content
ggplot(all_filter_for_ggridges, aes(x = GC, y = individual, height = count, fill = filtered)) +
  facet_wrap(~Population, scales = "free") +
  scale_y_discrete(expand = expand_scale(mult = c(0.05,0.21))) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  geom_density_ridges(stat = "identity", scale = 12, alpha = 0.5 ) +
  theme_classic() +
  theme(strip.background = element_rect(fill = "#bfbfbf", colour = "#bfbfbf"),
        strip.text = element_text(size = 14),
        axis.ticks.y.left = element_blank(), 
        axis.line.y.left = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(x = "GC content (%)", y = "Sample")
