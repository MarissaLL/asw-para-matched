library(scales)
library(tidyverse)

label_func <- function(x) {gsub("[[:alpha:]]", "", x)}

popstats_combined_Mm = read.csv("output_20180614/030_optim/stats_Mm/popstats_combined.csv")
samplestats_combined_Mm = read.csv("output_20180614/030_optim/stats_Mm/samplestats_combined.csv")

popstats_combined_n = read.csv("output_20180614/030_optim/stats_n/popstats_combined.csv")
samplestats_combined_n = read.csv("output_20180614/030_optim/stats_n/samplestats_combined.csv")


# Change format of the data for ggplot
samples_Mm_gathered <-  samplestats_combined_Mm %>% 
  gather(key, value, -c(sample,stacks_run,m,M,n,rep))

# Remove all points where M was also being varied at the same time from the plot of m
plot_m <- samples_Mm_gathered %>% 
  filter(M == "M2")

# Remove all points where m was also being varied at the same time from the plot of M
plot_M <- samples_Mm_gathered %>% 
  filter(m == "m3")

facet_labels <- c(assembled_loci = "Assembled loci",
                  polymorphic_loci = "Polymorphic loci",
                  snps = "SNPs")

# Add r0.8 points
popstats_Mm_gathered <- popstats_combined_Mm %>%
  gather(key, value, -c(stacks_run,m,M,n,rep))

# Remove M varied from plot of m for r0.8 values
m_r0.8 <- popstats_Mm_gathered %>% 
  filter(M == "M2")

# Remove m varied from plot of M for r0.8 values
M_r0.8 <- popstats_Mm_gathered %>% 
  filter(m == "m3")

# Set jitter and dodge positions for ggplot
jd <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)
d <- position_dodge(width = 0.9)


# Plot assembled loci, polymorphic loci and snps at different values of m
ggplot(plot_m, mapping = aes(x = m, 
                             y = value, 
                             colour = rep)) +
  geom_boxplot(position = d, 
               mapping = aes(colour = rep))+
  geom_point(position = jd)+
  geom_point(data = m_r0.8, 
             mapping = aes(x = m, 
                           y = value, 
                           group = rep), 
             position = d, 
             colour = 'black') +
  facet_wrap(~key, 
             nrow = 3, 
             scales = "free_y",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(strip.placement = "outside", 
        strip.background = element_blank()) +
  labs(x = "Value of m parameter" ) 

# Plot M
ggplot(plot_M, mapping = aes(x = M, 
                             y = value, 
                             colour = rep)) +
  geom_boxplot(position = d, 
               mapping = aes(colour = rep))+
  geom_point(position = jd)+
  geom_point(data = M_r0.8, 
             mapping = aes(x = M, 
                           y = value, 
                           group = rep), 
             position = d, 
             colour = 'black') +
  facet_wrap(~key, 
             nrow = 3, 
             scales = "free_y",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(strip.placement = "outside", 
        strip.background = element_blank()) +
  labs(x = "Value of M parameter")

# Gather n
samples_n_gathered <- samplestats_combined_n %>% 
  gather(key, value, -c(sample,stacks_run,m,M,n,rep))

popstats_n_gathered <- popstats_combined_n %>% 
  gather(key, value, -c(stacks_run,m,M,n,rep))


# Plot n
ggplot(samples_n_gathered, mapping = aes(x = n, 
                             y = value, 
                             colour = rep)) +
  geom_boxplot(position = d, 
               mapping = aes(colour = rep))+
  geom_point(position = jd)+
  geom_point(data = popstats_n_gathered, 
             mapping = aes(x = n, 
                           y = value, 
                           group = rep), 
             position = d, 
             colour = 'black') +
  facet_wrap(~key, 
             nrow = 3, 
             scales = "free_y",
             labeller = as_labeller(facet_labels)) +
  theme_classic() +
  theme(strip.placement = "outside", 
        strip.background = element_blank()) +
  labs(x = "Value of n parameter")

