#!/usr/bin/env Rscript

library(scales)
library(SNPRelate)
library(stringr)
library(tidyverse)



#############
# FUNCTIONS #
#############
get_SNP_set <- function(missing_rate, 
                        gds_obj, 
                        minor_allele_freq) {
  my_SNPs <- snpgdsSelectSNP(gdsobj = gds_obj, 
                             autosome.only = FALSE, 
                             remove.monosnp = TRUE, 
                             maf = minor_allele_freq, 
                             missing.rate = missing_rate, 
                             verbose = TRUE)
  my_SNPIDs <- unique(unlist(my_SNPs))
  return(my_SNPIDs)
}



# This gets the sample missing rate for the supplied list of SNPs, and turns it into a tibble with sample names and SNP missing rates
get_sample_missing_rates <- function(SNPlist, 
                                     gds_obj){
  my_sample_missing_rate <- snpgdsSampMissRate(gdsobj = gds_obj, 
                                               snp.id = SNPlist, 
                                               with.id = TRUE)
  
  return(
    tibble(individual = names(my_sample_missing_rate), 
           sample_missing_rate = my_sample_missing_rate))
}

# Counts the kept SNPs following filtering by the specified criteria
count_kept_SNPs <- function(missing_rate, 
                            gds_obj, 
                            individuals,
                            minor_allele_freq) {
  my_SNPs <- snpgdsSelectSNP(gdsobj = gds_obj, 
                             autosome.only = FALSE, 
                             remove.monosnp = TRUE, 
                             maf = minor_allele_freq, 
                             missing.rate = missing_rate, 
                             sample.id = individuals,
                             verbose = TRUE)
  my_SNPIDs <- unique(unlist(my_SNPs))
  return(length(my_SNPIDs))
}


###########
# GLOBALS #
###########


GDSfile <- 'output/060_pop_genet/snps.gds'
SNPfilterrates <- seq(0,1, 0.05)
MAF <- 0.05

para_info <- read_delim('output/010_config/tidy_sample_info.tsv', delim = '\t')

########
# MAIN #
########

# Loads the data
gds_data <- snpgdsOpen(GDSfile)



names(SNPfilterrates) <- SNPfilterrates


# Lists all the filtered SNPs at each of the different allowed missing rates 
SNPlists <- lapply(SNPfilterrates, 
                   get_SNP_set, 
                   gds_obj = gds_data, 
                   minor_allele_freq = MAF)


# Simplifies the list (containing sublists) down into a vector so we can use it
number_SNPs <- sapply(SNPlists, length)


# Produces a tibble of the number of SNPs kept at each allowed missing rate
missing_rate_SNPs <- tibble(SNP_missing_rate = names(number_SNPs), 
                            number_SNPs)


# Plots the number of SNPs kept at each allowed missing rate, with a loess line
ggplot(data = missing_rate_SNPs, 
       mapping = aes(x = SNP_missing_rate, 
                     y = number_SNPs, group = 1)) + 
  geom_smooth(se = FALSE) + 
  geom_point() + 
  geom_vline(xintercept = 5) +
  scale_y_log10(breaks = trans_breaks("log10", 
                                      function(x) 10^x),
                labels = trans_format("log10", 
                                      math_format(10^.x))) +
  theme_classic() + 
  labs(x = "Allowed SNP missing rate", 
       y = "Number of retained SNPs")



# Gets the missing rate in each sample, and makes it into a tibble along with the SNP missing rates
sample_missing_rate_list <- lapply(SNPlists, 
                                   get_sample_missing_rates, 
                                   gds_obj = gds_data)

sample_missing_rate <- bind_rows(sample_missing_rate_list, 
                                 .id = "SNP_missing_rate")



# Format data for plotting. Add a column giving the value of the 80th quantile for samples at each missing rate.
plot_data <- sample_missing_rate %>% 
  group_by(individual) %>% 
  mutate("population" = substr(individual, 0, 1)) %>% 
  group_by(SNP_missing_rate) %>% 
  mutate(q80 = quantile(sample_missing_rate, 0.80)) %>% 
  mutate(population = str_replace_all(individual, "[[:digit:]]", ""))


# Add parasitism info to samples
indiv_status <-  para_info$Parasitism
names(indiv_status) <- para_info$Individual
plot_data$Parasitism <- indiv_status[plot_data$individual]


# Shows how many parasitized/unparasitized individuals there are in each population in unfiltered data
plot_data %>% 
  group_by(population, Parasitism) %>% 
  summarise(n = length(unique(individual)))

# Returns the number of individuals in the dataset prior to filtering
length(unique(sample_missing_rate$individual))


# Plots the SNP and sample missing rates, by population and parasitism status, for missing rates from 0 to 1
# has multiple bars
ggplot(data = plot_data, 
       mapping = aes(x = SNP_missing_rate, 
                     y = sample_missing_rate,
                     colour = population)) +
  geom_boxplot(aes(colour = population),
               outlier.size = -1, 
               position = position_dodge(width = 0.9)) +
  geom_point(aes(shape = Parasitism,
                 group = population),
             position = position_jitterdodge(jitter.width = 0.15, 
                                             dodge.width = 0.9),
             alpha = 0.4) +
  labs(x = "Allowed SNP missing rate", 
       y = "Sample missing rate") +
  theme_classic() +
  scale_color_discrete(name = "Population") +
  scale_shape_discrete(name = "Parasitism", 
                       labels = c("Parasitised", 
                                  "Unparasitised"))

# Sample missing rate at SNP_missing_rate == "0.2 
ggplot(data = subset(plot_data, SNP_missing_rate == "0.2"), 
       mapping = aes(x = population, 
                     y = sample_missing_rate,
                     colour = population)) +
  geom_boxplot(outlier.size = -1, 
               position = position_dodge(width = 0.9)) +
  geom_point(aes(shape = Parasitism,
                 group = population),
             position = position_jitterdodge(jitter.width = 0.45, 
                                             dodge.width = 0.9),
             alpha = 0.5,
             size = 2) +
  geom_errorbar(width = 0.95, 
                aes(ymax = q80, 
                    ymin = q80), 
                colour = "black") +
  labs(x = "Population", 
       y = "Sample missing rate") +
  theme_classic() +
  scale_color_discrete(name = "Population") +
  scale_shape_discrete(name = "Parasitism", 
                       labels = c("Parasitised", 
                                  "Unparasitised"))


