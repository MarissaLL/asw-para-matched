#!/usr/bin/env Rscript

library(adegenet)
library(tidyverse)

###########
# Globals #
###########

plink_file <- 'output/060_pop_genet/plink.raw'
para_file <-  'output/010_config/tidy_sample_info.tsv'

pca_results <- 'output/070_pop_tests/pca_results.tsv'
dapc_8_results <- 'output/070_pop_tests/dapc_8_results.tsv'

########
# Main #
########

# Import the data 
SNP_data <- read.PLINK(plink_file)

# Check the number of NAs in the data compared to the overall amount of data
sum(is.na(tab(SNP_data, 
              NA.method = "asis")))
length(tab(SNP_data)) 

# Replace NAs with the mean allele frequency. Convert the output back into a geenlight. Add lost info (ploidy) back in
NA_means <- tab(SNP_data, 
                NA.method = "mean")
no_NA <- new("genlight", 
             NA_means)
ploidy(no_NA) <- 2

# Check there are no NAs remaining
sum(is.na(tab(no_NA, NA.method = "asis")))
length(tab(no_NA))

# Read in population and parasitism data for all (unfiltered) individuals 
para <- read_delim(para_file, delim = '\t') 

# Create a tibble of the filtered individuals (those in the genlight)
filtered_indivs <- indNames(no_NA) %>% 
  as.tibble() %>% 
  rename(Individual = value)

# Filter the population parasitism data to only include the filtered individuals
filtered_para <- semi_join(para, filtered_indivs) %>% 
  mutate(pop = str_replace_all(Individual, "[[:digit:]]", ""))


##### PCA #####

# Add population and parasitism data into the genlight

pop(no_NA) <- indiv_to_pop[no_NA$ind.names]

# Do PCA on the no_NA data. Only use 9 PCs so they can be plotted reasonably 
pca <- glPca(no_NA, nf = 9)

# Format the PCA information to plot with ggplot
PCA_scores <- as.data.frame(pca$scores)
PCA_scores ["Individual"] <- rownames(PCA_scores)
PCA_pop <- as_tibble(PCA_scores) 
PCA_plot <-  left_join(PCA_pop, filtered_para, by = "Individual")

# Write out PCA results
write_delim(PCA_plot, pca_results, delim = '\t')

### This is then plot.


##### DAPC between parasitism status #####
indiv_to_para <- filtered_para$Parasitism
names(indiv_to_para) <- filtered_para$Individual
no_NA_para <- no_NA
pop(no_NA_para) <- indiv_to_para[no_NA_para$ind.names]

# Check that the pop details match with the filenames (eg Rpoa10 individual is from Ruakura Poa grass)
data.frame(i = indNames(no_NA_para), 
           p = pop(no_NA_para))

dapc_unfit_para <- dapc(no_NA_para, n.pca = 50, n.da = 10)
scatter(dapc_unfit_para)

# Optim a-score
temp <- optim.a.score(dapc_unfit_para, n.pca = 1:50, smart = FALSE, n.sim = 50)
temp$mean

# Cross-validation
mat_para <- tab(no_NA_para, NA.method="mean")
grp_para <- pop(no_NA_para)
xval <- xvalDapc(mat_para, grp_para, n.pca.max = 10, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval[2:6]
xval2 <- xvalDapc(mat_para, grp_para, n.pca.max = 300, training.set = 0.9,
                  result = "overall", center = TRUE, scale = FALSE,
                  n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval2[2:6]


dapc_opt_2 <- dapc(no_NA_2, n.pca = 80, n.da = 10)

DAPC_scores_2 <- as.data.frame(dapc_opt_2$ind.coord)
DAPC_scores_2 ["Individual"] <- rownames(DAPC_scores_2)

DAPC_pop_2 <- as_tibble(DAPC_scores_2) %>% 
  left_join(filtered_para, by = "Individual") %>% 
  select(-c(Pasture, Head, Body, Parasitoid, sg_box, box_comment))

ggplot(data = DAPC_pop_2, 
       mapping = aes(x = Location, 
                     y = LD1, 
                     colour = Parasitism)) +
  geom_point(position = position_jitter(0.2), 
             size = 1.5) + 
  theme_grey() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(x = "Population", 
       y = "Linear discriminant score")

######### Extract SNPs driving LD 1
contrib <- loadingplot(dapc_opt_2$var.contr, axis=1,
                       lab.jitter=1, threshold = 0.005, xlab = "",
                       ylab = "", main = NULL)

########## DAPC between pops #################

# Set up a named THING??
indiv_to_pop <- filtered_para$pop
names(indiv_to_pop) <- filtered_para$Individual
no_NA_pop <- no_NA
pop(no_NA_pop) <- indiv_to_pop[no_NA_pop$ind.names]


dapc_unfit_pop <- dapc(no_NA_pop, n.pca = 50, n.da = 10)
scatter(dapc_unfit_pop)
