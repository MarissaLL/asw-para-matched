#!/usr/bin/env Rscript

# Note that this script also produces the cross-validation and a-score optimisation results. These are commented out now that the number of PCs to retain has been selected

library(adegenet)
library(tidyverse)

###########
# Globals #
###########

plink_file <- snakemake@input[["plink_file"]]
para_file <-  snakemake@input[["para_file"]]

pca_results <- snakemake@output[["pca_results"]]
dapc_para_results <- snakemake@output[["dapc_para_results"]]
dapc_pop_results <- snakemake@output[["dapc_pop_results"]]
dapc_invermay_results <- snakemake@output[["dapc_invermay_results"]]
dapc_ruakura_results <- snakemake@output[["dapc_ruakura_results"]]
dapc_rpoa_results <- snakemake@output[["dapc_rpoa_results"]]
dapc_lincoln_results <- snakemake@output[["dapc_lincoln_results"]]
dapc_island_results <- snakemake@output[["dapc_island_results"]]
log_file <- snakemake@log[[1]]

# dev
 # plink_file <- 'output/060_pop_genet/plink.raw'
 # para_file <-  'output/010_config/tidy_sample_info.tsv'
# 
# pca_results <- 'output/071_DAPC/pca_results.tsv'
# dapc_para_results <- 'output/071_DAPC/dapc_para_results.tsv'
# dapc_pop_results <- 'output/071_DAPC/dapc_pop_results.tsv'
# dapc_invermay_results <- 'output/071_DAPC/dapc_invermay_results.tsv'
# dapc_ruakura_results <- 'output/071_DAPC/dapc_ruakura_results.tsv'
# dapc_rpoa_results <- 'output/071_DAPC/dapc_rpoa_results.tsv'
# dapc_lincoln_results <- 'output/071_DAPC/dapc_lincoln_results.tsv'
# dapc_island_results <- 'output/071_DAPC/dapc_island_results.tsv'

########
# Main #
########
# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

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
indiv_to_pop <- filtered_para$pop
names(indiv_to_pop) <- filtered_para$Individual

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
# Add information about the parasitism status of each sample to the 'pop' slot of the genlight
indiv_to_para <- filtered_para$Parasitism
names(indiv_to_para) <- filtered_para$Individual
no_NA_para <- no_NA
pop(no_NA_para) <- indiv_to_para[no_NA_para$ind.names]

# Check that the pop details match with the filenames (though can't actually tell for this one)
data.frame(i = indNames(no_NA_para), 
           p = pop(no_NA_para))

dapc_unfit_para <- dapc(no_NA_para, n.pca = 50, n.da = 10)
# scatter(dapc_unfit_para)
# 
# 
# # Optim a-score
# temp <- optim.a.score(dapc_unfit_para, n.pca = 1:50, smart = FALSE, n.sim = 50)
# temp$mean
# 
# # Cross-validation
# mat_para <- tab(no_NA_para, NA.method="mean")
# grp_para <- pop(no_NA_para)
# xval <- xvalDapc(mat_para, grp_para, n.pca.max = 50, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
# xval[2:6]

# A-score optimisation and cross validation indicated that none of the discriminant functions perform better than chance, 
# and there isn't an optimal number of PCs to use, so just plot the DAPC that has already run (18, and 25 were the best)

# Format data to be able to plot it
DAPC_scores_para <- as.data.frame(dapc_unfit_para$ind.coord)
DAPC_scores_para ["Individual"] <- rownames(DAPC_scores_para)

DAPC_para <- as_tibble(DAPC_scores_para) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot parasitism DAPC without having to re-run DAPC
write_delim(DAPC_para, dapc_para_results, delim = '\t')



########## DAPC between 4 pops #################

# Add information about the population each sample came from to the 'pop' slot of the genlight
no_NA_pop <- no_NA
pop(no_NA_pop) <- indiv_to_pop[no_NA_pop$ind.names]

# # Check that the pop details match with the filenames (though can't actually tell for this one)
# data.frame(i = indNames(no_NA_pop), 
#            p = pop(no_NA_pop))
# 
# dapc_unfit_pop <- dapc(no_NA_pop, n.pca = 50, n.da = 10)
# scatter(dapc_unfit_pop)
# 
# 
# # Optim a-score
# temp <- optim.a.score(dapc_unfit_pop, n.pca = 1:50, smart = FALSE, n.sim = 50)
# temp$mean
# 
# # Cross-validation
# mat_pop <- tab(no_NA_pop, NA.method="mean")
# grp_pop <- pop(no_NA_pop)
# xval <- xvalDapc(mat_pop, grp_pop, n.pca.max = 50, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
# xval[2:6]

# Re-run the DAPC with the optimal number of PCs (10)
dapc_pop <- dapc(no_NA_pop, n.pca = 10, n.da = 10)

# Format data to be able to plot it
DAPC_scores_pop <- as.data.frame(dapc_pop$ind.coord)
DAPC_scores_pop ["Individual"] <- rownames(DAPC_scores_pop)

DAPC_pop <- as_tibble(DAPC_scores_pop) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot population DAPC without having to re-run DAPC
write_delim(DAPC_pop, dapc_pop_results, delim = '\t')


########### Para DAPC for individual pops ###########

# Separate the genlight by the four populations
pops <- seppop(no_NA_pop)

# Extract the four separate populations into their own genlights
invermay <- pops$I
ruakura <- pops$R
rpoa <- pops$Rpoa
lincoln <- pops$L

############## Invermay #############

# Change the pop slot to contain parasitism status
indiv_to_para_i <- filtered_para$Parasitism[filtered_para$Location == "Invermay"]
names(indiv_to_para_i) <- filtered_para$Individual[filtered_para$Location == "Invermay"]
pop(invermay) <- indiv_to_para_i[invermay$ind.names]

# Run a DAPC
# dapc_unfit_invermay <- dapc(invermay, n.pca = 50, n.da = 10)
# scatter(dapc_unfit_invermay)
# 
# # Optim a-score
# temp <- optim.a.score(dapc_unfit_invermay, n.pca = 1:50, smart = FALSE, n.sim = 50)
# temp$mean
# 
# # Cross-validation
# mat_i <- tab(invermay, NA.method="mean")
# grp_i <- pop(invermay)
# xval <- xvalDapc(mat_i, grp_i, n.pca.max = 50, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
# xval[2:6]

# Re-run the DAPC with the optimal number of PCs (4)
dapc_invermay <- dapc(invermay, n.pca = 4, n.da = 10)

# Format data to be able to plot it
DAPC_scores_i <- as.data.frame(dapc_invermay$ind.coord)
DAPC_scores_i ["Individual"] <- rownames(DAPC_scores_i)

DAPC_invermay <- as_tibble(DAPC_scores_i) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot population DAPC without having to re-run DAPC
write_delim(DAPC_invermay, dapc_invermay_results, delim = '\t')


############## Ruakura #############

# Change the pop slot to contain parasitism status
indiv_to_para_r <- filtered_para$Parasitism[filtered_para$pop == "R"]
names(indiv_to_para_r) <- filtered_para$Individual[filtered_para$pop == "R"]
pop(ruakura) <- indiv_to_para_r[ruakura$ind.names]

# Run a DAPC
# dapc_unfit_ruakura <- dapc(ruakura, n.pca = 50, n.da = 10)
# scatter(dapc_unfit_ruakura)
# 
# # Optim a-score
# temp <- optim.a.score(dapc_unfit_ruakura, n.pca = 1:50, smart = FALSE, n.sim = 50)
# temp$mean
# 
# # Cross-validation
# mat_r <- tab(ruakura, NA.method="mean")
# grp_r <- pop(ruakura)
# xval <- xvalDapc(mat_r, grp_r, n.pca.max = 50, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
# xval[2:6]

# Re-run the DAPC with the optimal number of PCs (6)
dapc_ruakura <- dapc(ruakura, n.pca = 6, n.da = 10)

# Format data to be able to plot it
DAPC_scores_r <- as.data.frame(dapc_ruakura$ind.coord)
DAPC_scores_r ["Individual"] <- rownames(DAPC_scores_r)

DAPC_ruakura <- as_tibble(DAPC_scores_r) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot population DAPC without having to re-run DAPC
write_delim(DAPC_ruakura, dapc_ruakura_results, delim = '\t')


############## Rpoa #############

# Change the pop slot to contain parasitism status
indiv_to_para_rpoa <- filtered_para$Parasitism[filtered_para$pop == "Rpoa"]
names(indiv_to_para_rpoa) <- filtered_para$Individual[filtered_para$pop == "Rpoa"]
pop(rpoa) <- indiv_to_para_rpoa[rpoa$ind.names]

# # Run a DAPC
# dapc_unfit_rpoa <- dapc(rpoa, n.pca = 50, n.da = 10)
# scatter(dapc_unfit_rpoa)
# 
# # Optim a-score
# temp <- optim.a.score(dapc_unfit_rpoa, n.pca = 1:50, smart = FALSE, n.sim = 50)
# temp$mean
# 
# # Cross-validation
# mat_rpoa <- tab(rpoa, NA.method="mean")
# grp_rpoa <- pop(rpoa)
# xval <- xvalDapc(mat_rpoa, grp_rpoa, n.pca.max = 50, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
# xval[2:6]

# Re-run the DAPC with the optimal number of PCs (6)
dapc_rpoa <- dapc(rpoa, n.pca = 6, n.da = 10)

# Format data to be able to plot it
DAPC_scores_rpoa <- as.data.frame(dapc_rpoa$ind.coord)
DAPC_scores_rpoa ["Individual"] <- rownames(DAPC_scores_rpoa)

DAPC_rpoa <- as_tibble(DAPC_scores_rpoa) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot population DAPC without having to re-run DAPC
write_delim(DAPC_rpoa, dapc_rpoa_results, delim = '\t')


############## Lincoln #############

# Change the pop slot to contain parasitism status
indiv_to_para_l <- filtered_para$Parasitism[filtered_para$Location =="Lincoln"]
names(indiv_to_para_l) <- filtered_para$Individual[filtered_para$Location == "Lincoln"]
pop(lincoln) <- indiv_to_para_l[lincoln$ind.names]

# # Run a DAPC
# dapc_unfit_l <- dapc(lincoln, n.pca = 50, n.da = 10)
# scatter(dapc_unfit_l)
# 
# # Optim a-score
# temp <- optim.a.score(dapc_unfit_l, n.pca = 1:50, smart = FALSE, n.sim = 50)
# temp$mean
# 
# # Cross-validation
# mat_l <- tab(lincoln, NA.method="mean")
# grp_l <- pop(lincoln)
# xval <- xvalDapc(mat_l, grp_l, n.pca.max = 50, training.set = 0.9,
#                  result = "groupMean", center = TRUE, scale = FALSE,
#                  n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
# xval[2:6]

# Re-run the DAPC with the optimal number of PCs (33)
dapc_l <- dapc(lincoln, n.pca = 33, n.da = 10)

# Format data to be able to plot it
DAPC_scores_l <- as.data.frame(dapc_l$ind.coord)
DAPC_scores_l ["Individual"] <- rownames(DAPC_scores_l)

DAPC_lincoln <- as_tibble(DAPC_scores_l) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot population DAPC without having to re-run DAPC
write_delim(DAPC_lincoln, dapc_lincoln_results, delim = '\t')

############ By Island ##############

# Change pop slot to contain island info
filtered_para <-  filtered_para %>% 
  mutate(island = case_when(Location == "Invermay" | Location == "Lincoln" ~ "South",
                                     Location == "Ruakura" ~ "North"))

indiv_to_island <- filtered_para$island 
names(indiv_to_island) <- filtered_para$Individual
no_NA_island <- no_NA
pop(no_NA_island) <- indiv_to_island[no_NA_island$ind.names]

# # Run a DAPC
dapc_unfit_is <- dapc(no_NA_island, n.pca = 50, n.da = 10)
scatter(dapc_unfit_is)
# 


# # Optim a-score
temp <- optim.a.score(dapc_unfit_is, n.pca = 1:50, smart = FALSE, n.sim = 50)
temp$mean
# 
# # Cross-validation
 mat_is <- tab(no_NA_island, NA.method="mean")
 grp_is <- pop(no_NA_island)
xval <- xvalDapc(mat_is, grp_is, n.pca.max = 50, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = 1:50, n.rep = 50, xval.plot = TRUE)
xval[2:6]

# Re-run the DAPC with the optimal number of PCs (1)
dapc_island <- dapc(no_NA_island, n.pca = 1, n.da = 10)

# Format data to be able to plot it
DAPC_scores_is <- as.data.frame(dapc_island$ind.coord)
DAPC_scores_is ["Individual"] <- rownames(DAPC_scores_is)

DAPC_island <- as_tibble(DAPC_scores_is) %>% 
  left_join(filtered_para, by = "Individual") 

# Write out data to plot population DAPC without having to re-run DAPC
write_delim(DAPC_island, dapc_island_results , delim = '\t')


## Calculate and plot SNP loadings for LD 1 ## 
contrib <- loadingplot(dapc_island$var.contr, axis=1, xlab = "SNP number",
                       ylab = "SNP Loading", lab.jitter=1, threshold = 0.0005, main = NULL)



# Log session info
sessionInfo()
