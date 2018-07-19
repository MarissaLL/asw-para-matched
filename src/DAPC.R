#!/usr/bin/env Rscript

library(dplyr)
library(adegenet)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(tidyverse)

###########
# Globals #
###########

plink_file <- 'output/060_pop_genet/plink.raw'
para_file <-  'sample_catalog.csv'


########
# Main #
########

# Import the data 
SNP_data <- read.PLINK(plink_file)


##### NAs #####

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


##### Add parasitism info to the genlight #####

# Read in and format the population and parasitism data for all (unfiltered) individuals 
para <- read.csv(para_file) %>% 
  as_tibble() %>% 
  filter(Head =="TRUE") %>% 
  mutate_at("Parasitoid", replace_na, "FALSE") %>%       
  mutate("Parasitism" = if_else(Parasitoid=="TRUE", "parasitized","non-parasitized")) %>% 
  mutate("pop_para" = paste(Location,Parasitism, sep = "_")) %>% 
  mutate("pop_para_past" = paste(Location, Parasitism, Pasture, sep = "_"))

# Create a tibble of the filtered individuals (those in the genlight)
filtered_indivs <- indNames(no_NA) %>% 
  as.tibble() %>% 
  rename(Individual = value)

# Filter the population parasitism data to only include the filtered individuals
filtered_para <- semi_join(para, filtered_indivs)

# Add population and parasitism data into the genlight
indiv_to_pop <- filtered_para$pop_para_past
names(indiv_to_pop) <- filtered_para$Individual
pop(no_NA) <- indiv_to_pop[no_NA$ind.names]


##### PCA #####

# Do PCA on the no_NA data
pca2 <- glPca(no_NA, nf = 50)

# Format the PCA information to plot with ggplot
PCA_scores <- as.data.frame(pca2$scores)
PCA_scores ["Individual"] <- rownames(PCA_scores)
PCA_pop <- as_tibble(PCA_scores) 
PCA_plot <-  left_join(PCA_pop, filtered_para, by = "Individual")
  
# Generate convex hulls for each group. Make sure there are no NA values when doing this
df <- PCA_plot
df <- na.omit(df)

hulls <- df %>% 
  group_by(pop_para_past) %>% 
  slice(chull(PC1, PC2))

# Plot PC1 and PC2 with ggplot
ggplot(data = PCA_plot, aes(x = PC1, y = PC2, colour = Location, shape = Parasitism)) +
  geom_polygon(data = hulls, aes(fill = pop_para_past), alpha = 0.1, colour = NA) +
  geom_point() +
  theme_classic()

# Plot more PCs using ggplot facets, after converting data to the appropriate form
PCA_plot_facets <- PCA_plot %>% 
  select(-c(Pasture, Head, Body, Parasitoid, sg_box, box_comment)) %>% 
  gather(component,score, -c(Individual, pop_para_past, Location, pop_para, Parasitism))


######### what is y??
ggplot(data = PCA_plot_facets, 
       mapping = aes(x = Location,
                     y = 
                     colour = Parasitism)) +
  facet_wrap(~ component) +
  geom_point(position = position_jitter(0.2), 
             size = 1.5) + 
  theme_grey() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(x = "Population", 
       y = "Principal component")

############################ DAPC ############################### Set xval to have lower

# DAPC between 8 groups

# Change population in the genlight to contain parasitism and location information
indiv_to_pop <- filtered_para$pop_para_past
names(indiv_to_pop) <- filtered_para$Individual
no_NA_8 <- no_NA
pop(no_NA_8) <- indiv_to_pop[no_NA_8$ind.names]

# Check that the pop details match with the filenames (eg Karamea11p individual has pop of Karamea_parasitized)
data.frame(i = indNames(no_NA_8), 
           p = pop(no_NA_8))


dapc_unfit_8 <- dapc(no_NA_8, n.pca = 50, n.da = 10)
scatter(dapc_unfit_8)

temp <- optim.a.score(dapc_unfit_8, n.pca = 1:50, smart = FALSE, n.sim = 50)
temp$mean

x <- no_NA_8
mat <- tab(no_NA_8, NA.method="mean")
grp <- pop(x)
par(mar=c(4,4,4,4))
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval[2:6]
xval2 <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                  result = "overall", center = TRUE, scale = FALSE,
                  n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval2[2:6]

dapc_opt_8 <- dapc(no_NA, n.pca = 80, n.da = 10)

DAPC_scores_8 <- as.data.frame(dapc_opt_8$ind.coord)
DAPC_scores_8 ["Individual"] <- rownames(DAPC_scores_8)

DAPC_pop_8 <- as_tibble(DAPC_scores_8) %>% 
  left_join(filtered_para, by = "Individual") %>% 
  select(-c(Pasture, Head, Body, Parasitoid, sg_box, box_comment))

#
DAPC_plot_8 <- DAPC_pop_8 %>% 
  gather(component,score, -c(Individual, pop_para_past, Location, pop_para, Parasitism))

# Plot the first 4 linear discriminants
ggplot(data = DAPC_plot_8, 
       mapping = aes(x = Location, 
                     y = score, 
                     colour = Parasitism)) +
  facet_wrap(~ component) +
  geom_point(position = position_jitter(0.2), 
             size = 1.5) + 
  theme_grey() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(x = "Population", 
       y = "Linear discriminant score")


#################### DAPC between 2 groups
indiv_to_para <- filtered_para$Parasitism
names(indiv_to_para) <- filtered_para$Individual
no_NA_2 <- no_NA
pop(no_NA_2) <- indiv_to_para[no_NA_2$ind.names]

# Check that the pop details match with the filenames (eg Karamea11p individual has pop of Karamea_parasitized) - can't actually tell here
data.frame(i = indNames(no_NA_2), 
           p = pop(no_NA_2))

dapc_unfit_2 <- dapc(no_NA_2, n.pca = 50, n.da = 10)
scatter(dapc_unfit_2)

# Optim a-score
temp <- optim.a.score(dapc_unfit_2, n.pca = 1:50, smart = FALSE, n.sim = 50)
temp$mean

# Cross-validation
x_2 <- no_NA_2
mat_2 <- tab(no_NA_2, NA.method="mean")
grp_2 <- pop(x_2)
xval <- xvalDapc(mat_2, grp_2, n.pca.max = 10, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval[2:6]
xval2 <- xvalDapc(mat_2, grp_2, n.pca.max = 300, training.set = 0.9,
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

############################## DAPC of para, separated by population ###############################

Lincoln_only <- no_NA[str_detect(indNames(no_NA), "L")]

Ruakura_only <- no_NA[str_detect(indNames(no_NA), "R")]

Invermay_only <- no_NA[str_detect(indNames(no_NA), "I")]


# Change population in the genlight to contain parasitism and location information
indiv_to_para <- filtered_para$Parasitism
names(indiv_to_para) <- filtered_para$Individual
no_NA_L <- Lincoln_only
pop(no_NA_L) <- indiv_to_para[no_NA_L$ind.names]

# Check that the pop details match with the filenames (eg Karamea11p individual has pop of Karamea_parasitized)
data.frame(i = indNames(no_NA_L), 
           p = pop(no_NA_L))


dapc_opt_L <- dapc(no_NA_L, n.pca = 80, n.da = 10)

DAPC_scores_L <- as.data.frame(dapc_opt_L$ind.coord)
DAPC_scores_L ["Individual"] <- rownames(DAPC_scores_L)

DAPC_pop_L <- as_tibble(DAPC_scores_L) %>% 
  left_join(filtered_para, by = "Individual") %>% 
  select(-c(Pasture, Head, Body, Parasitoid, sg_box, box_comment))

# Figure out if there is a way to make this look better
scatter(dapc_opt_L)

ggplot(data = DAPC_pop_L, 
       mapping = aes(x = Location, 
                     y = LD1, 
                     colour = Parasitism)) +
  geom_point(position = position_jitter(0.1), 
             size = 1.5) + 
  theme_grey() +
  theme(panel.background = element_rect(fill = "white")) +
  labs(x = "Population", 
       y = "Linear discriminant score")

contrib <- loadingplot(dapc_opt_L$var.contr, axis=1,
                       lab.jitter=1, threshold = 0.005, xlab = "",
                       ylab = "", main = NULL)
