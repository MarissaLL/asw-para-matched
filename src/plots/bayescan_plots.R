#!/usr/bin/env Rscript

library(boa)
library(coda)
library(scales)
library(tidyverse)

#############
# FUNCTIONS #
#############

list_sel_data <- function(data_dir){
  key_files <- list.files(data_dir,
                          recursive = FALSE,
                          pattern = ".sel",
                          full.names = TRUE) 
  sel_file_names <- key_files %>% str_replace("//", "/")
  return(sel_file_names)}

list_fst_data <- function(data_dir){
  key_files <- list.files(data_dir,
                          recursive = FALSE,
                          pattern = "_fst.txt",
                          full.names = TRUE) 
  fst_file_names <- key_files %>% str_replace("//", "/")
  return(fst_file_names)}

###########
# Globals #
###########

sel <- read.table('output/070_bayescan/compared_4pops.sel')
sel <- read.table('output/070_bayescan/compared_lincoln.sel')
sel <- read.table('output/070_bayescan/compared_ruakura_poa.sel')
sel <- read.table('output/070_bayescan/compared_para.sel')
sel <- read.table('output/070_bayescan/compared_2pops.sel')
sel <- read.table('output/070_bayescan/compared_ruakura.sel') # This one did not pass the Geweke diagnostic test
sel <- read.table('output/070_bayescan/compared_invermay.sel')
# sel <- read.table('output/070_bayescan/compared_2pops.sel')

fst_tbl <- read.table('output/070_bayescan/compared_4pops_fst.txt')
fst_tbl_l <- read.table('output/070_bayescan/compared_lincoln_fst.txt')
fst_tbl_rpoa <- read.table('output/070_bayescan/compared_ruakura_poa_fst.txt')
fst_tbl <- read.table('output/070_bayescan/compared_para_fst.txt')
fst_tbl <- read.table('output/070_bayescan/compared_2pops_fst.txt')
fst_tbl_r <- read.table('output/070_bayescan/compared_ruakura_fst.txt')
fst_tbl_i <- read.table('output/070_bayescan/compared_invermay_fst.txt')
# sel <- read.table('output/070_bayescan/compared_2pops_fst.txt')


# 

data_dir <- 'output/070_bayescan/'

########
# Main #
########

################## Trying to set up to analyse all the results in one go ###############
# List data files
sel_files <-  list_sel_data(data_dir)
fst_files <- list_fst_data(data_dir)

# Read in data files
sels <- lapply(sel_files, read.table)

# Evaluate convergence
chains <- lapply(sel_data, mcmc, thin = 10)


######################### One by one
chain <-  mcmc(sel, thin = 10)

plot(chain)
summary(chain)
autocorr.diag(chain)
effectiveSize(chain)

# Formal tests of convergence
geweke.diag(chain, frac1=0.1, frac2=0.5)
heidel.diag(chain, eps=0.1, pvalue=0.05)


# Plot alpha coefficients 
fst_tbl <- fst_tbl %>% 
  mutate(index = rownames(fst_tbl)) %>% 
  mutate(sig = if_else(qval <= 0.05, 'sig', 'non-sig'))

ggplot(fst_tbl, aes(x = index, y = alpha, colour = sig)) +
  geom_point(alpha = 0.3) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none") 
  

ggplot(fst_tbl, aes(x = qval, y = alpha, colour = sig)) +
  geom_point(alpha = 0.3)

sum(fst_tbl$sig == 'sig')
sum(fst_tbl$sig == 'non-sig')




##### Other extra stuff ##############
# Plot trace for each of the parameters. MAKE THIS LINES
sel <- sel %>% 
  mutate(index = rownames(sel))

ggplot(sel, aes(x = index, y = logL)) +
  geom_point()

ggplot(sel, aes(x = index, y = Fst1)) +
  geom_point()

ggplot(sel, aes(x = index, y = Fst2)) +
  geom_point()


# Plot posterior distribution of parameters of interest (e.g. Fst of pop1), include 95% Highest Probability Density Interval (HDPI)
hbd <- boa.hpd(sel[['Fst1']], 0.05)

ggplot(sel, aes(x = Fst1)) +
  geom_density() +
  geom_vline(aes(xintercept = hbd[1]))+
  geom_vline(aes(xintercept = hbd[2]))



##### Four pops together #############

lincoln_points <- fst_tbl_l %>% 
  mutate(index = rownames(fst_tbl_l)) %>% 
  mutate(sig = if_else(qval <= 0.05, 'sig', 'non-sig')) %>% 
  mutate(location = "Lincoln")

rpoa_points <- fst_tbl_rpoa %>% 
  mutate(index = rownames(fst_tbl_rpoa)) %>% 
  mutate(sig = if_else(qval <= 0.05, 'sig', 'non-sig')) %>% 
  mutate(location = "Ruakura (Poa)")

ruakura_points <- fst_tbl_r %>% 
  mutate(index = rownames(fst_tbl_r)) %>% 
  mutate(sig = if_else(qval <= 0.05, 'sig', 'non-sig')) %>% 
  mutate(location = "Ruakura")

invermay_points <- fst_tbl_i %>% 
  mutate(index = rownames(fst_tbl_i)) %>% 
  mutate(sig = if_else(qval <= 0.05, 'sig', 'non-sig')) %>% 
  mutate(location = "Invermay")

  
combined_results <- bind_rows(lincoln_points, rpoa_points, ruakura_points, invermay_points)

ggplot(combined_results, aes(x = index, y = alpha, colour = sig)) +
  geom_point(alpha = 0.3) +
  facet_grid(~location) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  scale_y_continuous(labels =function(n){format(n, scientific = FALSE)})


