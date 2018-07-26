library(boa)
library(coda)
library(tidyverse)


fst_tbl <- read.table('output/070_pop_tests/compared_pops_fst.txt')
sel <- read.table('output/070_pop_tests/compared_pops.sel', skip = 1)
AccRte <- read_table('output/070_pop_tests/compared_para_AccRte.txt')
props <-  read_delim('output/070_pop_tests/compared_para_prop.txt', delim = ' ', trim_ws = TRUE, skip = 4,  col_names = FALSE)
verif <- read.table('output/070_pop_tests/compared_para_Verif.txt', skip = 15)


# Evaluate convergence
chain <-  mcmc(sel, thin = 10)

plot(chain)
summary(chain)
autocorr.diag(chain) # What is too high?
effectiveSize(chain)

# Formal tests of convergence
geweke.diag(chain, frac1=0.1, frac2=0.5)
heidel.diag(chain, eps=0.1, pvalue=0.05)

# Run multiple chains to do this bit. Diagnostic based on multiple chains and comparing within and between chain variances
combined = mcmc.list(chain1,chain2)
plot(combined)
gelman.diag(combined)
gelman.plot(combined,ask)

# Fst
# posterior prob for model including sel
# log(posterior odds), including selection
# q-value including selection - is the p-value that has been FDR-adjusted
# alpha coefficient indicating strength and dir of selection - positive is diversifying selection
# Fst coefficient from posterior mean using model averaging


# Find PO threshold with FDR of no more than 5%

fst_tbl <- fst_tbl %>% 
  mutate(index = rownames(fst_tbl))

ggplot(fst_tbl, aes(x = qval, y = alpha)) +
  geom_point()

ggplot(fst_tbl, aes(x = index, y = alpha)) +
  geom_point() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Sel - the output from MCMC
# loglikelihood 
# Fst from each population

# Fis coefficient only when using dominant or AFLP data
# Need to turn on -all_trace to get alpha data - SEEMS LIKE I PROBABLY WANT THIS


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


# probably not valid way to plot
ggplot(sel, aes(x = logL, y = Fst1, colour = 'blue')) +
  geom_point()+
  geom_point(aes(y = Fst2, colour = 'green'))



# Verif
# This helps to verify that the input file was read correctly

# Do I actually want to plot any of this?
a1 <- verif$V6/(verif$V6+verif$V7)

ggplot(verif, aes(x = a1)) +
  geom_histogram(binwidth = 0.005)


# AccRte and prop 
# Optional files
# prop gives the results of the pilot runs
# AccRte gives evolution of the acceptance rate for the different model parameters
# Doesn't give much detail about them - leave this for now

ggplot(AccRte, aes(x = beta, y = ances)) +
  geom_point()

ggplot(props, aes(x = X7, y = X9)) +
  geom_point()

props <- props %>% 
  mutate(index = rownames(props))

ggplot(props, aes(x = index, y = X9)) + 
  geom_point()

# Unpick how it is doing this
plot_bayescan('output/070_pop_tests/compared_para_fst.txt')



