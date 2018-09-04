library(boa)
library(coda)
library(tidyverse)



sel <- read.table('output/070_pop_tests/compared_4pops.sel')
fst_tbl <- read.table('output/070_pop_tests/compared_4pops_fst.txt')

sel <- read.table('output/070_pop_tests/compared_para.sel')
fst_tbl <- read.table('output/070_pop_tests/compared_para_fst.txt')

sel <- read.table('output/070_pop_tests/compared_2pops.sel')
fst_tbl <- read.table('output/070_pop_tests/compared_2pops_fst.txt')


# Evaluate convergence
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
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(fst_tbl, aes(x = qval, y = alpha, colour = sig)) +
  geom_point(alpha = 0.3)

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








