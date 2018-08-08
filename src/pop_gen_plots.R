library(ggplot2)
library(SNPRelate)
library(tidyverse)
library(reshape2)
library(stringr)

# Fst
LR_fst <- read_delim('output/060_pop_genet/populations.fst_R-L.tsv', delim = '\t')
LR_phi <- read_delim('output/060_pop_genet/populations')

ggplot(LR_fst, aes(x = `Overall Pi`,
                   y = `Corrected AMOVA Fst`)) +
  geom_point() +

sig <- filter(LR_fst, LR_fst$`Corrected AMOVA Fst`<0.005)

# Between populations summary Fst
Fst_summary <- read_delim('output/060_pop_genet/populations.fst_summary.tsv', delim = '\t', col_names = TRUE) 


 Fst_table <- Fst_summary %>% 
  melt(id.vars = 'X1') %>% 
   mutate(value = if_else(X1 == variable, 0, as.numeric(value))) %>% 
   rbind(tibble(X1 = "I", variable = "I", value = 0)) %>% 
   filter(!is.na(value))


Fst_plot <- Fst_table %>% 
  mutate(X2 = variable, X3 = X1) %>% 
  transmute(X1 = X2, variable = X3, value) %>% 
  rbind(Fst_table) %>% 
  distinct() %>% 
  mutate(X1 = factor(X1, levels = pop_order), variable = factor(variable, levels = pop_order))

pop_order <- c("I", "L", "R", "Rpoa")


ggplot(Fst_plot, aes(x = X1, y = variable, fill = value)) +
  geom_raster() + 
  theme_classic() +
  theme(axis.line = element_blank(), axis.ticks = element_blank())+
  labs(x = "Population", y = "Population")







# Private alleles
pop_sumstats <- read_delim('output/060_pop_genet/populations.sumstats.tsv', delim = '\t', skip = 4)


ggplot(pop_sumstats, aes(x = `Chr`,
                         y = `Fis`)) +
  geom_point(position = position_jitter(0.48))



private_alleles <- pop_sumstats %>% 
  group_by(`Pop ID`) %>% 
  summarise(sum(Private == 1)) 


sum(pop_sumstats$Private == 1)
sum(pop_sumstats$Private == 0)



# Phi
pop_phistats <- read_delim('output/060_pop_genet/populations.phistats.tsv', delim = '\t', skip = 5)

ggplot(pop_phistats, aes(x = `# Locus ID`,
              y = `Fst'`)) +
  geom_point()


ggplot(pop_phistats, aes(x = `Chr`,
              y = `Fst'`)) +
  geom_point(position = position_jitter(0.48))



ggplot(pop_phistats, aes(x = `Chr`,
              y = `phi_st`)) +
  geom_point(position = position_jitter(0.48))



# identity by state

GDSfile <- 'abc_test/snps.gds'
GDSdata <-  snpgdsOpen(GDSfile)

IBS <- snpgdsIBS(gdsobj = GDSdata, autosome.only = FALSE)

IBS_plot <-  IBS$ibs
colnames(IBS_plot) <- IBS$sample.id 
rownames(IBS_plot) <- IBS$sample.id


hc <- hclust(as.dist(1 - IBS_plot), method = "complete")


indiv_order <- hc$labels[hc$order]

hs <- RColorBrewer::brewer.pal(6, "YlOrRd")

IBS_plot2 <- melt(IBS_plot) %>%
  mutate(Var1 = factor(Var1, levels = indiv_order), Var2 = factor(Var2, levels = indiv_order))
  
ggplot(IBS_plot2, aes(x = Var1, y = Var2, fill = value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = hs) +
  theme_classic() +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Individual", y = "Individual")


# Overall summary stats for populations

var_pos_pop_sumstats <-  read_delim('output/060_pop_genet/populations.sumstats_summary.tsv', delim = '\t', skip = 1, n_max = 4)
all_pos_pop_sumstats <- read_delim('output/060_pop_genet/populations.sumstats_summary.tsv', delim = '\t', skip = 7)


# format nicely for latex 
var_pos_pop_sumstats %>% 
  select(`# Pop ID`, Private, P, Obs_Het, Exp_Het, Pi, Fis ) %>% 
  xtable(digits = 5)



