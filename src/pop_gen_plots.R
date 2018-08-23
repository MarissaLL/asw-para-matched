library(ggplot2)
library(grid)
library(gtable)
library(SNPRelate)
library(tidyverse)
library(reshape2)
library(stringr)

# Fst
LR_fst <- read_delim('output/060_pop_genet/populations.fst_R-L.tsv', delim = '\t')
LR_phi <- read_delim('output/060_pop_genet/populations.phistats_R-L.tsv', delim = '\t')

ggplot(LR_fst, aes(x = `Overall Pi`,
                   y = `Corrected AMOVA Fst`)) +
  geom_point() 

sig <- filter(LR_fst, LR_fst$`Corrected AMOVA Fst`<0.005)

# Between populations summary Fst
Fst_summary <- read_delim('output/060_pop_genet/populations.fst_summary.tsv', delim = '\t', col_names = TRUE) 
pop_order <- c("I", "L", "R", "Rpoa")

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

ggplot(Fst_plot, aes(x = X1, y = variable, fill = value)) +
  geom_raster() + 
  theme_classic() +
  theme(axis.line = element_blank(), axis.ticks = element_blank())+
  labs(x = "Population", y = "Population")


# Private alleles
pop_sumstats <- read_delim('output/060_pop_genet/populations.sumstats.tsv', delim = '\t', skip = 4)

private_alleles <- pop_sumstats %>% 
  group_by(`Pop ID`) %>% 
  summarise(sum(Private == 1)) 

sum(pop_sumstats$Private == 1)
sum(pop_sumstats$Private == 0)

# Overall summary stats for populations

var_pos_pop_sumstats <-  read_delim('output/060_pop_genet/populations.sumstats_summary.tsv', delim = '\t', skip = 1, n_max = 4)
all_pos_pop_sumstats <- read_delim('output/060_pop_genet/populations.sumstats_summary.tsv', delim = '\t', skip = 7)


# Summary stat table formatted nicely for latex 
var_pos_pop_sumstats %>% 
  select(`# Pop ID`, Private, P, Obs_Het, Exp_Het, Pi, Fis ) %>% 
  xtable(digits = 5)





########## identity by state. ###############

GDSfile <- 'output/060_pop_genet/filtered_snps.gds'
GDSdata <-  snpgdsOpen(GDSfile)

IBS <- snpgdsIBS(gdsobj = GDSdata, autosome.only = FALSE)

IBS_plot <-  IBS$ibs
colnames(IBS_plot) <- IBS$sample.id 
rownames(IBS_plot) <- IBS$sample.id

# use ward.d2 or average
hc <- hclust(as.dist(1 - IBS_plot), method = "ward.D2")

indiv_order <- hc$labels[hc$order]



IBS_plot2 <- melt(IBS_plot) %>%
  mutate(Var1 = factor(Var1, levels = indiv_order), Var2 = factor(Var2, levels = indiv_order))

for_pop <- IBS_plot2 %>% 
  filter(Var2 == 'R2') %>% 
  mutate(pop = str_replace_all(Var1, "[[:digit:]]", ""))

para_info <-  read_delim('output/010_config/tidy_sample_info.tsv', delim = '\t')

for_para <- IBS_plot2 %>% 
  filter(Var2 == 'R2') 

indiv_status <-  para_info$Parasitism
names(indiv_status) <- para_info$Individual
for_para$para <- indiv_status[for_para$Var1]

# Colours for plot components. FIGURE OUT HOW TO APPLY THESE
hs <- RColorBrewer::brewer.pal(6, "YlOrRd")
pop_col <- RColorBrewer::brewer.pal(4, "PiYG")


# Order of plot margin units is c(top, right, bottom, left)

a <- ggplot(IBS_plot2, aes(x = Var1, y = Var2, fill = value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = hs) +
  theme_void() +
  theme(legend.position = "none", 
        plot.margin=unit(c(0,0,0,0), "cm"))

b <- ggplot(for_pop, aes(x = Var1, y = Var2, fill = pop)) + 
  geom_raster() +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))

c <- ggplot(for_para, aes(x = Var1, y = Var2, fill = para)) + 
  geom_raster() +
  theme_void() +
  theme(legend.position = "none",
        plot.margin=unit(c(-0.2, 0, 0.2, 0), "cm"))

d <- ggplot(for_pop, aes(x = Var2, y = Var1, fill = pop)) + 
  geom_raster() +
  theme_void() +
  theme(legend.position = "none", 
        plot.margin=unit(c(0, 0, 0, 0), "cm"))

e <- ggplot(for_para, aes(x = Var2, y = Var1, fill = para)) + 
  geom_raster() +
  theme_void() +
  theme(legend.position = "none", 
        plot.margin=unit(c(0, -0.2, 0, 0.2), "cm"))

grob_a <- ggplotGrob(a)
grob_b <- ggplotGrob(b)
grob_c <- ggplotGrob(c)
grob_d <- ggplotGrob(d)
grob_e <- ggplotGrob(e)
blank_one <- rectGrob(gp = gpar(fill = "white",col = "white"))

plot_pos <- matrix(list(grob_e, grob_d, grob_a, 
                        blank_one, blank_one, grob_b, 
                        blank_one, blank_one, grob_c), byrow = TRUE, nrow = 3, ncol = 3) 


g <- gtable_matrix(name = "IBS_plot",
                   grobs = plot_pos, 
                   widths = unit(c(1,1,25), "cm"),
                   heights = unit(c(25, 1, 1), "cm"))

grid.newpage()
grid.draw(g)








######################

# Plot phi against position. Not informative?
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

# Plot Fis against position. Not sure if this is informative either
ggplot(pop_sumstats, aes(x = `Chr`,
                         y = `Fis`)) +
  geom_point(position = position_jitter(0.48))

