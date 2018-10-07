#!/usr/bin/env Rscript

library(grid)
library(gtable)
library(reshape2)
library(SNPRelate)
library(tidyverse)

###########
# Globals #
###########

GDSfile <- 'output/060_pop_genet/snps.gds'
sample_info <- 'output/010_config/tidy_sample_info.tsv'

########
# MAIN #
########

para_info <-  read_delim(sample_info, delim = '\t')
GDSdata <-  snpgdsOpen(GDSfile)

#Calculate Identity by state
IBS <- snpgdsIBS(gdsobj = GDSdata, autosome.only = FALSE)

#Reformat data
IBS_data <-  IBS$ibs
colnames(IBS_data) <- IBS$sample.id 
rownames(IBS_data) <- IBS$sample.id

# Define order of individuals by clustering (ward.D agglomeration method)
hc <- hclust(as.dist(1 - IBS_data), method = "ward.D")
indiv_order <- hc$labels[hc$order]


IBS_melt <- melt(IBS_data)

# Add in parasitism info
indiv_status <-  para_info$Parasitism
names(indiv_status) <- para_info$Individual
IBS_melt$para <- indiv_status[IBS_melt$Var1]


# Re-order the individuals based on clustering results
IBS_plot <- IBS_melt %>%
  mutate(Var1 = factor(Var1, levels = indiv_order), Var2 = factor(Var2, levels = indiv_order))

for_pop <- IBS_plot %>% 
  filter(Var2 == 'R2') %>% 
  mutate(pop = str_replace_all(Var1, "[[:digit:]]", ""))

for_para <- IBS_plot %>% 
  filter(Var2 == 'R2') 



# Define colours for plot components
hs <- RColorBrewer::brewer.pal(6, "YlOrRd")
pop_col <- RColorBrewer::brewer.pal(4, "PiYG")
names(pop_col) <- c("I", "L", "R", "Rpoa")
para_col <- c("#99c2ff", "#003cb3")
names(para_col) <-  c("non-parasitized", "parasitized")
# pop_col <-  c("#D01C8B", "#F1B6DA", "#B8E186")
# pop_col <- c("salmon1", "salmon3", "#B8E186", "#4DAC26")

# Generate panels for the plot
a <- ggplot(IBS_plot, aes(x = Var1, y = Var2, fill = value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = hs) +
  theme_void() +
  theme(legend.position = "none")

b <- ggplot(for_pop, aes(x = Var1, y = Var2, fill = pop)) + 
  geom_raster() +
  theme_void() +
  scale_fill_manual(values = pop_col) +
  theme(legend.position = "none",
        plot.margin=unit(c(0, 0, 0, 0), "cm"))

c <- ggplot(for_para, aes(x = Var1, y = Var2, fill = para)) + 
  geom_raster() +
  theme_void() +
  scale_fill_manual(values = para_col) +
  theme(legend.position = "none",
        plot.margin=unit(c(-0.2, 0, 0.2, 0), "cm"))

d <- ggplot(for_pop, aes(x = Var2, y = Var1, fill = pop)) + 
  geom_raster() +
  theme_void() +
  scale_fill_manual(values = pop_col) +
  theme(legend.position = "none", 
        plot.margin=unit(c(0, 0, 0, 0), "cm"))

e <- ggplot(for_para, aes(x = Var2, y = Var1, fill = para)) + 
  geom_raster() +
  theme_void() +
  scale_fill_manual(values = para_col) +
  theme(legend.position = "none", 
        plot.margin=unit(c(0, -0.2, 0, 0.2), "cm"))

# Generate legends for the plot (same as panels a, b and c, but with legend allowed)
f <- ggplot(IBS_plot, aes(x = Var1, y = Var2, fill = value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = hs) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) + 
  labs(fill = "Identity by state")

g <- ggplot(for_pop, aes(x = Var1, y = Var2, fill = pop)) + 
  geom_raster() +
  scale_fill_manual(values = pop_col, labels = c("Invermay", "Lincoln", "Ruakura", "Ruakura (Poa)  ")) +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(fill = "Population")

h <- ggplot(for_para, aes(x = Var1, y = Var2, fill = para)) + 
  geom_raster() +
  theme_void() +
  scale_fill_manual(values = para_col, labels = c("non-parasitised", "parasitised")) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(fill = "Parasitism")

# Convert plot components into grobs, create additional blank grob
panel_a <- ggplotGrob(a)
panel_b <- ggplotGrob(b)
panel_c <- ggplotGrob(c)
panel_d <- ggplotGrob(d)
panel_e <- ggplotGrob(e)
blank_one <- rectGrob(gp = gpar(fill = "white",col = "white"))
legend_ibs <- gtable_filter(ggplotGrob(f), 'guide-box')
legend_pop <- gtable_filter(ggplotGrob(g), 'guide-box')
legend_para <- gtable_filter(ggplotGrob(h), 'guide-box')

# Combine components and plot them
combined_legends <- gtable_matrix(name = "legends",
                                  grobs = matrix(list(blank_one,legend_ibs,legend_pop, legend_para, blank_one), ncol = 1),
                                  widths = unit(3,"cm"),
                                  heights = unit(c(7,4,4,2,8), "cm"))

plot_pos <- matrix(list(panel_e, panel_d, panel_a, combined_legends,
                        blank_one, blank_one, panel_b, blank_one,
                        blank_one, blank_one, panel_c, blank_one), byrow = TRUE, nrow = 3, ncol = 4) 


g <- gtable_matrix(name = "IBS_plot",
                   grobs = plot_pos, 
                   widths = unit(c(1,1,20,5), "cm"),
                   heights = unit(c(20, 1, 1), "cm"))

grid.newpage()
grid.draw(g)

