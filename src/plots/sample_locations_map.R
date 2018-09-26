library(ggmap)
library(maps)
library(tidyverse)
#loc_file <- snakemake@input[[1]]
loc_file <- 'sampling_locations.txt'

# Load file containing the coordinates for sampling locations
sampling_loc <- read_delim(file = loc_file , delim = '\t', col_names = FALSE) %>% 
  rename(pop = X1, lat = X2, long = X3) 

# Load in map from maps package
nz_map <- map_data("nz")

# Manually specify label position to make it better
lab_pos <- c(1.3, -0.5, -0.3)
#loc_col <-  c("#B8E186", "#F1B6DA","#D01C8B") # Colours used on poster


# Plot the map
ggplot() + 
  geom_polygon(data = nz_map, 
               mapping = aes(x=long, 
                             y = lat, 
                             group = group), 
               fill = NA, colour = 'black' ) + 
  coord_fixed(1.4) + 
  geom_point(data = sampling_loc, 
             mapping = aes(x = long, 
                           y = lat), colour = "red", 
             size = 3.5) +
  geom_text(data = sampling_loc, 
            mapping = aes(x = long, 
                          y = lat, 
                          label = pop), 
            hjust= lab_pos,
            size = 5) +
  theme_classic() + 
  theme(axis.line.x.top = element_line(), 
        axis.line.y.right = element_line(), 
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.title = element_text(size = 18, margin = margin(r = 12)),
        axis.text = element_text(size = 12)) +
  labs(x = 'Longitude', 
       y = 'Latitude')
