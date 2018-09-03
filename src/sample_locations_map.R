library(ggplot2)
library(ggmap)
library(maps)

# Load in map from maps package
nz_map <- map_data("nz")

# Define sampling location points. 
sampling_loc <-  data.frame(
  long = c(175.310, 172.472, 170.348),
  lat = c(-37.774,-43.634, -45.878),
  names = c('Ruakura', 'Lincoln', 'Invermay'),
  stringsAsFactors = FALSE)  

# Manually specify label position to make it better
lab_pos <- c(1.6, -0.8, -0.5)

# Plot the map
# axis lines on top and right currently not working
ggplot() + 
  geom_polygon(data = nz_map, 
               mapping = aes(x=long, 
                             y = lat, 
                             group = group), 
               fill = NA, colour = 'black' ) + 
  coord_fixed(1.4) + 
  geom_point(data = sampling_loc, 
             mapping = aes(x = long, 
                           y = lat), 
             colour = 'red', 
             size = 3.5) +
  geom_text(data = sampling_loc, 
            mapping = aes(x = long, 
                          y = lat, 
                          label = names), 
            hjust= lab_pos) +
  theme_classic() + 
  theme(axis.line.x.top = element_line(), 
        axis.line.y.right = element_line(), 
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 )) +
  labs(x = 'Longitude', 
       y = 'Latitude')
