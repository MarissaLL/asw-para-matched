library(ggplot2)
library(ggmap)
library(maps)

# Load in map from maps package
nz_map <- map_data("nz")

# Define sampling location points. These are currently based on locations of Hamilton, Lincoln and Mosgiel, not precise locations.
sampling_loc <-  data.frame(
  long = c(175.32939, 172.48654, 170.34754),
  lat = c(-37.78259,-43.63980, -45.87758),
  names = c('Ruakura', 'Lincoln', 'Invermay'),
  stringsAsFactors = FALSE)  

# Plot the map
# axis lines on top and right currently not working
ggplot() + 
  geom_polygon(data = nz_map, aes(x=long, y = lat, group = group), fill = NA, colour = 'black' ) + 
  coord_fixed(1.45) + 
  geom_point(data = sampling_loc, aes(x = long, y = lat), colour = 'red', size = 3.5) +
  theme_classic() + 
  theme(axis.line.x.top = element_line(), axis.line.y.right = element_line()) +
  labs(x = 'Longitude', y = 'Latitude')


