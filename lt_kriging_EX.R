# https://rpubs.com/nabilabd/118172 Inspiration
# Turn tract data into a lattice, then krige at the rest of the lattice
# Ex with MA:
# download.file("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NVSS/USALEEP/CSV/MA_B.CSV", "MA_lt.csv")

library(tidyverse)
library(tidycensus)
library(tigris)
library(sf)
library(tmap)
library(gstat)

MA_lt <- read_csv("MA_lt.csv")
MA_shp <- tracts("MA") %>% st_as_sf() # Change this to 2010 tract boundaries

MA_ex <- MA_lt %>% dplyr::select("Tract ID", "Age Group", "e(x)", "se(e(x))") %>% 
  rename(age = "Age Group", GEOID = "Tract ID", ex = "e(x)", ex_se = "se(e(x))") %>%
  filter(age == "Under 1") %>%
  mutate(GEOID = as.character(GEOID)) %>%
  dplyr::select(-age)
MA_shp <- MA_shp %>% left_join(MA_ex)

tm_shape(MA_shp) + tm_fill(col = "ex") 

# Generate centroids
MA_pt <- MA_shp %>% st_centroid()

MA_known <- MA_pt %>% filter(!is.na(ex))
MA_unk <- MA_pt %>% filter(is.na(ex))

tm_shape(MA_shp) + 
  tm_borders(col = "grey") +
  tm_shape(MA_known) +
  tm_dots(col = "ex") + 
  tm_shape(MA_unk) + 
  tm_dots(col = "black")
# Create variogram
ex.vgm <- variogram(ex~1, MA_known) # Add standard errors?
# ex.vgm <- variogram(ex~ex_se, MA_known) # Add standard errors?
ex.fit <- fit.variogram(ex.vgm, model=vgm("Sph"))
plot(ex.vgm, ex.fit)

ex.kriged <- krige(ex ~ 1, MA_known, MA_unk, model=ex.fit)
ex.kriged$var1.sd <- sqrt(ex.kriged$var1.var)
ggplot(data = ex.kriged) + # Check for correlation in errors
  geom_point(aes(var1.pred, var1.sd))

tmap_mode("view")
tm_shape(MA_shp) +
  tm_borders(col = "grey") +
  tm_shape(MA_known) +
  tm_dots(col = "ex") +
  tm_shape(ex.kriged) +
  tm_dots(col = "var1.pred", size = "var1.sd")



# Trying as a grid
# Create a lattice of points and spatial join with polygons
MA_shp <- MA_shp %>% 
  st_transform(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
MA_grid <-  MA_shp %>%
  st_make_grid(cellsize = 0.1, what = "centers") %>% # grid of points
  st_intersection(MA_shp)   %>%                        # only within the polygon
  st_sf() 
MA_grid <- st_join(MA_grid, MA_shp, join = st_within)
tm_shape(MA_grid) + tm_dots(col = "ex")

# Subset known and unknown
MA_known <- MA_grid %>% filter(!is.na(ex))
MA_unk <- MA_grid %>% filter(is.na(ex))
# View the distribution
tm_shape(MA_shp) + 
  tm_borders(col = "grey") +
  tm_shape(MA_known) +
  tm_dots(col = "ex") + 
  tm_shape(MA_unk) + 
  tm_dots(col = "black")

# create variogram
# Ordinary model: ex~1
# ex.vgm <- variogram(ex~1, MA_known) 
# Universal model: ex~sd
ex.vgm <- variogram(ex~ex_se, MA_known) 
ex.fit <- fit.variogram(ex.vgm, model=vgm("Sph"))
plot(ex.vgm, ex.fit)

# Krige
# Ordinary:
# ex.kriged <- krige(ex ~ 1, MA_known, MA_unk, model=ex.fit)
# Universal:
ex.kriged <- krige(ex ~ ex_se, MA_known, MA_unk, model=ex.fit)
ex.kriged$var1.sd <- sqrt(ex.kriged$var1.var)
ggplot(data = ex.kriged) + # Check for correlation in errors
  geom_point(aes(var1.pred, var1.sd))
tm_shape(MA_shp) +
  tm_borders(col = "grey") +
  tm_shape(MA_known) +
  tm_dots(col = "ex") +
  tm_shape(ex.kriged) +
  tm_dots(col = "var1.pred", size = "var1.sd")






