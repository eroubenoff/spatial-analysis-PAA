#---- Tract Kriging -----------------------------------------------------------
#---- Whole dataset
#---- Ethan Roubenoff, 22 Jan 2020
#------------------------------------------------------------------------------

library(tidyverse)
library(tidycensus)
library(tigris)
library(sf)
library(tmap)
library(gstat)
library(raster)
library(fasterize)
library(magrittr)

setwd("~/kriging_PAA/")

# California example
CA_A <- read_csv("./data/CA_A.CSV")
CA_shp <- tracts("CA") %>% st_as_sf()
CA_shp <- left_join(CA_shp, CA_A, by = c("GEOID" = "Tract ID"))
CA_shp %<>% rename(e0 = "e(0)",
                   se.e0 = "se(e(0))")

# There are 7516 tracts in CA_A and 8057 tracts in CA_shp.  
# There are 541 "missing" tracts in the joined CA_shp, plotted below:
tmap_mode("view")
CA_shp %>% filter(is.na(`e(0)`)) %>% qtm()
# These are ostensibly tracts with populations too low to be estimated.  
# I won't sweat these for now.

#---- PROCEDURE ---------------------------------------------------------------
#---- Drop a random subset of tracts for prediction
#---- Create raster of the known pts
#---- Krige at the center of each unknown polygon


drop_tracts <- function(CA_shp, n){
  #' Takes the shapefile of CA tracts and randomly drops n of them
  #' Returns two objects:
  #'      - shapefile of known locations
  #'      - shapefile of 'unknown' locations
  #'      
  
  n_tracts <- nrow(CA_shp)
  dropped <- sample(1:n_tracts, n)
  
  CA_shp %<>% filter(!is.na(e0))
  
  CA_kno <- CA_shp %>% slice(-dropped)
  CA_unk <- CA_shp %>% slice(dropped)
  
  return(
    list("CA_kno" = CA_kno,
         "CA_unk" = CA_unk )
  )
}

obj <- drop_tracts(CA_shp, 100)
CA_kno <- obj$CA_kno
CA_unk <- obj$CA_unk
rm(obj)

tmap_mode("view")
n.row <- 100
n.col <- 100

#---- Simple kriging
CA_kno_rast <- raster(ncol = n.col, nrow = n.row, ext = extent(CA_kno))
CA_kno_rast <- fasterize(CA_kno, CA_kno_rast, field = "e0")
qtm(CA_kno_rast)
CA_unk_pts <- CA_unk %>% st_centroid()
ex.vgm <- variogram(layer~1, rasterToPoints(CA_kno_rast, spatial = T), width = 10)
ex.fit <- fit.variogram(ex.vgm, model=vgm("Sph"))
plot(ex.vgm, ex.fit)
ex.kriged <- krige(layer ~ 1, rasterToPoints(CA_kno_rast, spatial = T), CA_unk_pts, model=ex.fit)
ex.kriged$var1.sd <- sqrt(ex.kriged$var1.var)
ex.kriged

#---- Compare with the given numbers
ex.kriged %<>% st_join(CA_shp)
plot(ex.kriged$var1.pred, ex.kriged$e0)
mse <- sum((ex.kriged$var1.pred - ex.kriged$e0)^2) / nrow(ex.kriged)
mse
rmsd <- sqrt(mse)
rmsd
msevar <- sum((ex.kriged$var1.sd - ex.kriged$se.e0)^2)/nrow(ex.kriged )
msevar

#---- Cokriging EX with median hh income and pct white
# t <- load_variables(2014, "acs5", cache = T)
list()
CA_shp <- tidycensus::get_acs(geography = "tract", state = "CA", 
                              variables = c("med_inc" = "B21004_001","tot_pop" =  "B01003_001", "tot_white"= "B02001_002"), 
                              geometry = T, output = "wide")
CA_shp <- left_join(CA_shp, CA_A, by = c("GEOID" = "Tract ID"))
CA_shp %<>% rename(e0 = "e(0)",
                   se.e0 = "se(e(0))")

obj <- drop_tracts(CA_shp, 100)
CA_kno <- obj$CA_kno
CA_unk <- obj$CA_unk
rm(obj)

tmap_mode("view")
n.row <- 100
n.col <- 100

CA_kno_rast <- raster(ncol = n.col, nrow = n.row, ext = extent(CA_kno))
CA_kno %>% pivot_longer(cols = c(med_incE, med_incM, tot_popE, tot_popM, tot_whiteE, tot_whiteM, e0, se.e0))
CA_kno_rast <- fasterize(CA_kno, CA_kno_rast, field = "e0")
qtm(CA_kno_rast)
CA_unk_pts <- CA_unk %>% st_centroid()
ex.vgm <- variogram(layer~1 , rasterToPoints(CA_kno_rast, spatial = T), width = 10)
ex.fit <- fit.variogram(ex.vgm, model=vgm("Sph"))
plot(ex.vgm, ex.fit)
ex.kriged <- krige(layer ~ med_incE + tot_whiteE/tot_popE, rasterToPoints(CA_kno_rast, spatial = T), CA_unk_pts, model=ex.fit)
ex.kriged$var1.sd <- sqrt(ex.kriged$var1.var)
ex.kriged
