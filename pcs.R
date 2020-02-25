# Create principal components from the USAMDB (adapted from Monica)
library(tidyverse)
rm(list = ls())
setwd("~/kriging_PAA")

age.v <- c(0, 1, 5, 15, 25, 35, 45, 55, 65, 75, 85)  # Ages we need
age.i <- age.v + 1                                   # Indeces for those ages

# Load them all in!
state_names <- list.dirs("./lifetables/States", full.names = F, recursive = F)

for (s in state_names){
  f <- paste0("./lifetables/States/", s, "/", s, "_bltper_1x1.csv")
  
  # Create object if it exists
  if (!exists("lt")){
    lt <- suppressMessages(read_csv(f))
  } else {
    lt <- suppressMessages(bind_rows(lt, read_csv(f)))
  }
}

# Select years 1980-2010
lt <- lt %>% filter(Year > 1980 & Year < 2010) %>% select(PopName, Year, Age, mx)

lt_mat <- lt %>%  
  mutate(logmx = log(mx),
         StateYear = paste0(PopName, Year)) %>%
  select(-mx, -PopName, -Year) %>%
  pivot_wider(names_from = StateYear, values_from = logmx)  %>%
  select(-Age) %>%
  as.matrix() 

# If there are any values of logmx that are infinite, set them to 0 
lt_mat[is.infinite(lt_mat)] <- 0

# Do PCA
pcs <- base::svd(lt_mat)$u[,1:3]

plot(age.v, pcs[age.i,1], type = "l", ylim = c(-0.5, .25))
lines(age.v, pcs[age.i,2], type = "l")
lines(age.v, pcs[age.i,3], type = "l")

write.csv(pcs[age.i, 1:3], "pcs.csv")
