#---- Tract Bayesian Analysis -------------------------------------------------
#---- Whole dataset
#---- Ethan Roubenoff, 23 Jan 2020
#------------------------------------------------------------------------------
#---- Model 1 is just tracts grouped in states
#---- Will run a test on just CA before expanding to the whole US
# --- https://link.springer.com/article/10.1023/A:1009614003692
install.packages("nimble")
library(nimble)
library(tidyverse)
library(tigris)
library(sf)
rm(list = ls())

CA_A <- read_csv("./data/CA_A.CSV")
CA_A <- CA_A %>% rename(e0 = "e(0)",
                   se.e0 = "se(e(0))")
CA_shp <- tracts("CA", cb = F) %>% st_as_sf() %>% st_drop_geometry() # Remove this last bit to get the shapefile
CA_shp <- left_join(CA_shp, CA_A, by = c("GEOID" = "Tract ID"))
summary(CA_shp$e0)

# Our goal is to estimate the missing 541s

#### Simplest model: Tracts in States ####
# e0[t, s] ~ N(s.mu[s], s.sigma[s])
# s.mu[s] ~ U(0, 100)
# s.sigma[s] ~ U(0, 1)

S <- 1
e0 <- CA_shp[, c("e0", "GEOID")]
T.known <- CA_shp[!is.na(CA_shp$e0),] %>% rownames() %>% as.numeric() # Indeces of known tracts
T.unk <- CA_shp[is.na(CA_shp$e0),] %>% rownames() %>% as.numeric()  # Indeces of unknown tracts

m1_data <- list(
  S = 1,
  e0 = e0,
  T.known = as.array(T.known),
  T.unk = as.array(T.unk)
)

m1 <- nimbleCode({
  for (s in 1:S) {
    for (t in T.known) {
      e0[t] ~ dnorm(s.mu[s], s.sigma[s])
      # Priors 
      s.mu[s] ~ dunif(0, 100)
      s.sigma[s] ~ dunif(0, 1)
    }
    
    for (t in T.unk){
      # Prediction at unknowns
      e0[t] ~ dnorm(s.mu[s], s.sigma[s])
    }
  }
})

inits <- list(
  s.mu = array(0, c(S)),
  s.sigma = array(0, c(S))
)

m1_mcmc <- nimbleMCMC(code = m1,
                      inits = inits,
                      data = m1_data,
                      nchains = 2, 
                      niter = 1000,
                      summary = T,
                      WAIC = T,
                      monitors = c('s.mu', 's.sigma', 'e0'))

