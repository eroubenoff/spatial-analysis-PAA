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
T.known <- CA_shp[!is.na(CA_shp$e0),] %>% rownames() %>% as.numeric() %>% length() # Indeces of known tracts
T.unk <- CA_shp[is.na(CA_shp$e0),] %>% rownames() %>% as.numeric() %>% length() # Indeces of unknown tracts

m1_data <- list(
  e0.known = e0[!is.na(e0), "e0"]
)

m1_constants <- list(
  S = 1,
  T.known = T.known,
  T.unk = T.unk
)

m1_inits <- list(
  s.mu = array(0, c(S)),
  s.sigma = array(1, c(S)),
  e0.unk = array(0, c(T.unk))
)

m1 <- nimbleCode({
  for (s in 1:S) {
    for (t in 1:T.known) {
      e0.known[t] ~ dnorm(s.mu[s], s.sigma[s])
    }
    
      # Group means 
      s.mu[s] ~ dunif(0, 100)
      s.sigma[s] ~ dunif(0.001, 1)
      
    for (t in 1:T.unk){
      # Prediction at unknowns
      # In this super simple model each unknown one is just the state mean
      e0.unk[t] <- s.mu[s]
    }
  }
})



m1_mcmc <- nimbleMCMC(code = m1,
                      inits = m1_inits,
                      constants = m1_constants,
                      data = m1_data,
                      nchains = 2, 
                      niter = 1000,
                      summary = T,
                      WAIC = T,
                      monitors = c('s.mu', 's.sigma'))
m1_mcmc


#### Slightly more complex model: Tracts in Counties in States (no covariates or adjacency)
#### And let's make things a little more interesting and do 2 states: CA and DE (for contrast)
rm(list = ls())
CA_A <- read_csv("./data/CA_A.CSV")
CA_A <- CA_A %>% rename(e0 = "e(0)",
                        se.e0 = "se(e(0))")
CA_shp <- tracts("CA", cb = F) %>% st_as_sf() %>% st_drop_geometry() # Remove this last bit to get the shapefile
CA_shp <- left_join(CA_shp, CA_A, by = c("GEOID" = "Tract ID"))

DE_A <- read_csv("./data/DE_A.CSV")
DE_A <- DE_A %>% rename(e0 = "e(0)",
                        se.e0 = "se(e(0))")
DE_shp <- tracts("DE", cb = F) %>% st_as_sf() %>% st_drop_geometry()
DE_shp <- DE_shp %>% mutate(GEOID = as.numeric(GEOID)) 
DE_shp <- left_join(DE_shp, DE_A, by = c("GEOID" = "Tract ID"))

S <- 2
# Extract the e0 values for each. *_e0 is the data and *_T.known/unk
CA_e0 <- CA_shp[, c("e0", "GEOID", "COUNTYFP")]
CA_e0 <- CA_e0 %>% 
  mutate(GEOID = as.numeric(GEOID), COUNTYFP = as.numeric(COUNTYFP))
CA_counties <- CA_e0 %>% 
  pull(COUNTYFP) %>% 
  unique() %>% 
  sort()
CA_T.known <- CA_e0[!is.na(CA_e0$e0),] 
CA_T.unk <- CA_e0[is.na(CA_e0$e0),] 
DE_e0 <- DE_shp[, c("e0", "GEOID", "COUNTYFP")]
DE_e0 <- DE_e0 %>% 
  mutate(GEOID = as.numeric(GEOID), COUNTYFP = as.numeric(COUNTYFP))
DE_counties <- DE_e0 %>% 
  pull(COUNTYFP) %>% 
  unique() %>% 
  sort()
DE_T.known <- DE_e0[!is.na(DE_e0$e0),] 
DE_T.unk <- DE_e0[is.na(DE_e0$e0),] 

# Inspect the objects
head(CA_e0)
head(DE_e0)
(CA_counties)
(DE_counties)
head(CA_T.known)
head(CA_T.unk)
dim(CA_T.known)
dim(CA_T.unk)
head(DE_T.known)
head(DE_T.unk)
dim(DE_T.known)
dim(DE_T.unk)
# We're going to use the known counties to estimate the unknown counties

# Need to get in the right format
# indices: [state, county, tract]
S <- 2
C <- c(nrow(CA_counties), nrow(DE_counties))
T <- c(nrow(CA_T.known), nrow(CA_T.unk))
CA_T.known$STATE <- 1
CA_T.unk$STATE <- 1
DE_T.known$STATE <- 2
DE_T.unk$STATE <- 1
# reorder columns
s <- c("STATE", "COUNTYFP", "GEOID", "e0")
CA_T.known <- CA_T.known[,s]
CA_T.unk <- CA_T.unk[,s]
DE_T.known <- DE_T.known[,s]
DE_T.unk <- DE_T.unk[,s]

# Make this into a list of dfs, where 
CA_T.known %>% select(-STATE) %>% split(CA_T.known$COUNTYFP)


for (s in 1:S){
  for (c in e0.known[e0.known$STATE == s, "COUNTYFP"]){
    e0.known.array[s, c, ] <- e0.known %>% filter(STATE == s, COUNTYFP == c) %>% select(GEOID, e0) %>% as.matrix()
  }
}

# I'm getting stuck with the two states thing.  Let's just stick to one state for now

m2_data <- list(
  e0.known = CA_T.known %>% select(-STATE)
)

m2_constants <- list( 
  S = 1,
  n.T = CA_T.known$GEOID %>% unique(),
  n.C = CA_T.known$COUNTYFP %>% unique()
)

m2_inits <- list(
  s.mu = array(0, c(S)),
  s.sigma = array(1, c(S)),
  e0.unk = array(0, c(T.unk))
)

m2 <- nimbleCode({
    for (c in n.C){
      for (t in 1:T.known) {
      e0.known[t] ~ dnorm(s.mu[s], s.sigma[s])
      }
    
    # Group means 
    s.mu[s] ~ dunif(0, 100)
    s.sigma[s] ~ dunif(0.001, 1)
    
    for (t in 1:T.unk){
      # Prediction at unknowns
      # In this super simple model each unknown one is just the state mean
      e0.unk[t] <- s.mu[s]
    }
  }
})



