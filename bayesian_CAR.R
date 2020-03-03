###############################################################################
# Working on the full bayesian model (This time with autoregressive priors!)
# One state at a time
# tracts within states, no time smoothing

# Note about CAR modelling
# We can put CAR priors at the county level, the tract level, or both.
# I think it makes most sense to put it at the county level first for testing.

library(nimble)
library(tidyverse)
library(reshape2)
library(tictoc)
library(tigris)
library(sf)


setwd("~/kriging_PAA")

rm(list = ls())
gc()

message("~~~~~~~~~~~~~~~~~~~~~Starting Run~~~~~~~~~~~~~~~~~~~~~~~~~~")
message(Sys.time())
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

#### Create PCs ####

pcs <- as.matrix(read.csv("pcs.csv"))[,2:4]


##### Format Data ####
age.v <- c(0, 1, 5, 15, 25, 35, 45, 55, 65, 75, 85)  # Ages we need
n.v   <- c(1, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9)
CA_1 <- suppressMessages({read_csv("./data/CA_B_1.CSV")})
CA_2 <- suppressMessages({read_csv("./data/CA_B_2.CSV") %>% mutate(TRACT2KX = as.numeric(TRACT2KX))})
CA <- bind_rows(CA_1, CA_2)
rm(CA_1, CA_2)
CA <- CA %>% 
  mutate(mx = `nd(x)` / `nL(x)`)

CA <- CA %>%
  rename(tractl = `Tract ID`,
         state = STATE2KX,
         county = CNTY2KX,
         tract = TRACT2KX,
         age = `Age Group`,
         nqx = `nq(x)`,
         lx = `l(x)`,
         ndx = `nd(x)`,
         nLx = `nL(x)`,
         Tx  = `T(x)`,
         ex = `e(x)`,
         se_nqx = `se(nq(x))`,
         se_ex = `se(e(x))`,
         mx = mx)


CA <- CA %>% mutate(age = rep(age.v, nrow(CA)/11),
                    n = rep(n.v, nrow(CA)/11),
                    # tract = as.numeric(tract),
                    county = as.numeric(county),
                    ndx = as.numeric(ndx))

# For testing: limit to the bay area
# 1: Alameda, 13: Contra Costa, 55: Napa, 97: Sonoma, 41: Marin, 75: San Francisco, 95: Solano, 81: San Mateo, 85: Santa Clara
bayarea_counties <- c("Alameda" = 1, "Contra Costa" = 13, "Napa" = 55,
                      "Sonoma" = 97, "Marin" = 41, "San Francisco" = 75,
                      "Solano" = 95, "San Mateo" = 81, "Santa Clara" = 85) %>% sort

CA <- CA %>% filter(county %in% bayarea_counties)

# Go from tract numbers to tract index within county
CA <- CA %>% mutate(county_r = county) %>%
  group_by(county_r) %>%
  group_map(~mutate(.,tract_index=group_indices(.,tract))) %>% 
  bind_rows()


#### Creating the Spatial Adjacencies ####
# Process: download the shapefile and subset the counties
# Determine adjacency matrix

shp <- counties("CA", cb = T) %>% 
  st_as_sf() %>% 
  mutate(COUNTYFP = as.numeric(COUNTYFP))  %>% 
  filter(COUNTYFP %in% bayarea_counties) %>%
  arrange(COUNTYFP)        # Make sure they're in the same order (necessary for adj matrix)

adj <- st_touches(shp, sparse =  T) 
num <- sapply(adj, length)
adj <- adj %>% unlist
weights <- adj/adj


# Tract map (for reconstructing after simulation)
tract_map <- CA %>%
  select(state, county, tract, tractl, tract_index) %>%
  unique()

ndx <- CA %>% select(county, tract_index, age, ndx) %>% 
  acast(age~county~tract_index, fill = 0, value.var = "ndx")
nLx <- CA %>% select(county, tract_index, age, nLx) %>% 
  acast(age~county~tract_index, fill = 0, value.var = "nLx")

# Need to get number of tracts in each county
n.T <- CA %>% 
  select(tract, county) %>% 
  unique %>% 
  group_by(county) %>% 
  summarize(n = n()) %>%
  pull(n)

X <- length(age.v)
C <- length(n.T)
n.Tmax <- max(n.T)

data <- list(
  ndx = ndx,
  Yx = pcs
)
constants <- list(
  nLx = nLx,
  X = X,
  C = C,
  n.T = n.T,
  n.Tmax = n.Tmax,
  num = as.vector(num),
  adj = as.vector(adj),
  weights = as.vector(weights)
)




#### Model Definition (nonspatial) ####
# Nesting tracts in counties:
code <- nimbleCode({
  for (x in 1:X) {              # Age
    for (c in 1:C) {            # County
      for (t in 1:n.T[c]){        # Tract
        ndx[x, c, t] ~ dpois(mu[x, c, t])
        mu[x, c, t] <- mx[x, c, t] * nLx[x, c, t]
        mx[x, c, t] <- exp(logmx[x, c, t])
        logmx[x, c, t] <- beta[c, t, 1] * Yx[x, 1] + 
          beta[c, t, 2] * Yx[x, 2] +
          beta[c, t, 3] * Yx[x, 3] +
          u[x, c, t] + 
          s[c]               # Spatial term
      }
      for (t in n.T[c]+1:n.Tmax){ # Because counties have different numbers of tracts
        mx[x, c, t] <- 0
      }
    }
  }
  
  # Spatial term
  for(c in 1:C) {mu[c] <- 0}
  s[1:C] ~ dcar_normal(adj=adj, weights = weights, num = num, tau = tau, zero_mean = 1)
  tau ~ dgamma(0.001, 0.001)
  
  # Priors on beta
  for (c in 1:C){          # County
    for (t in 1:n.T[c]){   # Tract
      for(i in 1:3){       # Principal Component i
        beta[c, t, i] ~ dnorm(mu.beta[c, i], tau.beta[c, i])
      }
    }
    for (t in n.T[c]+1:n.Tmax){ # Because different numbers of tracts
      for(i in 1:3){ 
        beta[c, t, i] <- 0
      }
    }
    # Prior on precision
    for(i in 1:3){
      tau.beta[c,i] <- pow(sigma.beta[c,i], -2)
      sigma.beta[c,i] ~ dunif(0,40)
    }
  }
  
  # Priors on u
  for (c in 1:C){
    for (t in 1:n.T[c]){
      for (x in 1:X){
        u[x, c, t] ~ dnorm(0, tau.u[x, c])
      }
    }
  }
  
  # Priors on mu.beta
  for(c in 1:C){     # Counties
    for (i in 1:3){  # PCs
      mu.beta[c, i] ~ dnorm(0, tau.mu[c, i])
      tau.mu[c, i] <- pow(sigma.mu[c, i], -2)
      sigma.mu[c, i] ~ dunif(0, 40)
    }
  }
  
  # Priors on tau.u 
  for(x in 1:X){
    for (c in 1:C){
      tau.u[x, c] <- pow(sigma.u[x, c], -2)
      sigma.u[x, c] ~ dunif(0, 0.1)
    }
  }
})

inits <- list(
  sigma.beta = array(1, c(C, 3)),
  tau.beta = array(1, c(C, 3)),
  sigma.mu = array(1, c(C, 3)),
  tau.mu = array(1, c(C, 3)),
  sigma.u = array(.05, c(X, C)),
  tau.u = array(1, c(X, C)),
  mu.beta = array(0, c(C, 3)),
  u = array(0, c(X, C, n.Tmax)), 
  beta = array(0, c(C, n.Tmax, 3)), 
  logmx = array(0, c(X, C, n.Tmax)), 
  mx = array(0, c(X, C, n.Tmax)),
  mu = array(0, c(X, C, n.Tmax)),
  s = array(0, c(C))
)

message("~~~~~~~~~~~~~~~~~~~~~Initializing~~~~~~~~~~~~~~~~~~~~~~~~~~")
message(Sys.time())
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

mod <- nimbleModel(code = code, name = 'mod', constants = constants,
                   data = data, inits = inits, calculate = F, check = F)

message("~~~~~~~~~~~~~~~~~~~~~Running Chains~~~~~~~~~~~~~~~~~~~~~~~~")
message(Sys.time())
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")



monitors <- c("mx", "u")


# In one step:
tic()
mcmc_out <- nimbleMCMC(code = mod,
                       data = data,
                       constants = constants,
                       inits = inits,
                       monitors = monitors,
                       nchains = 2,
                       # nburnin = 2000,
                       niter = 1000,
                       check = F,
                       summary = T)
toc()

save(mcmc_out, file = paste0("./mcmc_runs/MCMC_OUT_spatial_", Sys.time(), ".RData"))





