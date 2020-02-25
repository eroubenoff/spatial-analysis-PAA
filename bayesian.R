###############################################################################
# Working on the full bayesian model
# One state at a time
# tracts within states, no time smoothing

library(nimble)
library(tidyverse)
library(reshape2)
library(tictoc)

rm(list = ls())
gc()


#### Create PCs ####

CA_1 <- read_csv("./data/CA_B_1.CSV")
CA_2 <- read_csv("./data/CA_B_2.CSV") %>% mutate(TRACT2KX = as.numeric(TRACT2KX))
CA <- bind_rows(CA_1, CA_2)
rm(CA_1, CA_2)
CA <- CA %>% 
  mutate(mx = `nd(x)` / `nL(x)`)

age.v <- c(0, 1, 5, 15, 25, 35, 45, 55, 65, 75, 85) 
n.v   <- c(1, 4, 9, 9, 9, 9, 9, 9, 9, 9, NA)

CA_mat <- CA %>%  
  select(`Age Group`, `Tract ID`, mx) %>%
  mutate(logmx = log(mx)) %>%
  select( -mx) %>%
  pivot_wider(names_from = `Tract ID`, values_from = logmx)  %>%
  select(-`Age Group`) %>%
  as.matrix() 

# Drop infinite values 
sum(sapply(CA_mat, is.infinite))

# Pull out the three largest PCs
pcs <- base::svd(CA_mat)$u[,1:3]

ggplot(as_data_frame(pcs)) + 
  geom_line(aes(x = age.v, y = V1), linetype = "solid") +
  geom_line(aes(x = age.v, y = V2), linetype = "dashed") +
  geom_line(aes(x = age.v, y = V3), linetype = "dotted") 

rm(CA_mat)



##### Format Data ####
CA_1 <- read_csv("./data/CA_B_1.CSV")
CA_2 <- read_csv("./data/CA_B_2.CSV") %>% mutate(TRACT2KX = as.numeric(TRACT2KX))
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

# For testing: limit to the first 4 counties
CA <- CA %>% filter(county %in% c(1, 3, 5, 7))

# Go from tract numbers to tract index within county
CA <- CA %>% mutate(county_r = county) %>%
  group_by(county_r) %>%
  group_map(~mutate(.,tract_index=group_indices(.,tract))) %>% 
  bind_rows()

# Tract map (for reconstructing after simulation)
tract_map <- CA %>%
  select(state, county, tract, tractl, tract_index) %>%
  unique()

ndx <- CA %>% select(county, tract_index, age, ndx) %>% acast(age~county~tract_index, fill = 0)
nLx <- CA %>% select(county, tract_index, age, nLx) %>% acast(age~county~tract_index, fill = 0)

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
  n.Tmax = n.Tmax
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
                          u[x, c, t]
      }
      for (t in n.T[c]+1:n.Tmax){ # Because counties have different numbers of tracts
        mx[x, c, t] <- 0
      }
    }
  }
  
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
  mu = array(0, c(X, C, n.Tmax))
)

mod <- nimbleModel(code = code, name = 'mod', constants = constants,
                         data = data, inits = inits)
# mod$initializeInfo()
# modMCMC <- buildMCMC(mod)
# # For testing only
# # runMCMC_samples <- runMCMC(modMCMC, niter = 5)
# Cmod <- compileNimble(mod)
# CmodMCMC <- compileNimble(modMCMC, project = mod)
# runMCMC_samples <- runMCMC(CmodMCMC, nburnin = 1000, niter = 10000)
# monitors <- c("y.xta", "mu.xta", "mx.xta", "logmx.xta", "beta.ta", "u.xta")


# In one step:
mcmc_out <- nimbleMCMC(code = mod,
           data = data,
           constants = constants,
           inits = inits,
           summary = T)






