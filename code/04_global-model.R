#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#             Gunnison Colony Count SSM Global Model - Nimble Code             #
#                                 02/02/2024                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(nimble)
library(dplyr)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Nimble model code ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
n_code <- nimbleCode({
  # Priors
  logN.est[1] ~ dnorm(9.7, sd = 1.75)    # Prior for initial population size
  logN.est[2] ~ dnorm(9.7, sd = 1.75)    # Prior for initial population size
  logN.est[3] ~ dnorm(9.7, sd = 1.75)    # Prior for initial population size
  logN.est[4] ~ dnorm(9.7, sd = 1.75)    # Prior for initial population size
  
  beta ~ dnorm(0, sd = 1.75)
  
  beta.yr ~ dnorm(0, sd = 1.75)
  
  beta.DD1 ~ dnorm(0, sd = 1.75)
  beta.DD4 ~ dnorm(0, sd = 1.75)
  beta.water1 ~ dnorm(0, sd = 1.75)
  beta.temp0 ~ dnorm(0, sd = 1.75)
  beta.temp.pow0 ~ dnorm(0, sd = 1.75)
  beta.bridge1 ~ dnorm(0, sd = 1.75)
  
  sigma.proc ~ dgamma(0.1, 0.1)    # Prior for sd of state process
  sigma.obs ~ dgamma(0.1, 0.1)    # Prior for sd of observation process
  tau.obs <- pow(sigma.obs, -2)
  
  # Likelihood
  # State process
  for (t in 1:(T-lag)){
    z.proc[t] ~ dnorm(0, sd = 1)
    logN.est[t + 4] <- beta +
      
      beta.yr * year[t + 4] +
      
      (1 + beta.DD1) * logN.est[t + 3] +
      beta.DD4 * logN.est[t] +
      
      beta.water1 * waterlvl[t + 3] +
      
      beta.temp0 * min_temp[t + 4] +
      beta.temp.pow0 * pow(min_temp[t + 4], 2) +
      
      beta.bridge1 * landbridge[t + 3] +
      
      z.proc[t] * sigma.proc
  }
  
  # Observation process
  for (t in 1:T){
    y[t] ~ dnorm(logN.est[t], tau.obs)
  }
  
  # Population sizes on real scale
  for (t in 1:T){
    N.est[t] <- exp(logN.est[t])
  }
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data cleaning & initialize model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~ Read in data ----
awpe <- read.csv("data/colony_count_simple.csv", na.strings = "", header = TRUE)
waterlvl <- readRDS("data/waterlevels.RDS")
min_temp <- readRDS("data/min_temp.RDS")

#~ Filter to 1985 breakpoint data ----
awpe_break <- filter(awpe, year %in% 1985:2020)
waterlvl_break <- filter(waterlvl, year %in% 1985:2020)
min_temp_break <- filter(min_temp, year %in% 1985:2020)

counts <- awpe_break[,2] # colony counts
year <- (awpe_break[,1] - min(awpe_break[,1]))/100

lag <- 4 # max number of years to test for time lagged effects

#~ Calculate residuals
waterlvl.r <- resid(lm(mean ~ year, data = waterlvl_break))
min_temp.r <- resid(lm(mean ~ year, data = min_temp_break))

n_model <- nimbleModel(code = n_code,
                       constants = list(T = length(year), lag = lag),
                       inits = list(beta = runif(1, -0.5, 0.5),
                                    beta.yr = runif(1, -0.25, 0),
                                    beta.DD1 = runif(1, -0.25, 0),
                                    beta.DD4 = runif(1, -0.25, 0),
                                    beta.water1 = runif(1, -0.25, 0),
                                    beta.temp0 = runif(1, -0.25, 0),
                                    beta.temp.pow0 = runif(1, -0.25, 0),
                                    beta.bridge1 = runif(1, -0.25, 0),
                                    sigma.proc = runif(1, 0, 1),
                                    sigma.obs = runif(1, 0, 1),
                                    logN.est = c(runif(lag, 8, 9),
                                                 rep(NA, length(year) - lag)),
                                    z.proc = rnorm(length(year) - lag),
                                    y = rep(NA, length(year))),
                       data = list(y = log(counts),
                                   year = year,
                                   waterlvl = waterlvl.r,
                                   landbridge = waterlvl_break$landbridge,
                                   min_temp = min_temp.r))

n_model$initializeInfo()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Configure model, build MCMC, compile ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
n_config <- configureMCMC(n_model, print = TRUE, useConjugacy = FALSE,)
n_config$addMonitors("N.est")

target <- c("beta",
            "beta.yr",
            "beta.DD1",
            "beta.DD4",
            "beta.water1",
            "beta.temp0",
            "beta.temp.pow0",
            "beta.bridge1")

n_config$removeSamplers(target)
n_config$addSampler(target, type = 'AF_slice') 

#~ Check samplers
n_config$printSamplers()
n_mcmc <- buildMCMC(n_config)

comp_mod <- compileNimble(n_model)
comp_mcmc <- compileNimble(n_mcmc)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fit model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
niter <- 100000
nburnin <- 10000
thin <- 1

system.time(n_samples <- runMCMC(comp_mcmc,
                                 niter = niter,
                                 nburnin = nburnin,
                                 thin = thin, 
                                 nchains = 3,
                                 summary = TRUE,
                                 samplesAsCodaMCMC = FALSE))

awpe.ssm.nimble <- MCMCvis::MCMCsummary(n_samples$samples)

ifelse(max(awpe.ssm.nimble$Rhat) >= 1.1,
       "Rhat indicates convergence failure",
       "Rhat indicates successful convergence")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save outputs ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MCMCvis::MCMCdiag(n_samples$samples,
                  file_name =
                    "output/model-output_post-break_full-with-year.txt")

MCMCvis::MCMCtrace(n_samples$samples,
                   params = target,
                   iter = niter,
                   filename =
                     "output/model-output_post-break_full-with-year.pdf")

write_rds(awpe.ssm.nimble,
          "output/model-output_post-break_full-with-year.RDS")

write_rds(n_samples$samples,
          "output/posteriors/posterior-samples_post-break_full-with-year.RDS")
