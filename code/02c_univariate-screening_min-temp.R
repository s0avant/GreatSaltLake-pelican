#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                  Min Spring Temp With Time Lag - Nimble Code                 #
#                                 05/19/2023                                   #
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
  
  beta0 ~ dnorm(0, sd = 1.75) # no lag
  beta1 ~ dnorm(0, sd = 1.75) # one year lag
  beta3 ~ dnorm(0, sd = 1.75) # three year lag
  beta4 ~ dnorm(0, sd = 1.75) # four year lag
  
  beta.pow0 ~ dnorm(0, sd = 1.75) # no lag
  beta.pow1 ~ dnorm(0, sd = 1.75) # one year lag
  beta.pow3 ~ dnorm(0, sd = 1.75) # three year lag
  beta.pow4 ~ dnorm(0, sd = 1.75) # four year lag
  
  sigma.proc ~ dunif(0, 10)    # Prior for sd of state process
  sigma.obs ~ dunif(0, 10)    # Prior for sd of observation process
  tau.obs <- pow(sigma.obs, -2)
  
  # Likelihood
  # State process
  for (t in 1:(T-lag)){
    z.proc[t] ~ dnorm(0, sd = 1)
    logN.est[t + 4] <- beta + 
      beta.yr * year[t] +
      beta0 * min_temp[t + 4] +
        beta.pow0 * pow(min_temp[t + 4], 2) + # no lag
      beta1 * min_temp[t + 3] +
        beta.pow1 * pow(min_temp[t + 3], 2) + # one year lag
      beta3 * min_temp[t + 1] +
        beta.pow3 * pow(min_temp[t + 1], 2) + # three year lag
      beta4 * min_temp[t] +
        beta.pow4 * pow(min_temp[t], 2) + # four year lag
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
min_temp <- readRDS("data/min_temp.RDS")

#~ Filter to 1985 breakpoint data ----
awpe_break <- filter(awpe, year %in% 1985:2020)
min_temp_break <- filter(min_temp, year %in% 1985:2020)

counts <- awpe_break[,2] # colony counts
year <- (awpe_break[,1] - min(awpe_break[,1]))/100

lag <- 4 # max number of years to test for time lagged effects

#~ Calculate residuals to remove possible year effect ----
min_temp.r <- resid(lm(mean ~ year, data = min_temp_break))

n_model <- nimbleModel(code = n_code,
                       constants = list(T = length(year), lag = lag),
                       inits = list(beta = runif(1, -0.5, 0.5),
                                    beta.yr = runif(1, -0.25, 0),
                                    beta0 = runif(1, -0.25, 0),
                                    beta1 = runif(1, -0.25, 0),
                                    beta3 = runif(1, -0.25, 0),
                                    beta4 = runif(1, -0.25, 0),
                                    beta.pow0 = runif(1, -0.25, 0),
                                    beta.pow1 = runif(1, -0.25, 0),
                                    beta.pow3 = runif(1, -0.25, 0),
                                    beta.pow4 = runif(1, -0.25, 0),
                                    sigma.proc = runif(1, 0, 1),
                                    sigma.obs = runif(1, 0, 1),
                                    logN.est = c(runif(lag, 8, 9),
                                                 rep(NA, length(year) - lag)),
                                    z.proc = rnorm(length(year) - lag),
                                    y = rep(NA, length(year))),
                       data = list(y = log(counts),
                                   min_temp = min_temp.r,
                                   year = year
                                   ))

n_model$initializeInfo()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Configure model, build MCMC, compile ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
n_config <- configureMCMC(n_model, print = TRUE, useConjugacy = FALSE,)
n_config$addMonitors("N.est")

target <- c("beta",
            "beta.yr",
            "beta0",
            "beta1",
            "beta3",
            "beta4",
            "beta.pow0",
            "beta.pow1",
            "beta.pow3",
            "beta.pow4")

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

ifelse(max(awpe.ssm.nimble$Rhat) >= 1.12,
       "Rhat indicates convergence failure",
       "Rhat indicates successful convergence")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save outputs ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MCMCvis::MCMCdiag(n_samples$samples, round = 2,
                  file_name = "output/model-output_post-break_lag_min-temp.txt")

MCMCvis::MCMCtrace(n_samples$samples,
                   params = target,
                   filename =
                     "output/model-output_post-break_lag_min-temp.pdf")

write_rds(awpe.ssm.nimble, "output/model-output_post-break_lag_min-temp.RDS")

write_rds(n_samples$samples,
        "output/posteriors/posterior-samples_post-break_lag_min-temp.RDS")
