#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Gunnison Colony Count SSM Global Model Prediction - Nimble Code       #
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
  
  beta.DD1 ~ dnorm(0, sd = 1.75)
  beta.DD4 ~ dnorm(0, sd = 1.75)
  beta.water1 ~ dnorm(0, sd = 1.75)
  beta.temp0 ~ dnorm(0, sd = 1.75)
  beta.temp.pow0 ~ dnorm(0, sd = 1.75)
  beta.bridge1 ~ dnorm(0, sd = 1.75)
  
  sigma.proc ~ dgamma(0.1, 0.1)    # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  
  sigma.obs ~ dgamma(0.1, 0.1)    # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)
  
  # Likelihood
  # State process
  for (t in 1:(T-lag)){
    z.proc[t] ~ dnorm(0, sd = 1)
    logN.est[t + 4] <- beta +
      
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
    y[t] ~ dnorm(logN.est[t], sd = sigma.obs)
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

awpe_NA <- awpe
awpe_NA$count[awpe_NA$year >= 2013] <- NA # Set to NA to test model prediction capability

#~ Filter to 1985 breakpoint data ----
awpe_break <- filter(awpe_NA, year %in% 1985:2020)
waterlvl_break <- filter(waterlvl, year  %in% 1985:2020)
min_temp_break <- filter(min_temp, year  %in% 1985:2020)

counts <- awpe_break[,2] # colony counts
year <- (awpe_break[,1] - min(awpe_break[,1]))/100

lag <- 4 # max number of years to test for time lagged effects

#~ Normalize variables
waterlvl_center <- (waterlvl_break$mean - mean(waterlvl_break$mean))/
  sd(waterlvl_break$mean)
min_temp_center <- (min_temp_break$mean - mean(min_temp_break$mean))/
  sd(min_temp_break$mean)

n_model <- nimbleModel(code = n_code,
                       constants = list(T = length(year), lag = lag),
                       inits = list(beta = runif(1, -0.5, 0.5),
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
                                    y = c(rep(NA, length(counts[is.na(counts) == FALSE])),
                                          runif(length(year) -
                                                  length(counts[is.na(counts) == FALSE]),
                                                8, 9))), 
                       data = list(y = log(counts),
                                   # year = year,
                                   waterlvl = waterlvl_center,
                                   landbridge = waterlvl_break$landbridge,
                                   min_temp = min_temp_center))

n_model$initializeInfo()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Configure model, build MCMC, compile ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
n_config <- configureMCMC(n_model, print = TRUE, useConjugacy = FALSE,)
n_config$addMonitors("N.est")

target <- c("beta",
            "beta.DD1",
            "beta.DD4",
            "beta.water1",
            "beta.temp0",
            "beta.temp.pow0",
            "beta.bridge1")

n_config$removeSamplers(target = target)
n_config$addSampler(target = target, type = 'AF_slice') 

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
MCMCvis::MCMCdiag(n_samples$samples,
                  file_name = "output/global-model_prediction-capability.txt")

MCMCvis::MCMCtrace(n_samples$samples,
                   params = target,
                   iter = niter,
                   filename =
                     "output/global-model_prediction-capability.pdf")

write_rds(awpe.ssm.nimble, "output/global-model_prediction-capability.RDS")

write_rds(n_samples$samples,
          "output/posteriors/posterior-samples_global-model_prediction-capability.RDS")
