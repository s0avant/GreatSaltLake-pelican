#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                    Gunnison Colony Counts w/post-hoc PVA                     #
#                                 05/19/2023                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in data & minor data cleaning ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
env_data <- read_csv("data/model-data_normalizedonly-for-PVA.csv") %>%
  filter(year %in% 1985:2040)

global_estimates <-
  read_rds("output/model-output_post-break_full-without-year_for-PVA.RDS")

posteriors <-
  read_rds("output/posteriors/posterior-samples_post-break_full-without-year_for-PVA.RDS")

#~ Last minute data cleaning ----
years_to_estimate <- nrow(env_data) - 36 # n years to predict in the future

lag <- 4 # max n years of lagged effects

env_data_future <- env_data %>%
  filter(year > 2020 - lag)

beta_subset <- global_estimates %>%
  select(mean, "2.5%", "97.5%") %>%
  filter(grepl("beta", rownames(.)))

posteriors_allchains <- rbind(posteriors$chain1,
                              posteriors$chain2, 
                              posteriors$chain3)

#~ Extract all betas from each chain of posterior object ----
betas_chain <- posteriors_allchains %>%
  as.data.frame() %>%
  select(contains("beta"))

colnames(betas_chain) <- rownames(beta_subset)

sigma_proc_chain <- posteriors_allchains[,"sigma.proc"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create additional objects for model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
year <- env_data_future$year

iterations <- nrow(betas_chain)

#~ Generate multiple columns of inits because we have 2 DD covariates
logN.est_init <- cbind(posteriors_allchains[, "logN.est[33]"],
                       posteriors_allchains[, "logN.est[34]"],
                       posteriors_allchains[, "logN.est[35]"],
                       posteriors_allchains[, "logN.est[36]"])


logN.est <- cbind(logN.est_init,
                  matrix(NA, nrow = iterations, ncol = years_to_estimate))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Run PVA loops ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for(t in 1:years_to_estimate){
  logN.est[, t + 4] <-
    betas_chain$beta +
    
    (1 + betas_chain$beta.DD1) * logN.est[, t + 3] +
    betas_chain$beta.DD4 * logN.est[, t] +
    
    betas_chain$beta.water1 * env_data_future$water_mean[t + 3] +
    
    betas_chain$beta.bridge1 * env_data_future$landbridge_absent[t + 3] +
    
    betas_chain$beta.temp0 * env_data_future$temp_predicted[t + 4] +
    betas_chain$beta.temp.pow0 * (env_data_future$temp_predicted[t + 4])^2 +
    
    rnorm(n = iterations, mean = 0, sd = sigma_proc_chain)
}

# Model name for easy saving
model_name <-
  "all_predicted"
  # "waterlvl-predicted_landbridge-predicted_temp-mean"
  # "waterlvl-predicted_landbridge-absent_temp-mean"
# "waterlvl-predicted_landbridge-absent_temp-predicted"
# "waterlvl-mean_landbridge-absent_temp-predicted"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate objects for full dataframe ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
logN.est_df <- logN.est[, -c(1:4)] # Removes init columns

logN.est_posteriors <- posteriors_allchains[, grepl("logN.est", colnames(posteriors_allchains))]

logN.est_posteriors_full <- cbind(logN.est_posteriors,
                                  logN.est_df)

colnames(logN.est_posteriors_full) <-
  c(colnames(logN.est_posteriors), paste0("logN.est[", 37:ncol(logN.est_posteriors_full), "]"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate objects for simplified dataframe ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~ Generate credible intervals ----
logN_est <- data.frame(N = apply(logN.est_df, 2, mean), 
                       LCI = apply(logN.est_df, 2, function(x) quantile(x, 0.025)), 
                       UCI = apply(logN.est_df, 2, function(x) quantile(x, 0.975)))


#~ Create objects from historical population estimates ----
awpe_mean <- global_estimates[grepl("logN.est", rownames(global_estimates)),
                              "mean"]

awpe_LCI <- global_estimates[grepl("logN.est", rownames(global_estimates)),
                             "2.5%"]

awpe_UCI <- global_estimates[grepl("logN.est", rownames(global_estimates)),
                             "97.5%"]

df <- data.frame(year = 1985:2040,
                 mean = exp(c(awpe_mean, logN_est$N)),
                 LCI = exp(c(awpe_LCI, logN_est$LCI)),
                 UCI = exp(c(awpe_UCI, logN_est$UCI)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot population estimates ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Summary of all chains
allchains <- ggplot(df) +
  geom_ribbon(aes(x = year, ymin = LCI, ymax = UCI), fill = "grey80") +
  geom_line(aes(x = year, y = mean)) +
  geom_point(aes(x = year, y = mean))

allchains

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate % change from 2020 to 2040 for each scenario ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
start_end_obs <- df %>%
  filter(year %in% c(1985, 2020)) %>%
  select(mean)

start_end_pred <- df %>%
  filter(year %in% c(2020, 2040)) %>%
  select(mean)

df$percent_change_annual <- NA

for(i in 2:length(df$year)){
  df$percent_change_annual[i] <- ((df$mean[i] / df$mean[i-1])^(1/1) - 1) * 100
}

percent_change_obs <- (start_end_obs[2,] - start_end_obs[1,]) /
  start_end_obs[1,] * 100

percent_change_pred <- (start_end_pred[2,] - start_end_pred[1,]) /
  start_end_pred[1,] * 100

df2 <- rbind(df, c("percent change (1985 - 2020)", percent_change_obs, " ",
                   "percent change (2020 - 2040)", percent_change_pred))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save outputs ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ggsave(allchains +
         ggtitle("Gunnison Island Pelican Population Estimates",
                 subtitle = model_name) +
         theme_classic(),
       width = 10,
       height = 10,
       filename = 
         paste0("output/plots/PVA_", model_name, ".png"))

write_csv(df2,
          file = 
            paste0("output/PVA/PVA_", model_name, ".csv"))

write_rds(logN.est_posteriors_full,
          file = 
            paste0("output/PVA/PVA-posteriors_", model_name, ".RDS"))

