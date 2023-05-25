#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                        Model Selection (80% CI) Code                         #
#                                 05/19/2023                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(tidyverse)
library(bayestestR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# List model files & build objects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
model_names <- c("Density_Dependence",
                 "PDO",
                 "SOI",
                 "Water_Levels",
                 "Minimum_Spring_Temperature",
                 "Landbridge_Presence")

model_summary <- c("output/model-output_post-break_lag_DD.RDS",
                   "output/model-output_post-break_lag_PDO.RDS",
                   "output/model-output_post-break_lag_SOI.RDS",
                   "output/model-output_post-break_lag_waterlvl.RDS",
                   "output/model-output_post-break_lag_min-temp.RDS",
                   "output/model-output_post-break_lag_landbridge.RDS")

posterior_samples <- c("output/posteriors/posterior-samples_post-break_lag_DD.RDS",
                       "output/posteriors/posterior-samples_post-break_lag_PDO.RDS",
                       "output/posteriors/posterior-samples_post-break_lag_SOI.RDS",
                       "output/posteriors/posterior-samples_post-break_lag_waterlvl.RDS",
                       "output/posteriors/posterior-samples_post-break_lag_min-temp.RDS",
                       "output/posteriors/posterior-samples_post-break_lag_landbridge.RDS")

CI80 <- data.frame(model = NA,
                   beta = NA,
                   mean = NA,
                   low = NA,
                   high = NA,
                   overlap_zero = NA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Begin for loop to calculate 80% CI & generate plots ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
for(i in 1:length(posterior_samples)){
  # Read in model files
  awpe.ssm.nimble <- read_rds(model_summary[i])
  n_samples <- read_rds(posterior_samples[i])
  
  posterior_beta <- n_samples %>%
    as.data.frame() %>%
    select(contains("beta"))
  
  df_CI <- as.data.frame(matrix(NA, nrow = 2, ncol = ncol(posterior_beta)))
  rownames(df_CI) <- c("CI_low", "CI_high")
  colnames(df_CI) <- colnames(posterior_beta)
  
  for(j in 1:ncol(df_CI)){
    df_CI[1, j] <- quantile(posterior_beta[, j], 0.1)
    df_CI[2, j] <- quantile(posterior_beta[, j], 0.9)
  }
  
  beta <- n_samples$chain1 %>%
    as.data.frame() %>%
    select(contains("beta")) %>%
    colnames()
  beta <- beta[2:length(beta)] # Keeps covariates only
  
  for(k in beta){
    cred_int <- df_CI %>%
      select(contains(k)) %>%
      rowMeans() %>%
      as.numeric()
    
    posterior_mean <- awpe.ssm.nimble %>%
      filter(rownames(awpe.ssm.nimble) %in% k) %>%
      select("mean")
    
    density <- posterior_beta %>%
      select(contains(k)) %>%
      estimate_density(extend = TRUE) %>%
      ggplot(aes(x = x, y = y)) + 
      geom_area(aes(fill = Parameter), alpha = 0.4) +
      theme_classic() +
      geom_vline(xintercept = cred_int[1], color = "royalblue") +
      geom_vline(xintercept = cred_int[2], color = "royalblue") +
      ggtitle(paste0(model_names[i], " ", k))
    
    temp <- data.frame(model = model_names[i],
                       beta = k,
                       mean = posterior_mean,
                       low = cred_int[1],
                       high = cred_int[2],
                       overlap_zero = ifelse(cred_int[1] > 0 & cred_int[2] > 0 |
                                               cred_int[1] < 0 & cred_int[2] < 0,
                                             "no", "yes"
                       ))
    rownames(temp) <- NULL
    
    CI80 <- rbind(CI80, temp)
    
    # ggsave(paste0("output/plots/posterior_CI80_",
    #               model_names[i], "_", k, ".pdf"),
    #        width = 8, height = 6)
    
  }
}

#~ Extract all posteriors that don't overlap zero ----
CI80 <- CI80[-1, ]

sig <- CI80 %>%
  filter(overlap_zero == "no" & beta != "beta.yr")

sig

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save files ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
write_csv(CI80, "output/CI80_all.csv")
write_csv(sig, "output/CI80_no-zero-overlap.csv")
