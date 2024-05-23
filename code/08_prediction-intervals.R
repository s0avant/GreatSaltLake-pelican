#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                       PVA Prediction Interval Creation                       #
#                                 02/02/2024                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in data & check structure ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
awpe_ssm <-
  readRDS("output/model-output_post-break_full-without-year_for-PVA.RDS")
prediction <-
  readRDS("output/global-model_prediction-capability.RDS")

awpe_ssm_posterior <-
  readRDS("output/posteriors/posterior-samples_post-break_full-without-year_for-PVA.RDS")
prediction_posterior <-
  readRDS("output/posteriors/posterior-samples_global-model_prediction-capability.RDS")

str(awpe_ssm_posterior)
str(prediction_posterior)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean data objects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
prediction_posterior_allchains <- rbind(prediction_posterior$chain1,
                                        prediction_posterior$chain2, 
                                        prediction_posterior$chain3)

awpe_ssm_allchains <- rbind(awpe_ssm_posterior$chain1[1:90000,],
                            awpe_ssm_posterior$chain2[1:90000,],
                            awpe_ssm_posterior$chain3[1:90000,])

str(prediction_posterior_allchains)
str(awpe_ssm_allchains)

prediction_sigma_df <- cbind(prediction_posterior_allchains,
                             sigma.obs_ssm = awpe_ssm_allchains[,"sigma.obs"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate prediction interval ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
iterations <- nrow(prediction_posterior_allchains)

# Select only logN.est posterior estimates
logN.est <- prediction_sigma_df %>%
  as.data.frame() %>%
  select(starts_with(c("logN.est"))) %>%
  colnames()

prediction_int <- matrix(NA, nrow = iterations, ncol = length(logN.est))

# Take random draw from normal dist. 
for(i in 1:length(logN.est)){
  prediction_int[, i] <- rnorm(n = iterations, 
                               mean = prediction_sigma_df[,logN.est[i]],
                               sd = prediction_sigma_df[,"sigma.obs_ssm"])
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
lpi_upi <- data.frame(LPI = apply(prediction_int, 2, function(x) quantile(x, 0.025)), 
                      UPI = apply(prediction_int, 2, function(x) quantile(x, 0.975)))

df <- data.frame(year = 1985:2020,
                 counts = awpe_ssm$mean[grepl("logN.est", rownames(awpe_ssm))],
                 counts_LCI = awpe_ssm$`2.5%`[grepl("logN.est", rownames(awpe_ssm))],
                 counts_UCI = awpe_ssm$`97.5%`[grepl("logN.est", rownames(awpe_ssm))],
                 predicted = prediction$mean[grepl("logN.est", rownames(prediction))],
                 pred_LCI = prediction$`2.5%`[grepl("logN.est", rownames(prediction))],
                 pred_UCI = prediction$`97.5%`[grepl("logN.est", rownames(prediction))],
                 pred_LPI = lpi_upi$LPI,
                 pred_UPI = lpi_upi$UPI)


pi_plot_exp <- df %>%
  mutate(counts = exp(counts),
         counts_LCI = exp(counts_LCI),
         counts_UCI = exp(counts_UCI),
         predicted = exp(predicted), 
         pred_LCI = exp(pred_LCI),
         pred_UCI = exp(pred_UCI),
         pred_LPI = exp(pred_LPI),
         pred_UPI = exp(pred_UPI)) %>%
  
  ggplot() +
  geom_ribbon(aes(x = year,
                  ymin = pred_LCI,
                  ymax = pred_UCI,
                  fill = "Training data credible"), alpha = 0.3) +
  geom_ribbon(aes(x = year,
                  ymin = pred_LPI,
                  ymax = pred_UPI,
                  fill = "Prediction"), alpha = 0.2) +
  geom_ribbon(aes(x = year,
                  ymin = counts_LCI,
                  ymax = counts_UCI,
                  fill = "Full time-series credible"), alpha = 0.4) +
  geom_line(aes(x = year, y = predicted, color = "Training data")) +
  geom_line(aes(x = year, y = counts, color = "Full time-series")) +
  geom_vline(aes(xintercept = 2013), linetype = "dashed", alpha = 0.5) +
  geom_vline(aes(xintercept = 2020), linetype = "dashed", alpha = 0.5) +
  xlab("Year") + ylab("Abundance") +
  annotate("text", x = 1998.5, y = 1000, label = "Training data", alpha = 0.6) +
  annotate("text", x = 2016.5, y = 1000, label = "Hold-out data", alpha = 0.6) +
  ggtitle("Global state space model prediction performance",
          subtitle = "Predictions made with full time-series versus 80% training data") +
  labs(fill = "95% interval", color = "Mean estimates") +
  theme_classic() +
  scale_color_manual(values = c("orange3", "grey40")) +
  scale_fill_manual(values = c("orange2", "grey40", "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = c(.14, .2),
        legend.background = element_rect(fill = "transparent"))

pi_plot_log <- df %>%
  ggplot() +
  geom_ribbon(aes(x = year,
                  ymin = pred_LCI,
                  ymax = pred_UCI,
                  fill = "Training data credible"), alpha = 0.3) +
  geom_ribbon(aes(x = year,
                  ymin = pred_LPI,
                  ymax = pred_UPI,
                  fill = "Prediction"), alpha = 0.2) +
  geom_ribbon(aes(x = year,
                  ymin = counts_LCI,
                  ymax = counts_UCI,
                  fill = "Full time-series credible"), alpha = 0.4) +
  geom_line(aes(x = year, y = predicted, color = "Training data")) +
  geom_line(aes(x = year, y = counts, color = "Full time-series")) +
  geom_vline(aes(xintercept = 2013), linetype = "dashed", alpha = 0.5) +
  geom_vline(aes(xintercept = 2020), linetype = "dashed", alpha = 0.5) +
  xlab("Year") + ylab("log(Abundance)") +
  annotate("text", x = 1998.5, y = 7.5, label = "Training data", alpha = 0.6) +
  annotate("text", x = 2016.5, y = 7.5, label = "Hold-out data", alpha = 0.6) +
  ggtitle("Global state space model prediction performance",
          subtitle = "Predictions made with full time-series versus 80% training data") +
  labs(fill = "95% interval", color = "Mean estimates") +
  theme_classic() +
  scale_color_manual(values = c("orange3", "grey40")) +
  scale_fill_manual(values = c("orange2", "grey40", "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = c(.14, .2),
        legend.background = element_rect(fill = "transparent"))

ggsave(plot = pi_plot_exp, "output/plots/predicted_interval_plot.png",
       height = 8,
       width = 10,
       units = "in")

ggsave(plot = pi_plot_log, "output/plots/log_predicted_interval_plot.png",
       height = 8,
       width = 10,
       units = "in")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calculate RMSE & bias metric ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
df_subset <- df %>%
  filter(year >= 2013)

rmse_rmspe <- data.frame(rmse = sqrt(mean((df_subset$counts -
                                             df_subset$predicted)^2)),
                         rmspe = sqrt(mean(100 * abs(abs(df_subset$counts -
                                                           df_subset$predicted)/df_subset$counts))^2))

write_csv(rmse_rmspe, "output/rmse_rmspe.csv")

pred_obs <- (df_subset$predicted - df_subset$counts)/df_subset$counts
mean_pred_obs <- mean(pred_obs)
df2 <- data.frame(year = c(2013:2020),
                  pred_obs = pred_obs)

write_csv(df2, "output/bias_predicted.csv")

bias_plot <- ggplot(df2) +
  geom_point(aes(x = year, y = pred_obs), shape = 1) +
  geom_hline(aes(yintercept = mean_pred_obs), linetype = "dashed") +
  theme_light() +
  ggtitle("Distance: predicted from observed",
          subtitle = "'(predicted - observed)/observed'") +
  xlab("Year") + ylab("Value") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(plot = bias_plot, "output/plots/bias_plot.png",
       height = 5,
       width = 5,
       units = "in")
