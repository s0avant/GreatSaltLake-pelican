#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                Data Cleaning and Prediction Code For PVA Models              #
#                                 05/19/2023                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
awpe <- read.csv("data/colony_count_simple.csv", na.strings = "", header = TRUE)
waterlvl <- readRDS("data/waterlevels.RDS")
min_temp <- readRDS("data/min_temp.RDS")

#~ Last minute data cleaning ----
years_to_estimate <- 20

future <- data.frame(year = 2021:(2020 + years_to_estimate),
                     count = rep(NA, years_to_estimate))
awpe <- rbind(awpe, future)
rownames(awpe) <- NULL

year <- seq(1985, 2020)

#~ Filter to 1985 breakpoint data ----
counts <- filter(awpe, year %in% 1985:2020)
waterlvl_break <- filter(waterlvl, year %in% 1985:2020) %>% pull(mean)
waterlvlft_break <- filter(waterlvl, year %in% 1985:2020) %>% pull(mean_ft)
min_temp_break <- filter(min_temp, year %in% 1985:2020) %>% pull(mean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Predict future environmental data values based on current data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
newdata <- data.frame(year = seq(2021, 2020 + years_to_estimate))

#~ Water levels & landbridge ----
waterlvl.lm <- lm(waterlvl_break ~ year) # lm of observed data
water_future <- as.numeric(predict(waterlvl.lm, newdata, type = "response"))
plot(c(waterlvl_break, water_future))

water_predicted <- c(waterlvl_break, water_future)
water_mean <- c(waterlvl_break, rep(mean(waterlvl_break), nrow(newdata)))
water_min <- c(waterlvl_break, rep(min(waterlvl_break), nrow(newdata)))
water_max <- c(waterlvl_break, rep(max(waterlvl_break), nrow(newdata)))

#~ Use imperial units to calculate land bridge, because that's how it was done originally
waterlvlft.lm <- lm(waterlvlft_break ~ year) # lm of observed data
waterft_future <- as.numeric(predict(waterlvlft.lm, newdata, type = "response"))
waterft_predicted <- c(waterlvlft_break, waterft_future)
plot(c(waterlvlft_break, waterft_future))

landbridge_predicted <- case_when(waterft_predicted > 4194.983 ~ 0,
                                  waterft_predicted <= 4194.983 ~ 1)
landbridge_break <- case_when(waterlvlft_break > 4194.983 ~ 0,
                              waterlvlft_break <= 4194.983 ~ 1)

#~ Normalize variables
c.water_predicted <- (water_predicted - mean(water_predicted))/
  sd(water_predicted)
c.water_mean <- (water_mean - mean(water_mean))/sd(water_mean)
c.water_min <- (water_min - mean(water_min))/sd(water_min)
c.water_max <- (water_max - mean(water_max))/sd(water_max)

par(mfrow = c(2,2))
plot(c.water_predicted)
plot(c.water_mean)
plot(c.water_min)
plot(c.water_max)
par(mfrow = c(1,1))

#~ Min temperature ----
temp.lm <- lm(min_temp_break ~ year) # lm of observed data
temp_future <- as.numeric(predict(temp.lm, newdata, type = "response"))
plot(c(min_temp_break, temp_future))

temp_predicted <- c(min_temp_break, temp_future)
temp_mean <- c(min_temp_break, rep(mean(min_temp_break), nrow(newdata)))
temp_min <- c(min_temp_break, rep(min(min_temp_break), nrow(newdata)))
temp_max <- c(min_temp_break, rep(max(min_temp_break), nrow(newdata)))

#~ Normalize variables ----
c.temp_predicted <- (temp_predicted - mean(temp_predicted))/
  sd(temp_predicted)
c.temp_mean <- (temp_mean - mean(temp_mean))/sd(temp_mean)
c.temp_min <- (temp_min - mean(temp_min))/sd(temp_min)
c.temp_max <- (temp_max - mean(temp_max))/sd(temp_max)

par(mfrow = c(2,2))
plot(c.temp_predicted)
plot(c.temp_mean)
plot(c.temp_min)
plot(c.temp_max)
par(mfrow = c(1,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build data frames for PVA ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
waterlvl_df <- data.frame(water_predicted = c.water_predicted,
                          water_mean = c.water_mean,
                          water_min = c.water_min,
                          water_max = c.water_max)

landbridge_df <- data.frame(landbridge_predicted = landbridge_predicted,
                            landbridge_present =
                              c(landbridge_break, rep(1, years_to_estimate)),
                            landbridge_absent =
                              c(landbridge_break, rep(0, years_to_estimate)))

temp_df <- data.frame(temp_predicted = c.temp_predicted,
                      temp_mean = c.temp_mean,
                      temp_min = c.temp_min,
                      temp_max = c.temp_max)

year_scaled <- (c(year, newdata$year) - min(c(year, newdata$year)))/100

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create input objects for predicted, min, max, & mean vars ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
year_df <- data.frame(year = seq(min(year), max(newdata), 1),
                      year_scaled = year_scaled)

env_vars_df <- cbind(year_df, waterlvl_df, landbridge_df, temp_df)

env_vars_df

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save output ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
write_csv(env_vars_df, "data/model-data_normalizedonly-for-PVA.csv")
