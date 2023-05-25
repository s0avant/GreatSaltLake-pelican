#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                             Data Cleaning Code                               #
#                                 05/19/2023                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(tidyverse)
library(lubridate)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Read in data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
waterlvl <- read.csv("data/water-levels_raw.csv", sep = "")
PDO <- read_csv("data/PDO_raw.csv")
SOI <- read_csv("data/SOI_raw.csv")
min_temp <- read.csv("data/min-temp_raw.csv", header = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Clean data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~ Reformat dates ----
PDO <- PDO %>% mutate(date_ymd = ym(Date))
SOI <- SOI %>% mutate(date_ymd = ym(Date))

#~ Calculate annual means from monthly data ----
PDO <- PDO %>%
  mutate(year = year(date_ymd)) %>%
  group_by(year) %>%
  mutate(mean_PDO = mean(Value))

SOI <- SOI %>%
  mutate(year = year(date_ymd)) %>%
  group_by(year) %>%
  mutate(mean_SOI = mean(Value))
  
min_temp_metric <- min_temp %>%
  mutate(mean = rowMeans(.[, c("Apr", "May", "Jun", "Jul")])) %>%
  select(year = Year, mean) %>%
  mutate(mean = (mean-32)*(5/9))

#~ Collapse data to annual format for study period ----
annual_PDO <- PDO %>%
  select(year, mean_PDO) %>%
  filter(year %in% 1963:2020) %>%
  distinct()

annual_SOI <- SOI %>%
  select(year, mean_SOI) %>%
  filter(year %in% 1963:2020) %>%
  distinct()

#~ Create binary landbridge presence/absence data
waterlvl <- waterlvl %>%
  select(year_nu, month_nu, mean_va)

waterlvl_fall <- waterlvl %>%
  filter(month_nu %in% c(9:12)) %>%
  group_by(year_nu) %>%
  mutate(min_fall = min(mean_va)) %>%
  select(year_nu, min_fall) %>%
  unique() %>%
  ungroup() %>%
  select(min_fall) %>%
  slice(-(nrow(.)-1))

# Adjust years assigned for fall water level, because this won't have an effect
## on breeding AWPE/chicks until the next season
waterlvl_fall <- rbind(NA, waterlvl_fall)

waterlvl_breeding <- waterlvl %>%
  filter(month_nu %in% c(4:7)) %>%
  group_by(year_nu) %>%
  mutate(mean_breeding = mean(mean_va),
         min_breeding = min(mean_va)) %>%
  select(year_nu, mean_breeding, min_breeding) %>%
  unique()

waterlvl_season <-  waterlvl_breeding %>%
  cbind(waterlvl_fall) %>%
  mutate(landbridge_breeding =
           case_when(min_breeding > 4194.5 ~ 0, # As per J. Neill pers comm
                     min_breeding <= 4194.5 ~ 1))  %>%
  mutate(landbridge_fall =
           case_when(min_fall > 4194.5 ~ 0,
                     min_fall <= 4194.5 ~ 1))

waterlvl_season <- waterlvl_season %>%
  select(year_nu, mean_breeding, landbridge_breeding)

waterlvl_metric <- waterlvl_season %>%
  mutate(mean_m = mean_breeding * 0.3048)

colnames(waterlvl_metric) <- c("year", "mean_ft", "landbridge", "mean")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Save cleaned data objects ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
write_rds(waterlvl_metric, "data/waterlevels.RDS")
write_rds(annual_PDO, "data/PDO_annual.RDS")
write_rds(annual_SOI, "data/SOI_annual.RDS")
write_rds(min_temp_metric, "data/min_temp.RDS")
