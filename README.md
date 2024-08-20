# GreatSaltLake-pelican
Public repository for data and code peer review of Great Salt Lake pelican population drivers.

## Data file description
"colony_count_simple.csv"
- Gunnison Island annual American white pelican *(Pelecanus erythrorhynchos)* colony count dataset
- Source: Great Salt Lake Ecosystem Program - Utah Division of Wildlife Resources

"min-temp_raw.csv"
- Monthly minimum temperature data (ºF) collected at the Salt Lake City International Airport weather station
- Retrieved from the National Weather Service [https://w2.weather.gov/climate/xmacis.pp?wfo=slc]

"PDO_raw.csv"
- Monthly Pacific Decadal Oscillation (PDO) index data
- Retrieved from the National Centers for Environmental Information [https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat]
"SOI_raw.csv"

- Monthly Southern Oscillation Index (SOI) data
- Retrieved from the NWS Climate Prediction Center [https://www.cpc.ncep.noaa.gov/data/indices/soi]
"water-levels_raw.csv"

- Monthly water level data (feet) collected at the USGS water level gauge nearest to Gunnison Island (USGS Saline, UT)
- Retrieved from USGS: Water Data for the Nation [https://waterdata.usgs.gov/monitoring-location/10010100/]

## Code file description
"01a_clean-data_nimble.R"
- Script to prepare raw data for Nimble models

"01b_clean-data_PVA.R"
- Script to project environmental data into the future for use in Population Viability Analyisis (PVA)

"02a_univariate-screening_DD.R"
- Script that runs univariate density dependence Nimble model with time lags

"02b_univariate-screening_landbridge.R"
- Script that runs univariate land bridge presence/absence Nimble model with time lags

"02c_univariate-screening_min-temp.R"
- Script that runs univariate mean minimum spring temperature (April-July) Nimble model with time lags

"02d_univariate-screening_PDO.R"
- Script that runs univariate Pacific Decadal Oscillation Nimble model with time lags

"02e_univariate-screening_SOI.R"
- Script that runs univariate Southern Oscillation Index Nimble model with time lags

"02f_univariate-screening_waterlevel.R"
- Script that runs univariate Great Salt Lake water level Nimble model with time lags

"03_model-selection_C180.R"
- Script to calculate 80% credible intervals of each covariate within univariate screening models

"04_global-model.R"
- Script that runs final "global" model with only covariates with 80% credible intervals not overlapping zero retained

"05_PVA-model.R"
- Script to run PVA Nimble model

"06_PVA-forecast.R"
- Post-hoc population projection under selected management and environmental scenarios

"07_prediction-capability.R"
- Script assessing model fit and prediction capability

"08_prediction-intervals.R"
- Script to generate prediction intervals
