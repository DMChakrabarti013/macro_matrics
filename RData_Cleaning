# D.M. Chakrabarti
# March 29, 2025
# Project on Unconventional Monetary Policy
# This file - data import, merging, cleaning, organization, summary

###### Packages ########
if (!require(readxl)) install.packages("readxl")
if(!require(readr)) install.packages("readr")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(tidyr)) install.packages("tidyr")
if(!require(stringr)) install.packages("stringr")
if(!require(purrr)) install.packages("purrr")
if(!require(mFilter)) install.packages("mFilter")
if(!require(lubridate)) install.packages("lubridate")
if(!require(zoo)) install.packages("zoo")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(moments)) install.packages("moments")
if(!require(tseries)) install.packages("tseries")
if(!require(xtable)) install.packages("xtable")
if (!require(urca)) install.packages("urca")
if (!require(readxl)) install.packages("readxl")
if (!require(readr)) install.packages("readr")
if (!require(openxlsx)) install.packages("openxlsx")

library(readr)
library(dplyr)
library(readxl)
library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)
library(mFilter)
library(lubridate)
library(zoo)
library(ggplot2)
library(moments)
library(tseries)
library(xtable)
library(urca)
library(readxl)
library(readr)
library(openxlsx)


####### reshaping and importing GDP data #########


# function to import and reshape a single GDP file
import_state_gdp <- function(file_path) {
  # Read the CSV file while skipping the first 3 rows so row 4 becomes the header.
  df <- read_csv(file_path, skip = 3, show_col_types = FALSE)
  
  # Convert columns 5 to end to numeric.
  # We assume the first 4 columns are text (GeoFips, GeoName, LineCode, Description)
  df[ , 5:ncol(df)] <- lapply(df[ , 5:ncol(df)], function(x) as.numeric(x))
  
  # Check for any conversion issues; if many NAs appear, use problems(df) to inspect.
  if(any(is.na(df[ , 5:ncol(df)]))) {
    warning("Some columns from 5 onward contain NA values after conversion. Check the parsing issues with problems(df).")
  }
  
  # Filter for the row corresponding to the series of interest ("All industry total")
  df_total <- df %>% filter(Description == "All industry total")
  
  # Pivot the time period columns (all columns except the first four) to long format.
  df_long <- df_total %>%
    pivot_longer(
      cols = -c(GeoFips, GeoName, LineCode, Description),
      names_to = "Time",
      values_to = "GDP"
    )
  
  # Extract state name from the file name: remove the fixed suffix, e.g., "_quarterly_gdp.csv"
  state_name <- str_remove(basename(file_path), "_quarterly_gdp\\.csv")
  df_long <- df_long %>% mutate(State = state_name)
  
  return(df_long)
}

gdp_files <- list.files(path = "/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/gdp_real_2017",
                        pattern = "\\.csv$", full.names = TRUE)

all_gdp_data <- map_dfr(gdp_files, import_state_gdp)

glimpse(all_gdp_data)
head(all_gdp_data)

# converting Time to quarterly format for R
all_gdp_data <- all_gdp_data %>%
  mutate(YearQuarter = as.yearqtr(Time, format = "%Y:Q%q"))

head(all_gdp_data)

# hp filter for cyclical component
all_gdp_data_cycle <- all_gdp_data %>%
  group_by(GeoName) %>%
  arrange(YearQuarter) %>%
  group_modify(~ {
    hp_out <- hpfilter(.x$GDP, freq = 1600) # quarterly data
    .x %>% mutate(cyclicalGDP = hp_out$cycle)
  })


# summary statistics for each state
state_gdp_stats <- all_gdp_data_cycle %>%
  group_by(GeoName) %>%
  summarize(
    Mean = mean(cyclicalGDP, na.rm = TRUE),
    Median = median(cyclicalGDP, na.rm = TRUE),
    SD = sd(cyclicalGDP, na.rm = TRUE),
    Min = min(cyclicalGDP, na.rm = TRUE),
    Max = max(cyclicalGDP, na.rm = TRUE),
    IQR = IQR(cyclicalGDP, na.rm = TRUE),
    Skewness = skewness(cyclicalGDP, na.rm = TRUE),
    Kurtosis = kurtosis(cyclicalGDP, na.rm = TRUE),
    ACF1 = if(n() > 1) cor(cyclicalGDP[2:n()], cyclicalGDP[1:(n()-1)], use = "complete.obs") else NA,
    PP_p = tryCatch(pp.test(cyclicalGDP$p.value, error = function(e) NA),
    ZA_stat = tryCatch({
      za_test <- ur.za(cyclicalGDP, model = "both", lag = 1)
      za_test@teststat
    }, error = function(e) NA)
  ) %>%
  ungroup()

state_gdp_stats_table <- xtable(state_gdp_stats)
print(state_gdp_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

p_mean <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, Mean), y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Mean of Cyclical GDP by State", x = "State", y = "Mean")
print(p_mean)

p_sd <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, SD), y = SD)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Standard Deviation of Cyclical GDP by State", x = "State", y = "SD")
print(p_sd)

p_skew <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, Skewness), y = Skewness)) +
  geom_point(color = "purple", size = 3) +
  coord_flip() +
  labs(title = "Skewness of Cyclical GDP by State", x = "State", y = "Skewness")
print(p_skew)

p_kurt <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, Kurtosis), y = Kurtosis)) +
  geom_point(color = "orange", size = 3) +
  coord_flip() +
  labs(title = "Kurtosis of Cyclical GDP by State", x = "State", y = "Kurtosis")
print(p_kurt)

p_acf <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, ACF1), y = ACF1)) +
  geom_point(color = "brown", size = 3) +
  coord_flip() +
  labs(title = "Lag-1 Autocorrelation of Cyclical GDP by State", x = "State", y = "ACF(1)")
print(p_acf)

p_PP <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, PP_p), y = PP_p)) +
  geom_point(color = "red", size = 3) +
  coord_flip() +
  labs(title = "PP Test p-values for Cyclical GDP by State", x = "State", y = "PP p-value") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  annotate("text", x = 1, y = 0.06, label = "0.05 threshold", hjust = 0)
print(p_PP)

p_ZA <- ggplot(state_gdp_stats, aes(x = reorder(GeoName, ZA_stat), y = ZA_stat)) +
  geom_point(color = "blue", size = 3) +
  coord_flip() +
  labs(title = "Zivot-Andrews Test Statistic for Cyclical GDP by State", x = "State", y = "ZA Statistic")
print(p_ZA)


############# importing and exploring interest rates data #########

composite_rates_data <- read_csv("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/interest_rates/composite_rates.csv")

composite_rates_data <- composite_rates_data %>%
  mutate(Date = as.Date(as.yearmon(month_time, '%b-%y'))) %>%
  arrange(Date)

composite_rates_data <- composite_rates_data %>%
  mutate(composite_rates_diff = c(NA, diff(composite_rates)))

# summary statistics
interest_rate_stats <- composite_rates_data %>%
  summarize(
    Mean = mean(composite_rates, na.rm = TRUE),
    Median = median(composite_rates, na.rm = TRUE),
    SD = sd(composite_rates, na.rm = TRUE),
    Min = min(composite_rates, na.rm = TRUE),
    Max = max(composite_rates, na.rm = TRUE),
    IQR = IQR(composite_rates, na.rm = TRUE),
    Skewness = skewness(composite_rates, na.rm = TRUE),
    Kurtosis = kurtosis(composite_rates, na.rm = TRUE),
    ACF1_diff = if(sum(!is.na(composite_rates)) > 1) {
      cor(composite_rates[!is.na(composite_rates)][-1],
          composite_rates[!is.na(composite_rates)][-length(composite_rates[!is.na(composite_rates)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p_diff = tryCatch(pp.test(na.omit(composite_rates))$p.value, error = function(e) NA)
    
  )

interest_rate_stats_table <- xtable(interest_rate_stats)
print(interest_rate_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

p_interest <- ggplot(composite_rates_data, aes(x = composite_rates_data$Date, y = composite_rates_data$composite_rates)) +
  geom_line(color = "blue") +
  labs(title = "Composite Interest Rate", 
       x = "Date", 
       y = "Rate") +
  theme_minimal()
print(p_interest)


######## importing and analyzing financial stress ###############
financial_stress_data <- read_csv("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/financial stress/STLFSI4.csv")

financial_stress_stats <- financial_stress_data %>%
  summarize(
    Mean = mean(STLFSI4, na.rm = TRUE),
    Median = median(STLFSI4, na.rm = TRUE),
    SD = sd(STLFSI4, na.rm = TRUE),
    Min = min(STLFSI4, na.rm = TRUE),
    Max = max(STLFSI4, na.rm = TRUE),
    IQR = IQR(STLFSI4, na.rm = TRUE),
    Skewness = skewness(STLFSI4, na.rm = TRUE),
    Kurtosis = kurtosis(STLFSI4, na.rm = TRUE),
    ACF1 = if(sum(!is.na(STLFSI4)) > 1) {
      cor(STLFSI4[!is.na(STLFSI4)][-1],
          STLFSI4[!is.na(STLFSI4)][-length(STLFSI4[!is.na(STLFSI4)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(STLFSI4))$p.value, error = function(e) NA)
    
  )

financial_stress_stats_table <- xtable(financial_stress_stats)
print(financial_stress_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

######## importing and analyzing VIX ###############

vix_history_data <- read_csv("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/market volatility/VIX_History.csv")
vvix_history_data <- read_csv("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/market volatility/VVIX_History.csv")

vix_history_open_stats <- vix_history_data %>%
  summarize(
    Mean = mean(OPEN, na.rm = TRUE),
    Median = median(OPEN, na.rm = TRUE),
    SD = sd(OPEN, na.rm = TRUE),
    Min = min(OPEN, na.rm = TRUE),
    Max = max(OPEN, na.rm = TRUE),
    IQR = IQR(OPEN, na.rm = TRUE),
    Skewness = skewness(OPEN, na.rm = TRUE),
    Kurtosis = kurtosis(OPEN, na.rm = TRUE),
    ACF1 = if(sum(!is.na(OPEN)) > 1) {
      cor(OPEN[!is.na(OPEN)][-1],
          OPEN[!is.na(OPEN)][-length(OPEN[!is.na(OPEN)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(OPEN))$p.value, error = function(e) NA)
    
  )
vix_history_open_stats_table <- xtable(vix_history_open_stats)
print(vix_history_open_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)


vix_history_high_stats <- vix_history_data %>%
  summarize(
    Mean = mean(HIGH, na.rm = TRUE),
    Median = median(HIGH, na.rm = TRUE),
    SD = sd(HIGH, na.rm = TRUE),
    Min = min(HIGH, na.rm = TRUE),
    Max = max(HIGH, na.rm = TRUE),
    IQR = IQR(HIGH, na.rm = TRUE),
    Skewness = skewness(HIGH, na.rm = TRUE),
    Kurtosis = kurtosis(HIGH, na.rm = TRUE),
    ACF1 = if(sum(!is.na(HIGH)) > 1) {
      cor(HIGH[!is.na(HIGH)][-1],
          HIGH[!is.na(HIGH)][-length(HIGH[!is.na(HIGH)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(HIGH))$p.value, error = function(e) NA)
    
  )
vix_history_high_stats_table <- xtable(vix_history_high_stats)
print(vix_history_high_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)



vix_history_low_stats <- vix_history_data %>%
  summarize(
    Mean = mean(LOW, na.rm = TRUE),
    Median = median(LOW, na.rm = TRUE),
    SD = sd(LOW, na.rm = TRUE),
    Min = min(LOW, na.rm = TRUE),
    Max = max(LOW, na.rm = TRUE),
    IQR = IQR(LOW, na.rm = TRUE),
    Skewness = skewness(LOW, na.rm = TRUE),
    Kurtosis = kurtosis(LOW, na.rm = TRUE),
    ACF1 = if(sum(!is.na(LOW)) > 1) {
      cor(LOW[!is.na(LOW)][-1],
          LOW[!is.na(LOW)][-length(LOW[!is.na(LOW)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(LOW))$p.value, error = function(e) NA)
    
  )
vix_history_low_stats_table <- xtable(vix_history_low_stats)
print(vix_history_low_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

vix_history_close_stats <- vix_history_data %>%
  summarize(
    Mean = mean(CLOSE, na.rm = TRUE),
    Median = median(CLOSE, na.rm = TRUE),
    SD = sd(CLOSE, na.rm = TRUE),
    Min = min(CLOSE, na.rm = TRUE),
    Max = max(CLOSE, na.rm = TRUE),
    IQR = IQR(CLOSE, na.rm = TRUE),
    Skewness = skewness(CLOSE, na.rm = TRUE),
    Kurtosis = kurtosis(CLOSE, na.rm = TRUE),
    ACF1 = if(sum(!is.na(CLOSE)) > 1) {
      cor(CLOSE[!is.na(CLOSE)][-1],
          CLOSE[!is.na(CLOSE)][-length(CLOSE[!is.na(CLOSE)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(CLOSE))$p.value, error = function(e) NA)
    
  )
vix_history_close_stats_table <- xtable(vix_history_close_stats)
print(vix_history_close_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)




vvix_history_stats <- vvix_history_data %>%
  summarize(
    Mean = mean(VVIX, na.rm = TRUE),
    Median = median(VVIX, na.rm = TRUE),
    SD = sd(VVIX, na.rm = TRUE),
    Min = min(VVIX, na.rm = TRUE),
    Max = max(VVIX, na.rm = TRUE),
    IQR = IQR(VVIX, na.rm = TRUE),
    Skewness = skewness(VVIX, na.rm = TRUE),
    Kurtosis = kurtosis(VVIX, na.rm = TRUE),
    ACF1 = if(sum(!is.na(VVIX)) > 1) {
      cor(VVIX[!is.na(VVIX)][-1],
          VVIX[!is.na(VVIX)][-length(VVIX[!is.na(VVIX)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(VVIX))$p.value, error = function(e) NA)
    
  )
vvix_history_stats_table <- xtable(vvix_history_stats)
print(vvix_history_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

####### importing and exploring economic policy uncertainty ###########
policy_index_data <- read_csv("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/uncertainty/USEPUINDXD.csv")
policy_uncertainty_data_1 <- read_excel("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/uncertainty/US_Policy_Uncertainty_Data.xlsx", sheet = 1)
policy_uncertainty_data_2 <- read_excel("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/uncertainty/US_Policy_Uncertainty_Data.xlsx", sheet = 2)

policy_index_stats <- policy_index_data %>%
  summarize(
    Mean = mean(USEPUINDXD, na.rm = TRUE),
    Median = median(USEPUINDXD, na.rm = TRUE),
    SD = sd(USEPUINDXD, na.rm = TRUE),
    Min = min(USEPUINDXD, na.rm = TRUE),
    Max = max(USEPUINDXD, na.rm = TRUE),
    IQR = IQR(USEPUINDXD, na.rm = TRUE),
    Skewness = skewness(USEPUINDXD, na.rm = TRUE),
    Kurtosis = kurtosis(USEPUINDXD, na.rm = TRUE),
    ACF1 = if(sum(!is.na(USEPUINDXD)) > 1) {
      cor(USEPUINDXD[!is.na(USEPUINDXD)][-1],
          USEPUINDXD[!is.na(USEPUINDXD)][-length(USEPUINDXD[!is.na(USEPUINDXD)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(USEPUINDXD))$p.value, error = function(e) NA)
    
  )
policy_index_stats_table <- xtable(policy_index_stats)
print(policy_index_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

# state policy uncertainty
state_policy_uncertainty_data <- read_excel("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/exogenous shocks/uncertainty/State_Policy_Uncertainty.xlsx")

state_policy_uncertainty_long <- state_policy_uncertainty_data %>%
  pivot_longer(
    cols = -c(year,month),
    names_to = c("EPUType", "State"),
    names_pattern = "EPU_([^0-9]+)([A-Z]{2})",
    values_to = "EPU_value"
  )

state_policy_uncertainty_stats <- state_policy_uncertainty_long %>%
  group_by(State, EPUType) %>%
  summarize(
    Mean = mean(EPU_value, na.rm = TRUE),
    Median = median(EPU_value, na.rm = TRUE),
    SD = sd(EPU_value, na.rm = TRUE),
    Min = min(EPU_value, na.rm = TRUE),
    Max = max(EPU_value, na.rm = TRUE),
    Skewness = skewness(EPU_value, na.rm = TRUE),
    Kurtosis = kurtosis(EPU_value, na.rm = TRUE),
    IQR = IQR(EPU_value, na.rm = TRUE),
    ACF1 = if(sum(!is.na(EPU_value)) > 1) {
      cor(EPU_value[!is.na(EPU_value)][-1],
          EPU_value[!is.na(EPU_value)][-length(EPU_value[!is.na(EPU_value)])],
          use = "complete.obs")
    } else NA,
    PP_p = tryCatch(pp.test(na.omit(EPU_value))$p.value, error = function(e) NA),
    .groups = "drop"
  )
state_policy_uncertainty_stats_table <- xtable(state_policy_uncertainty_stats)
print(state_policy_uncertainty_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

######## importing and exploring inflaiton data #########
inflation_data <- read_csv("/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/time varying covariates/inflation/CPIAUCSL.csv")

inflation_stats <- inflation_data %>%
  summarize(
    Mean = mean(CPIAUCSL, na.rm = TRUE),
    Median = median(CPIAUCSL, na.rm = TRUE),
    SD = sd(CPIAUCSL, na.rm = TRUE),
    Min = min(CPIAUCSL, na.rm = TRUE),
    Max = max(CPIAUCSL, na.rm = TRUE),
    IQR = IQR(CPIAUCSL, na.rm = TRUE),
    Skewness = skewness(CPIAUCSL, na.rm = TRUE),
    Kurtosis = kurtosis(CPIAUCSL, na.rm = TRUE),
    ACF1 = if(sum(!is.na(CPIAUCSL)) > 1) {
      cor(CPIAUCSL[!is.na(CPIAUCSL)][-1],
          CPIAUCSL[!is.na(CPIAUCSL)][-length(CPIAUCSL[!is.na(CPIAUCSL)])],
          use = "complete.obs")
    } else NA,
    # Compute the PP test p-value for the differenced series
    PP_p = tryCatch(pp.test(na.omit(CPIAUCSL))$p.value, error = function(e) NA)
    
  )
inflation_stats_table <- xtable(inflation_stats)
print(inflation_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)


####### import and analyze statewise unemployment data ########
import_employment <- function(file_path) {
  df <- read_excel(file_path, skip = 10)
  
  state_name <- str_remove(basename(file_path), "_employment\\.xlsx")
  
  df <- df %>%
    mutate(State = state_name)
  
  df <- df %>%
    mutate(Date = as.Date(paste(Year, Period, "01", sep = "-"), format = "%Y-%b-%d"))
  
  return(df)
}

files <- list.files(path = "/Users/dwaipayanchakrabarti/Desktop/Rutgers/Research/macrometrics/data/baseline covariates W/employment",
                    pattern = "\\.xlsx$", full.names = TRUE)

employment_panel <- map_dfr(files, import_employment)

state_unemployment_stats <- employment_panel %>%
  group_by(State) %>%
  summarize(
    Mean = mean(`unemployment rate`, na.rm = TRUE),
    Median = median(`unemployment rate`, na.rm = TRUE),
    SD = sd(`unemployment rate`, na.rm = TRUE),
    Min = min(`unemployment rate`, na.rm = TRUE),
    Max = max(`unemployment rate`, na.rm = TRUE),
    IQR = IQR(`unemployment rate`, na.rm = TRUE),
    Skewness = skewness(`unemployment rate`, na.rm = TRUE),
    Kurtosis = kurtosis(`unemployment rate`, na.rm = TRUE),
    ACF1 = if(n() > 1) cor(`unemployment rate`[2:n()], `unemployment rate`[1:(n()-1)], use = "complete.obs") else NA,
    PP_p = tryCatch(pp.test(`unemployment rate`)$p.value, error = function(e) NA),
  ) %>%
  ungroup()

state_unemployment_stats_table <- xtable(state_unemployment_stats)
print(state_unemployment_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)

state_unemployment_stats <- employment_panel %>%
  group_by(State) %>%
  summarize(
    Mean = mean(`unemployment rate`, na.rm = TRUE),
    Median = median(`unemployment rate`, na.rm = TRUE),
    SD = sd(`unemployment rate`, na.rm = TRUE),
    Min = min(`unemployment rate`, na.rm = TRUE),
    Max = max(`unemployment rate`, na.rm = TRUE),
    IQR = IQR(`unemployment rate`, na.rm = TRUE),
    Skewness = skewness(`unemployment rate`, na.rm = TRUE),
    Kurtosis = kurtosis(`unemployment rate`, na.rm = TRUE),
    ACF1 = if(n() > 1) cor(`unemployment rate`[2:n()], `unemployment rate`[1:(n()-1)], use = "complete.obs") else NA,
    PP_p = tryCatch(pp.test(`unemployment rate`)$p.value, error = function(e) NA),
  ) %>%
  ungroup()

state_unemployment_stats_table <- xtable(state_unemployment_stats)
print(state_unemployment_stats_table, type = "latex", include.rownames = FALSE, longtable = TRUE)







  
