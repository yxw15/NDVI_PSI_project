# Set working directory and load the data
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# Load required packages
library(dplyr)
library(tidyr)
library(terra)

# 1. Data Preparation --------------------------------------------------

# Filter for Pine species, select desired columns, and convert year to numeric
df_filtered <- all_results_df %>%
  filter(species == "Pine") %>% 
  select(x, y, year, soil_water_potential) %>%
  mutate(year = as.numeric(year))

# Pivot the data so that each unique year becomes its own column (wide format)
df_wide <- df_filtered %>%
  pivot_wider(names_from = year, values_from = soil_water_potential)

# Convert the wide dataframe into a SpatRaster.
# The first two columns (x, y) are interpreted as coordinates.
r <- rast(df_wide, type = "xyz")
print(r)

# 2. Define Analysis Functions -----------------------------------------

# Function to compute the spectral exponent for each pixel's time series.
# It detrends the series, computes the FFT, and fits a regression in log-log space.
spectral_exponent_fun <- function(ts) {
  # Return NA if the entire time series is missing
  if (all(is.na(ts))) return(NA)
  
  # Create a time vector assuming equally spaced time steps (e.g., 1, 2, ..., T)
  time <- 1:length(ts)
  
  # Detrend: Fit a linear model (ts ~ time) and extract residuals.
  mod <- lm(ts ~ time)
  resid <- mod$residuals
  
  n <- length(resid)
  # Compute the FFT of the detrended time series (residuals)
  fft_vals <- fft(resid)
  power <- (Mod(fft_vals))^2  # Power spectrum
  
  # Create a frequency vector (assuming dt = 1)
  freq <- (0:(n-1)) / n
  
  # Use only positive frequencies (skip zero frequency)
  pos_ind <- 2:floor(n/2)
  log_freq <- log10(freq[pos_ind])
  log_power <- log10(power[pos_ind])
  
  # Fit a linear regression to the log-log plot to get the spectral exponent (slope)
  fit <- lm(log_power ~ log_freq)
  slope <- coef(fit)[2]
  
  return(slope)
}

# Function to detrend each pixel's time series by removing the linear trend.
# It returns the residuals, which represent fluctuations around the long-term trend.
detrend_fun <- function(ts) {
  # Return a vector of NAs if the time series is entirely missing
  if (all(is.na(ts))) {
    return(rep(NA, length(ts)))
  }
  
  # Create a time vector assuming equally spaced time steps
  time <- 1:length(ts)
  
  # Fit a linear model to remove the trend
  model <- lm(ts ~ time)
  
  # The residuals represent the detrended data
  detrended <- model$residuals
  
  return(detrended)
}

# 3. Apply the Analysis Functions ---------------------------------------

# Apply the spectral exponent function to each pixel's time series.
# This produces a single-layer SpatRaster where each pixel's value is the spectral exponent.
spectral_exponent_rast <- app(r, spectral_exponent_fun)
plot(spectral_exponent_rast, main = "Spectral Exponent per Pixel")

# Apply the detrending function to each pixel's time series.
# This produces a new SpatRaster stack with the same number of layers as 'r',
# where each layer corresponds to the detrended (residual) value for that time step.
detrended_rast <- app(r, detrend_fun)
plot(detrended_rast, main = "Detrended Time Series per Pixel")
