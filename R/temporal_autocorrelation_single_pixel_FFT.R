# Spectral Analysis of Soil Water Potential Time Series
# Author: Yixuan Wang
# Date: Sys.Date() (update as needed)

# This script demonstrates how to perform spectral analysis on a soil water potential time series
# extracted from a SpatRaster. The steps include:
#   1. Data Preparation: Filtering the dataset for "Pine" species, selecting the required columns,
#      converting the 'year' column to numeric, and reshaping the data to wide format (each year as a separate column).
#   2. Conversion to SpatRaster: Creating a SpatRaster using the spatial coordinates.
#   3. Random Pixel Analysis: Randomly selecting one pixel with a valid time series.
#   4. Detrending: Removing the long-term trend via a simple linear regression.
#   5. Spectral Analysis: Computing the Fast Fourier Transform (FFT) of the detrended series,
#      calculating the power spectrum, and constructing the frequency vector.
#   6. Log-Log Transformation and Regression: Transforming frequency and power values using logarithms,
#      and fitting a linear model to estimate the spectral exponent.
#
# The spectral exponent is defined by the equation:
#   log10(P(f)) = a + b * log10(f)
# where P(f) is the power at frequency f, a is the intercept, and b is the spectral exponent.
# A more negative value of b indicates that low-frequency (long-term) fluctuations dominate.

# Load required libraries
library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)

# Set working directory and load the data
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# -----------------------------
# Data Preparation
# -----------------------------

# Filter the dataset for "Pine" species, select columns x, y, year, and soil_water_potential,
# and convert 'year' to numeric.
df_filtered <- all_results_df %>%
  filter(species == "Pine") %>% 
  select(x, y, year, soil_water_potential) %>%
  mutate(year = as.numeric(year))

# Reshape the data: pivot wider so that each year becomes a separate column.
df_wide <- df_filtered %>%
  pivot_wider(names_from = year, values_from = soil_water_potential)

# Convert the data frame to a SpatRaster.
# The first two columns (x and y) are used as spatial coordinates.
r <- rast(df_wide, type = "xyz")
print(r)

# -----------------------------
# Random Pixel Analysis
# -----------------------------

# Set seed for reproducibility.
set.seed(123)

# Convert the SpatRaster to a matrix where each row represents a pixel's time series.
r_vals <- values(r)

# Identify indices of pixels that are not all NA.
valid_cells <- which(apply(r_vals, 1, function(x) !all(is.na(x))))

# Randomly select one pixel from the valid cells.
sample_cell <- sample(valid_cells, 1)
print(paste("Selected Pixel:", sample_cell))

# -----------------------------
# Extract and Plot the Original Time Series
# -----------------------------

# Extract the time series for the selected pixel.
ts_df <- terra::extract(r, sample_cell)
ts <- as.numeric(ts_df[,-1])  # Remove the ID column

# Create a time vector (assume equally spaced time steps, e.g., one observation per year).
time <- 1:length(ts)

# Plot the original time series.
plot(time, ts, type = "b", pch = 16, col = "darkgreen",
     main = paste("Pixel", sample_cell, "- Original Time Series"),
     xlab = "Time", ylab = "Soil Water Potential")

# -----------------------------
# Detrending the Time Series
# -----------------------------

# Fit a linear regression to capture the trend in the time series.
plot(time, ts, type = "b", pch = 16, col = "darkgreen",
     main = paste("Pixel", sample_cell, "- Original Time Series"),
     xlab = "Time", ylab = "Soil Water Potential")
mod <- lm(ts ~ time)
abline(mod, col = "red", lwd = 2)

# Extract the residuals, which represent the detrended time series.
resid <- mod$residuals

# Plot the detrended time series (residuals).
plot(time, resid, type = "b", pch = 16, col = "blue",
     main = paste("Pixel", sample_cell, "- Detrended (Residuals)"),
     xlab = "Time", ylab = "Residuals")

# -----------------------------
# Detrended Data and Plot
# -----------------------------
# Compute the detrended soil water potential explicitly
detrended_psi <- ts - predict(mod)

# Plot the detrended time series (soil water potential)
plot(time, detrended_psi, type = "b", pch = 16, col = "blue",
     main = paste("Pixel", sample_cell, "- Detrended Soil Water Potential"),
     xlab = "Time", ylab = "Detrended Soil Water Potential")

# -----------------------------
# FFT and Power Spectrum on Detrended Data
# -----------------------------

# Compute the Fast Fourier Transform (FFT) of the detrended time series.
n <- length(detrended_psi)
fft_vals <- fft(detrended_psi)

# Calculate the power spectrum: the squared magnitude of the Fourier coefficients.
power <- (Mod(fft_vals))^2

# Create a frequency vector assuming a time step (Δt) of 1.
freq <- (0:(n - 1)) / n

# Plot the full power spectrum.
plot(freq, power, type = "l", col = "blue",
     main = "Power Spectrum (Detrended Data)",
     xlab = "Frequency (cycles per unit time)",
     ylab = "Power")

# -----------------------------
# Log-Log Transformation and Estimation of Spectral Exponent for Detrended Data
# -----------------------------

# Use only the positive frequencies (ignoring the zero frequency).
pos_ind <- 2:floor(n/2)
log_freq <- log10(freq[pos_ind])
log_power <- log10(power[pos_ind])

# Fit a linear regression to the log-log transformed data.
fit <- lm(log_power ~ log_freq)
slope <- coef(fit)[2]  # Spectral exponent (slope) for detrended data

# Plot the log-log plot with the fitted regression line.
plot(log_freq, log_power, pch = 16, col = "purple",
     main = paste("Pixel", sample_cell, "- Log-Log Plot (Detrended Data)"),
     xlab = "log10(Frequency)",
     ylab = "log10(Power)")
abline(fit, col = "blue", lwd = 2)
mtext(paste("Detrended Slope =", round(slope, 3)), side = 3)

# Print the computed spectral exponent for detrended data.
print(paste("Pixel", sample_cell, "spectral exponent (detrended slope):", round(slope, 3)))

# -----------------------------
# FFT and Power Spectrum on Original Data
# -----------------------------

# Compute the Fast Fourier Transform (FFT) of the original time series.
n_orig <- length(ts)
fft_vals_orig <- fft(ts)

# Calculate the power spectrum for the original data.
power_orig <- (Mod(fft_vals_orig))^2

# Create a frequency vector for the original data.
freq_orig <- (0:(n_orig - 1)) / n_orig

# Plot the full power spectrum for the original data.
plot(freq_orig, power_orig, type = "l", col = "red",
     main = "Power Spectrum (Original Data)",
     xlab = "Frequency (cycles per unit time)",
     ylab = "Power")

# -----------------------------
# Log-Log Transformation and Estimation of Spectral Exponent for Original Data
# -----------------------------

# Use only the positive frequencies (ignoring the zero frequency) for the original data.
pos_ind_orig <- 2:floor(n_orig/2)
log_freq_orig <- log10(freq_orig[pos_ind_orig])
log_power_orig <- log10(power_orig[pos_ind_orig])

# Fit a linear regression to the log-log transformed original data.
fit_orig <- lm(log_power_orig ~ log_freq_orig)
slope_orig <- coef(fit_orig)[2]  # Spectral exponent (slope) for original data

# Plot the log-log plot with the fitted regression line for original data.
plot(log_freq_orig, log_power_orig, pch = 16, col = "red",
     main = paste("Pixel", sample_cell, "- Log-Log Plot (Original Data)"),
     xlab = "log10(Frequency)",
     ylab = "log10(Power)")
abline(fit_orig, col = "blue", lwd = 2)
mtext(paste("Original Slope =", round(slope_orig, 3)), side = 3)

# Print the computed spectral exponent for original data.
print(paste("Pixel", sample_cell, "spectral exponent (original slope):", round(slope_orig, 3)))

# -----------------------------
# Conclusion
# -----------------------------
# In this script we:
#   1. Filtered and reshaped the dataset, then converted it to a SpatRaster.
#   2. Randomly selected one pixel with a valid time series.
#   3. Removed the long-term trend from the time series using linear regression.
#   4. Performed spectral analysis by computing the FFT and the power spectrum on both detrended and original data.
#   5. Transformed the frequency and power data to logarithmic scales and fitted linear regressions,
#      where the slopes of the regressions give the spectral exponents.
#
# The spectral exponent provides insight into the temporal dynamics of the soil water potential time series.
# A more negative spectral exponent indicates dominance by long-term, low-frequency fluctuations.
# Comparing the slopes from detrended and original data shows how trends affect the spectral analysis.
