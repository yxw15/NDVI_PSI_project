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

# Set seed for reproducibility
set.seed(123)

# -----------------------------------------------------------------
# STEP 1: Identify One Random Pixel with a Valid Time Series
# -----------------------------------------------------------------

# Convert the SpatRaster to a matrix where rows are pixels
r_vals <- values(r)

# Identify indices of pixels that are not all NA
valid_cells <- which(apply(r_vals, 1, function(x) !all(is.na(x))))

# Randomly select 1 pixel from the valid cells
sample_cell <- sample(valid_cells, 1)
print(paste("Selected Pixel:", sample_cell))

# -----------------------------------------------------------------
# STEP 2: Extract and Plot the Original Time Series
# -----------------------------------------------------------------

# Extract the time series for the selected pixel.
# terra::extract returns a data.frame with an ID column, so we remove it.
ts_df <- terra::extract(r, sample_cell)
ts <- as.numeric(ts_df[,-1])

# Create a time vector (assuming equally spaced time steps, e.g., 1, 2, 3, ...)
time <- 1:length(ts)

# Plot the original time series
plot(time, ts, type = "b", pch = 16, col = "darkgreen",
     main = paste("Pixel", sample_cell, "- Original Time Series"),
     xlab = "Time", ylab = "Soil Water Potential")

# -----------------------------------------------------------------
# STEP 3: Detrend the Time Series (Remove Long-Term Trend)
# -----------------------------------------------------------------

# Fit a linear regression (ts ~ time) to capture the trend
mod <- lm(ts ~ time)

# Overlay the fitted trend line on the original plot
abline(mod, col = "red", lwd = 2)

# Extract the residuals, which represent the detrended data
resid <- mod$residuals

# Plot the detrended time series (residuals)
plot(time, resid, type = "b", pch = 16, col = "blue",
     main = paste("Pixel", sample_cell, "- Detrended (Residuals)"),
     xlab = "Time", ylab = "Residuals")

# -----------------------------------------------------------------
# STEP 4: Compute the FFT and Power Spectrum of the Detrended Data
# -----------------------------------------------------------------

n <- length(resid)
fft_vals <- fft(resid)                # Compute the Fast Fourier Transform
power <- (Mod(fft_vals))^2            # Calculate the power spectrum

# Create a frequency vector (assuming a time step dt = 1)
freq <- (0:(n - 1)) / n

# Plot 1: Power Spectrum vs Frequency (full range)
plot(freq, power, type = "l", col = "blue",
     main = "Power Spectrum",
     xlab = "Frequency (cycles per unit time)",
     ylab = "Power")

# -----------------------------------------------------------------
# STEP 5: Focus on Positive Frequencies and Log-Transform Data
# -----------------------------------------------------------------

# Use only positive frequencies (skip the zero frequency)
pos_ind <- 2:floor(n/2)
log_freq <- log10(freq[pos_ind])
log_power <- log10(power[pos_ind])

# -----------------------------------------------------------------
# STEP 6: Fit a Linear Regression to the Log-Log Data to Obtain the Spectral Exponent
# -----------------------------------------------------------------

fit <- lm(log_power ~ log_freq)
slope <- coef(fit)[2]  # The spectral exponent

# Plot the log-log plot with the fitted line
plot(log_freq, log_power, pch = 16, col = "purple",
     main = paste("Pixel", sample_cell, "- Log-Log Plot"),
     xlab = "log10(Frequency)", ylab = "log10(Power)")
abline(fit, col = "blue", lwd = 2)
mtext(paste("Slope =", round(slope, 3)), side = 3)

# Print the computed spectral exponent
print(paste("Pixel", sample_cell, "spectral exponent (slope):", round(slope, 3)))
