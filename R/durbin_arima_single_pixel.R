# -----------------------------
# Load required libraries
# -----------------------------
library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(lmtest)     # For Durbin-Watson test
library(orcutt)     # For Cochrane-Orcutt procedure
library(forecast)   # For auto.arima

# -----------------------------
# Set working directory and load the data
# -----------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# -----------------------------
# Data Preparation
# -----------------------------
# Filter dataset for "Pine", select columns and convert year to numeric
df_filtered <- all_results_df %>%
  filter(species == "Pine") %>% 
  select(x, y, year, soil_water_potential) %>%
  mutate(year = as.numeric(year))

# Reshape the data: pivot wider so that each year becomes a separate column
df_wide <- df_filtered %>%
  pivot_wider(names_from = year, values_from = soil_water_potential)

# Convert the data frame to a SpatRaster (using x and y as coordinates)
r <- rast(df_wide, type = "xyz")
print(r)

# -----------------------------
# Random Pixel Analysis
# -----------------------------
# Set seed for reproducibility and select a random pixel with valid data
set.seed(123)
r_vals <- values(r)
valid_cells <- which(apply(r_vals, 1, function(x) !all(is.na(x))))
sample_cell <- sample(valid_cells, 1)
print(paste("Selected Pixel:", sample_cell))

# Extract the time series for the selected pixel
ts_df <- terra::extract(r, sample_cell)
ts <- as.numeric(ts_df[,-1])  # Remove the ID column
time <- 1:length(ts)

# -----------------------------
# Plot Original Time Series and ACF
# -----------------------------
# Save original time series plot
png(filename = "results_temporal/compare/original_ts.png", width = 800, height = 600)
plot(time, ts, type = "b", pch = 16, col = "darkgreen",
     main = paste("Pixel", sample_cell, "- Original Time Series"),
     xlab = "Time", ylab = "Soil Water Potential")
dev.off()

# Save ACF plot of original time series
png(filename = "results_temporal/compare/original_acf.png", width = 800, height = 600)
acf(ts, main = "ACF - Original Time Series")
dev.off()

# -----------------------------
# Test for Autocorrelation using Linear Model
# -----------------------------
model_lm <- lm(ts ~ time)
dw_result <- dwtest(model_lm)
print(dw_result)

# -----------------------------
# Method 1: Cochrane-Orcutt Procedure
# -----------------------------
model_co <- cochrane.orcutt(model_lm)
summary(model_co)
detrended_co <- resid(model_co)

png(filename = "results_temporal/compare/cochrane_ts.png", width = 800, height = 600)
plot(time, detrended_co, type = "b", pch = 16, col = "blue",
     main = "Detrended Time Series (Cochrane-Orcutt)",
     xlab = "Time", ylab = "Residuals")
dev.off()

png(filename = "results_temporal/compare/cochrane_acf.png", width = 800, height = 600)
acf(detrended_co, main = "ACF - Cochrane-Orcutt Residuals")
dev.off()

# -----------------------------
# Method 2: ARMA(1,1) Model
# -----------------------------
model_arma <- arima(ts, order = c(1, 0, 1))
detrended_arma <- residuals(model_arma)

png(filename = "results_temporal/compare/arma_ts.png", width = 800, height = 600)
plot(time, detrended_arma, type = "b", pch = 16, col = "blue",
     main = "Detrended Time Series (ARMA(1,1))",
     xlab = "Time", ylab = "Residuals")
dev.off()

png(filename = "results_temporal/compare/arma_acf.png", width = 800, height = 600)
acf(detrended_arma, main = "ACF - ARMA(1,1) Residuals")
dev.off()

# -----------------------------
# Method 3: Auto ARIMA Model
# -----------------------------
model_auto <- auto.arima(ts, seasonal = FALSE)
summary(model_auto)
detrended_auto <- residuals(model_auto)

png(filename = "results_temporal/compare/autoarima_ts.png", width = 800, height = 600)
plot(time, detrended_auto, type = "b", pch = 16, col = "blue",
     main = "Detrended Time Series (auto.arima)",
     xlab = "Time", ylab = "Residuals")
dev.off()

png(filename = "results_temporal/compare/autoarima_acf.png", width = 800, height = 600)
acf(detrended_auto, main = "ACF - auto.arima Residuals")
dev.off()

# -----------------------------
# Method 4: FFT-based Detrending
# -----------------------------
# Compute FFT of the original time series
fft_vals <- fft(ts)
n <- length(ts)
freq <- (0:(n-1)) / n

# Choose a cutoff frequency for the low-frequency components (trend)
cutoff <- 0.2  # Adjust as needed based on your data characteristics
fft_vals_filtered <- fft_vals
# Zero out the components with frequency less than the cutoff (including the DC component)
fft_vals_filtered[freq < cutoff] <- 0

# Compute the inverse FFT to obtain the FFT-based detrended series
detrended_fft <- Re(fft(fft_vals_filtered, inverse = TRUE)) / n

png(filename = "results_temporal/compare/fft_detrended_ts.png", width = 800, height = 600)
plot(time, detrended_fft, type = "b", pch = 16, col = "blue",
     main = "Detrended Time Series (FFT-based)",
     xlab = "Time", ylab = "Detrended Value")
dev.off()

png(filename = "results_temporal/compare/fft_detrended_acf.png", width = 800, height = 600)
acf(detrended_fft, main = "ACF - FFT-based Detrended Series")
dev.off()

# -----------------------------
# Comparison of Methods
# -----------------------------
# Combine detrended series from all methods into one data frame for comparison
# Note: For consistency, we refer to the "detrended" results as follows:
#   - Cochrane-Orcutt: residuals from linear regression with autocorrelation adjustment
#   - ARMA(1,1): residuals from ARMA model
#   - auto.arima: residuals from automatically selected ARIMA model
#   - FFT-based: detrended series via FFT filtering
residuals_df <- data.frame(
  Time = time,
  CochraneOrcutt = detrended_co,
  ARMA11 = detrended_arma,
  AutoARIMA = detrended_auto,
  FFT = detrended_fft
)

# Calculate summary statistics for detrended methods
residuals_stats <- data.frame(
  Method = c("Cochrane-Orcutt", "ARMA(1,1)", "auto.arima", "FFT-based"),
  Mean = c(mean(detrended_co, na.rm = TRUE),
           mean(detrended_arma, na.rm = TRUE),
           mean(detrended_auto, na.rm = TRUE),
           mean(detrended_fft, na.rm = TRUE)),
  SD = c(sd(detrended_co, na.rm = TRUE),
         sd(detrended_arma, na.rm = TRUE),
         sd(detrended_auto, na.rm = TRUE),
         sd(detrended_fft, na.rm = TRUE))
)
print(residuals_stats)

# Compute summary statistics for the original time series
original_mean <- mean(ts, na.rm = TRUE)
original_sd <- sd(ts, na.rm = TRUE)

# Create text annotation including original data and detrended methods
stats_text <- paste0(
  "Original: Mean = ", round(original_mean, 2), ", SD = ", round(original_sd, 2), "\n",
  "Cochrane-Orcutt: Mean = ", round(residuals_stats$Mean[1], 2),
  ", SD = ", round(residuals_stats$SD[1], 2), "\n",
  "ARMA(1,1): Mean = ", round(residuals_stats$Mean[2], 2),
  ", SD = ", round(residuals_stats$SD[2], 2), "\n",
  "auto.arima: Mean = ", round(residuals_stats$Mean[3], 2),
  ", SD = ", round(residuals_stats$SD[3], 2), "\n",
  "FFT-based: Mean = ", round(residuals_stats$Mean[4], 2),
  ", SD = ", round(residuals_stats$SD[4], 2)
)

# Create a composite plot combining time series and ACF for all methods
# 5 rows (Original, Cochrane-Orcutt, ARMA(1,1), auto.arima, FFT-based) and 2 columns (Time Series, ACF)
png(filename = "results_temporal/compare/composite_with_fft_regression.png", width = 800, height = 1000)
par(mfrow = c(5,2), mar = c(4,4,3,1))

# Row 1: Original Data
plot(time, ts, type = "b", pch = 16, col = "purple",
     main = paste("Pixel", sample_cell, "- Original Time Series"),
     xlab = "Time", ylab = "Soil Water Potential")
abline(lm(ts ~ time), col = "red", lwd = 2)  # Add regression line
acf(ts, main = "ACF - Original")

# Row 2: Cochrane-Orcutt
plot(time, detrended_co, type = "b", pch = 16, col = "blue",
     main = "Detrended (Cochrane-Orcutt)",
     xlab = "Time", ylab = "Residuals")
abline(lm(detrended_co ~ time), col = "red", lwd = 2)  # Add regression line
acf(detrended_co, main = "ACF - Cochrane-Orcutt")

# Row 3: ARMA(1,1)
plot(time, detrended_arma, type = "b", pch = 16, col = "blue",
     main = "Detrended (ARMA(1,1))",
     xlab = "Time", ylab = "Residuals")
abline(lm(detrended_arma ~ time), col = "red", lwd = 2)  # Add regression line
acf(detrended_arma, main = "ACF - ARMA(1,1)")

# Row 4: auto.arima
plot(time, detrended_auto, type = "b", pch = 16, col = "blue",
     main = "Detrended (auto.arima)",
     xlab = "Time", ylab = "Residuals")
abline(lm(detrended_auto ~ time), col = "red", lwd = 2)  # Add regression line
acf(detrended_auto, main = "ACF - auto.arima")

# Row 5: FFT-based
# Left panel: FFT-based detrended time series with annotation added on lower left corner
plot(time, detrended_fft, type = "b", pch = 16, col = "blue",
     main = "Detrended (FFT-based)",
     xlab = "Time", ylab = "Detrended Value")
abline(lm(detrended_fft ~ time), col = "red", lwd = 2)  # Add regression line
usr <- par("usr")  # Get current plot boundaries: [xmin, xmax, ymin, ymax]
text_x <- usr[1] + 0.05 * (usr[2] - usr[1])  # 5% from the left edge
text_y <- usr[3] + 0.05 * (usr[4] - usr[3])  # 5% from the bottom edge
text(x = text_x, y = text_y, labels = stats_text, adj = c(0,0), cex = 1.5)

# Right panel: FFT-based detrended ACF
acf(detrended_fft, main = "ACF - FFT-based")

dev.off()
