library(terra)
library(tseries)
library(forecast)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_NDVI_PSI_TDiff_species.RData")

# Define the path
Beech_PSI_mask_folder <- "results/Beech/PSI_mask"

# List all .tif files that match PSI_<year>_mask.tif pattern
tif_files <- list.files(Beech_PSI_mask_folder, pattern = "^PSI_\\d{4}_mask\\.tif$", full.names = TRUE)

# Extract the years from filenames
tif_years <- as.numeric(gsub("^.*PSI_(\\d{4})_mask\\.tif$", "\\1", tif_files))

# Sort files and years by year
sorted_idx <- order(tif_years)
tif_files <- tif_files[sorted_idx]
tif_years <- tif_years[sorted_idx]

Beech_PSI_stack <- rast(tif_files)

# Rename layers by year
names(Beech_PSI_stack) <- as.character(tif_years)

# Optional: check the result
print(Beech_PSI_stack)

# Calculate mean value for each layer (year)
beech_yearly_means <- global(Beech_PSI_stack, fun = "mean", na.rm = TRUE)

# Add year names as a column
beech_yearly_means$year <- as.numeric(names(Beech_PSI_stack))

# Print results
print(beech_yearly_means)

# Create a time series object
psi_ts <- ts(beech_yearly_means$mean, start = min(beech_yearly_means$year), frequency = 1)

# Run Augmented Dickey-Fuller test
adf_result <- adf.test(psi_ts)

# Print result
print(adf_result)

acf(psi_ts, main = "ACF of PSI Mean Time Series")

Box.test(psi_ts, lag = 6, type = "Ljung-Box")

pacf(psi_ts, main = "Partial ACF of PSI Mean Time Series")


auto_model <- auto.arima(diff_psi_ts)
summary(auto_model)

# Extract residuals
res <- residuals(auto_model)

# Check autocorrelation
acf(res, main = "ACF of Residuals")  # Should be all inside dashed lines

# Formal test: Ljung-Box
Box.test(res, lag = 10, type = "Ljung-Box")  # p > 0.05 = no autocorrelation

adf.test(psi_ts)
acf(psi_ts, main = "ACF of PSI Mean Time Series")


fft_detrend_fun <- function(ts) {
  if (all(is.na(ts))) return(rep(NA, length(ts)))
  n <- length(ts)
  fft_vals <- fft(ts)
  freq <- (0:(n - 1)) / n
  fft_vals_filtered <- fft_vals
  cutoff <- 0.15  # Adjust cutoff frequency as needed
  fft_vals_filtered[freq < cutoff] <- 0
  detrended_ts <- Re(fft(fft_vals_filtered, inverse = TRUE)) / n
  detrended_ts
}

psi_detrended <- fft_detrend_fun(psi_ts)
plot(psi_ts, type = "l", main = "Original vs. Detrended", col = "black")
lines(psi_detrended, col = "red")
legend("topright", legend = c("Original", "FFT Detrended"), col = c("black", "red"), lty = 1)

acf(psi_detrended, main = "ACF of detrended")  # Should be all inside dashed lines


