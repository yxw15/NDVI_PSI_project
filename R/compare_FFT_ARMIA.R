# Example: Compare residual summary statistics for "Pine" soil_water_potential

# Load the FFT-detrended raster for Pine (each layer is one year)
r_fft <- rast(file.path(output_folder_raster_fft, "Pine_soil_water_potential_fft_detrended.tif"))
names(r_fft) <- years  # set layer names as years

# Load the auto.arima-detrended raster for Pine
r_arima <- rast(file.path(output_folder_raster_arima, "Pine_soil_water_potential_arima_detrended.tif"))
names(r_arima) <- years

# Extract the pixel time series as a matrix (rows = pixels, columns = years)
vals_fft <- values(r_fft)
vals_arima <- values(r_arima)

# Calculate per-pixel mean and standard deviation (ignoring NAs)
pixel_means_fft <- rowMeans(vals_fft, na.rm = TRUE)
pixel_sd_fft <- apply(vals_fft, 1, sd, na.rm = TRUE)

pixel_means_arima <- rowMeans(vals_arima, na.rm = TRUE)
pixel_sd_arima <- apply(vals_arima, 1, sd, na.rm = TRUE)

# Combine results into a data frame for plotting
df_res <- data.frame(
  Method = rep(c("FFT", "ARIMA"), each = length(pixel_means_fft)),
  Mean = c(pixel_means_fft, pixel_means_arima),
  SD = c(pixel_sd_fft, pixel_sd_arima)
)

# Boxplot for per-pixel Mean
boxplot(Mean ~ Method, data = df_res, 
        main = "Distribution of Pixel Means\n(Pine soil_water_potential)",
        ylab = "Mean Value", col = c("#F0E442", "#0072B2"))

# Boxplot for per-pixel Standard Deviation
boxplot(SD ~ Method, data = df_res, 
        main = "Distribution of Pixel Standard Deviations\n(Pine soil_water_potential)",
        ylab = "Standard Deviation", col = c("#F0E442", "#0072B2"))

# Alternatively, using ggplot2 for a combined view:
library(ggplot2)

# Plot histogram of Means
ggplot(df_res, aes(x = Mean, fill = Method)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  labs(title = "Histogram of Pixel Means (Pine soil_water_potential)",
       x = "Mean", y = "Count") +
  scale_fill_manual(values = c("FFT" = "#F0E442", "ARIMA" = "#0072B2"))

# Plot histogram of Standard Deviations
ggplot(df_res, aes(x = SD, fill = Method)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  labs(title = "Histogram of Pixel SDs (Pine soil_water_potential)",
       x = "Standard Deviation", y = "Count") +
  scale_fill_manual(values = c("FFT" = "#F0E442", "ARIMA" = "#0072B2"))
