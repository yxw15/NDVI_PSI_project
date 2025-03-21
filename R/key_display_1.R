# ------------------------------------------------------------------------------
# NDVI/PSI Temporal Analysis Script (Revised & Organized)
# Author: Yixuan Wang (modified by user)
# Date: Sys.Date()
#
# Description:
# This script performs a comprehensive temporal analysis of NDVI/PSI data. 
# It processes both the original and FFT-detrended datasets to explore:
#
#   1. Correlation Analysis:
#      - Computes correlation matrices and visualizes them using corrplot.
#
#   2. Time Series & Scatter Plot Analysis:
#      - Generates multi-panel time series plots and correlation scatter plots.
#
#   3. Temporal Autocorrelation (ACF) Analysis:
#      - Produces ACF plots for averaged data as well as species-level data.
#
#   4. Power Spectrum Analysis:
#      - Applies Fast Fourier Transform (FFT) to compute power spectra 
#        and estimates spectral slopes.
#
# The script is divided into two main sections:
#   - Section 1: Analysis on Original Data.
#   - Section 2: Analysis on Detrended Data.
#
# External functions are sourced from:
#   - R/functions_NDVI.R
#   - R/functions_temporal_autocorrelation.R
#
# All key displays and figures are saved in the designated output directory.
# ------------------------------------------------------------------------------

# ------------------------------
# Setup: Working Directory, Data, Functions, and Output Folder
# ------------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Source required functions
source("R/functions_NDVI.R")
source("R/functions_temporal_autocorrelation.R")

# Define output directory for key displays and create it if necessary
output_dir <- "results/key_displays"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load required libraries
library(ggplot2)
library(patchwork)
if (!require("corrplot")) {
  install.packages("corrplot")
  library(corrplot)
}

# ------------------------------------------------------------------------------
# Section 1: Analysis on Original Data
# ------------------------------------------------------------------------------

# Load the combined original data
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

## 1A. All Data: Correlation Matrix
# Compute averaged data and create correlation matrix
df <- df_average_year(all_results_df)
cor_matrix <- cor(df[, c("avg_quantile", "avg_psi", "avg_tdiff")], use = "complete.obs")
cat("Correlation matrix for averaged (original) data:\n")
print(cor_matrix)

# Save the correlation matrix plot
png(filename = file.path(output_dir, "correlation_matrix_original.png"), 
    width = 800, height = 600)
corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         order = "hclust",
         addCoef.col = "white",    # White coefficient labels
         tl.col = "black", 
         number.cex = 2,           # Coefficient text size set to 2
         tl.srt = 45, 
         diag = FALSE,
         title = "Correlation Matrix (Original)", 
         mar = c(0, 0, 1, 0))
dev.off()

## 1B. All Data: Time Series & Correlation Scatter Plots
# Generate and save multi-panel time series and scatter plots
plot_time_series_NDVI_PSI_TDiff_avg(
  df_all = all_results_df, 
  plot_ts_path = file.path(output_dir, "time_series_three_panel_original.png")
)

plot_correlation_NDVI_PSI_TDiff_corr(
  df_all = all_results_df, 
  plot_corr_path = file.path(output_dir, "correlation_two_panel_original.png")
)

## 1C. All Data: ACF Plots
# Create ACF plots for averaged data across three variables
png(filename = file.path(output_dir, "acf_plots_original.png"), width = 1200, height = 400)
par(mfrow = c(1, 3))
acf(df$avg_quantile, main = "ACF: avg_quantile (Original)")
acf(df$avg_psi, main = "ACF: avg_psi (Original)")
acf(df$avg_tdiff, main = "ACF: avg_tdiff (Original)")
par(mfrow = c(1, 1))
dev.off()

## 1D. All Data: Power Spectrum Analysis
# Define a function to compute and plot the power spectrum via FFT
compute_power_spectrum_gg <- function(ts_data, variable_name) {
  n <- length(ts_data)
  fft_vals <- fft(ts_data)
  power <- (Mod(fft_vals))^2
  freq <- (0:(n - 1)) / n
  pos_ind <- 2:floor(n / 2)
  log_freq <- log10(freq[pos_ind])
  log_power <- log10(power[pos_ind])
  df_fft <- data.frame(log_freq = log_freq, log_power = log_power)
  fit <- lm(log_power ~ log_freq, data = df_fft)
  slope <- coef(fit)[2]
  p <- ggplot(df_fft, aes(x = log_freq, y = log_power)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(x = "log10(Frequency)", y = "log10(Power)",
         title = paste("Power Spectrum for", variable_name)) +
    annotate("text", x = min(log_freq), y = max(log_power),
             label = paste("Slope:", round(slope, 2)),
             hjust = 0, vjust = 1, size = 5) +
    theme_minimal()
  cat("Spectral slope for", variable_name, slope, "\n")
  return(list(slope = slope, plot = p))
}

# Compute power spectra for each variable and combine the plots
ps_quantile <- compute_power_spectrum_gg(df$avg_quantile, "avg_quantile")
ps_psi      <- compute_power_spectrum_gg(df$avg_psi, "avg_psi")
ps_tdiff    <- compute_power_spectrum_gg(df$avg_tdiff, "avg_tdiff")
combined_power_original <- ps_quantile$plot / ps_psi$plot / ps_tdiff$plot
ggsave(filename = file.path(output_dir, "power_spectrum_original_all.png"), 
       plot = combined_power_original, width = 8, height = 12, dpi = 300)

## 1E. Species-Level Analysis (Original Data)
# Generate species-level time series and correlation scatter plots
plot_time_series_NDVI_PSI_TDiff_species_avg(df_all = all_results_df, out_dir = output_dir)
plot_correlation_NDVI_PSI_TDiff_species_corr(df_all = all_results_df, out_dir = output_dir)

# Create species-level ACF plots
df_species <- df_average_year_species(all_results_df)
unique_species <- unique(df_species$species)
for (sp in unique_species) {
  df_sp <- subset(df_species, species == sp)
  png(filename = file.path(output_dir, paste0("acf_plots_original_", sp, ".png")), 
      width = 1200, height = 400)
  par(mfrow = c(1, 3))
  acf(df_sp$avg_quantile, main = paste("ACF:", sp, "avg_quantile (Original)"))
  acf(df_sp$avg_psi, main = paste("ACF:", sp, "avg_psi (Original)"))
  acf(df_sp$avg_tdiff, main = paste("ACF:", sp, "avg_tdiff (Original)"))
  par(mfrow = c(1, 1))
  dev.off()
}

# Species-level power spectrum analysis
for (sp in unique_species) {
  df_sp <- subset(df_species, species == sp)
  ps_quantile_sp <- compute_power_spectrum_gg(df_sp$avg_quantile, paste("avg_quantile -", sp, "(Original)"))
  ps_psi_sp      <- compute_power_spectrum_gg(df_sp$avg_psi, paste("avg_psi -", sp, "(Original)"))
  ps_tdiff_sp    <- compute_power_spectrum_gg(df_sp$avg_tdiff, paste("avg_tdiff -", sp, "(Original)"))
  combined_power_sp <- ps_quantile_sp$plot / ps_psi_sp$plot / ps_tdiff_sp$plot
  ggsave(filename = file.path(output_dir, paste0("power_spectrum_original_", sp, ".png")),
         plot = combined_power_sp, width = 8, height = 12, dpi = 300)
}

## 1F. Combined Species-Specific Correlation Plots (Original Data)
# Define species order for the combined plot (lowercase for matching)
species_order <- c("oak", "beech", "spruce", "pine")
png(filename = file.path(output_dir, "species_correlation_matrix.png"), width = 1600, height = 1200)
par(mfrow = c(2, 2), bg = "white")
for (sp in species_order) {
  sp_data <- subset(df_species, tolower(species) == sp)
  if(nrow(sp_data) == 0) next
  cor_matrix_sp <- cor(sp_data[, c("avg_quantile", "avg_psi", "avg_tdiff")], use = "complete.obs")
  corrplot(cor_matrix_sp, 
           method = "color", 
           type = "upper", 
           order = "hclust",
           addCoef.col = "white",   # White coefficient labels
           number.cex = 2,          # Coefficient text size set to 2
           tl.col = "black",        # Black text labels
           tl.srt = 0, 
           diag = FALSE,
           main = paste("Correlation Matrix (", sp, ")", sep = ""),
           mar = c(0,0,1,0))
}
par(mfrow = c(1, 1))
dev.off()

# ------------------------------------------------------------------------------
# Section 2: Analysis on Detrended Data
# ------------------------------------------------------------------------------

# Load the FFT-detrended combined data
load("results/Data/All_Quantiles_PSI_TDiff_species_year_fft_detrended.RData")

# Source NDVI functions again if necessary
source("R/functions_NDVI.R")

## 2A. Create Detrended Data for All Variables
df_det <- df_average_year(final_df_clean)

## 2B. All Data: Correlation Matrix (Detrended)
cor_matrix_det <- cor(df_det[, c("avg_quantile", "avg_psi", "avg_tdiff")], use = "complete.obs")
cat("Correlation matrix for detrended averaged data:\n")
print(cor_matrix_det)
png(filename = file.path(output_dir, "correlation_matrix_detrended.png"), 
    width = 800, height = 600)
corrplot(cor_matrix_det, 
         method = "color", 
         type = "upper", 
         order = "hclust",
         addCoef.col = "white",   # White coefficient labels
         number.cex = 2,          # Coefficient text size set to 2
         tl.col = "black", 
         tl.srt = 45, 
         diag = FALSE,
         title = "Correlation Matrix (Detrended)", 
         mar = c(0,0,1,0))
dev.off()

## 2C. All Data: Time Series & Correlation Scatter Plots (Detrended)
# Determine NDVI column and label for detrended data
ndvi_column <- if ("avg_quantile" %in% names(df_det)) "avg_quantile" else "avg_proportion"
ndvi_label <- if (ndvi_column == "avg_quantile") "Average NDVI (Quantiles)" else "Average NDVI (Proportions)"

# Time Series Plots
p1_det <- ggplot(df_det, aes(x = as.numeric(year), y = .data[[ndvi_column]])) +
  geom_point(color = "orange", size = 3, alpha = 0.7) +
  geom_line(color = "orange") +
  labs(title = "NDVI over Years (Detrended)", x = "Year", y = ndvi_label) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

p2_det <- ggplot(df_det, aes(x = as.numeric(year), y = avg_psi)) +
  geom_point(color = "purple", size = 3, alpha = 0.7) +
  geom_line(color = "purple") +
  labs(title = "Soil Water Potential over Years (Detrended)", x = "Year", y = "Average PSI") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

p3_det <- ggplot(df_det, aes(x = as.numeric(year), y = avg_tdiff)) +
  geom_point(color = "green4", size = 3, alpha = 0.7) +
  geom_line(color = "green4") +
  labs(title = "Transpiration Deficit over Years (Detrended)", x = "Year", y = "Average TDiff") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

combined_ts_det <- p1_det / p2_det / p3_det
ggsave(filename = file.path(output_dir, "time_series_three_panel_detrended.png"),
       plot = combined_ts_det, width = 10, height = 15, dpi = 300)

# Correlation Scatter Plots (Detrended)
correlation_psi_det <- cor(df_det$avg_psi, df_det[[ndvi_column]], use = "complete.obs")
plot_psi_det <- ggplot(df_det, aes(x = avg_psi, y = .data[[ndvi_column]])) +
  geom_point(color = "purple", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  annotate("text", x = min(df_det$avg_psi, na.rm = TRUE),
           y = max(df_det[[ndvi_column]], na.rm = TRUE),
           label = paste("Correlation:", round(correlation_psi_det, 2)),
           hjust = 0, vjust = 1.5, size = 5, color = "blue") +
  labs(title = "NDVI vs Soil Water Potential (Detrended)", x = "Average PSI", y = ndvi_label) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

correlation_tdiff_det <- cor(df_det$avg_tdiff, df_det[[ndvi_column]], use = "complete.obs")
plot_tdiff_det <- ggplot(df_det, aes(x = avg_tdiff, y = .data[[ndvi_column]])) +
  geom_point(color = "green4", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  annotate("text", x = min(df_det$avg_tdiff, na.rm = TRUE),
           y = max(df_det[[ndvi_column]], na.rm = TRUE),
           label = paste("Correlation:", round(correlation_tdiff_det, 2)),
           hjust = 0, vjust = 1.5, size = 5, color = "blue") +
  labs(title = "NDVI vs Transpiration Deficit (Detrended)", x = "Average TDiff", y = ndvi_label) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

combined_corr_det <- plot_psi_det + plot_tdiff_det + plot_layout(ncol = 2, guides = "collect")
ggsave(filename = file.path(output_dir, "correlation_two_panel_detrended.png"),
       plot = combined_corr_det, width = 15, height = 7, dpi = 300)

# ACF Plots for detrended all data
png(filename = file.path(output_dir, "acf_plots_detrended.png"), width = 1200, height = 400)
par(mfrow = c(1, 3))
acf(df_det$avg_quantile, main = "ACF: avg_quantile (Detrended)")
acf(df_det$avg_psi, main = "ACF: avg_psi (Detrended)")
acf(df_det$avg_tdiff, main = "ACF: avg_tdiff (Detrended)")
par(mfrow = c(1, 1))
dev.off()

# Power Spectrum for detrended all data
ps_quantile_det <- compute_power_spectrum_gg(df_det$avg_quantile, "avg_quantile (Detrended)")
ps_psi_det      <- compute_power_spectrum_gg(df_det$avg_psi, "avg_psi (Detrended)")
ps_tdiff_det    <- compute_power_spectrum_gg(df_det$avg_tdiff, "avg_tdiff (Detrended)")
combined_power_det <- ps_quantile_det$plot / ps_psi_det$plot / ps_tdiff_det$plot
ggsave(filename = file.path(output_dir, "power_spectrum_detrended_all.png"),
       plot = combined_power_det, width = 8, height = 12, dpi = 300)

## 2D. Species-Level Analysis (Detrended)
# Detrend species-level data using the provided detrending function
df_species_det <- df_species
df_species_det$avg_quantile <- detrend_fun(df_species_det$avg_quantile)
df_species_det$avg_psi      <- detrend_fun(df_species_det$avg_psi)
df_species_det$avg_tdiff    <- detrend_fun(df_species_det$avg_tdiff)

unique_species <- unique(df_species_det$species)
for (sp in unique_species) {
  df_sp_det <- subset(df_species_det, species == sp)
  
  # --- Time Series for Species (Detrended)
  p1_sp_det <- ggplot(df_sp_det, aes(x = as.numeric(year), y = avg_quantile)) +
    geom_point(color = "orange", size = 3, alpha = 0.7) +
    geom_line(color = "orange") +
    labs(title = paste("NDVI over Years (Detrended) -", sp),
         x = "Year", y = "Average NDVI (Quantiles)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  p2_sp_det <- ggplot(df_sp_det, aes(x = as.numeric(year), y = avg_psi)) +
    geom_point(color = "purple", size = 3, alpha = 0.7) +
    geom_line(color = "purple") +
    labs(title = paste("Soil Water Potential (Detrended) -", sp),
         x = "Year", y = "Average PSI") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  p3_sp_det <- ggplot(df_sp_det, aes(x = as.numeric(year), y = avg_tdiff)) +
    geom_point(color = "green4", size = 3, alpha = 0.7) +
    geom_line(color = "green4") +
    labs(title = paste("Transpiration Deficit (Detrended) -", sp),
         x = "Year", y = "Average TDiff") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  combined_ts_sp_det <- p1_sp_det / p2_sp_det / p3_sp_det
  ggsave(filename = file.path(output_dir, paste0("time_series_three_panel_detrended_", sp, ".png")),
         plot = combined_ts_sp_det, width = 10, height = 15, dpi = 300)
  
  # --- Correlation Scatter Plots for Species (Detrended)
  correlation_psi_sp_det <- cor(df_sp_det$avg_psi, df_sp_det[[ndvi_column]], use = "complete.obs")
  plot_psi_sp_det <- ggplot(df_sp_det, aes(x = avg_psi, y = .data[[ndvi_column]])) +
    geom_point(color = "purple", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", x = min(df_sp_det$avg_psi, na.rm = TRUE), 
             y = max(df_sp_det[[ndvi_column]], na.rm = TRUE), 
             label = paste("Correlation:", round(correlation_psi_sp_det, 2)),
             hjust = 0, vjust = 1.5, size = 5, color = "blue") +
    labs(title = paste("NDVI vs PSI (Detrended) -", sp),
         x = "Average PSI", y = ndvi_label) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  correlation_tdiff_sp_det <- cor(df_sp_det$avg_tdiff, df_sp_det[[ndvi_column]], use = "complete.obs")
  plot_tdiff_sp_det <- ggplot(df_sp_det, aes(x = avg_tdiff, y = .data[[ndvi_column]])) +
    geom_point(color = "green4", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", x = min(df_sp_det$avg_tdiff, na.rm = TRUE), 
             y = max(df_sp_det[[ndvi_column]], na.rm = TRUE), 
             label = paste("Correlation:", round(correlation_tdiff_sp_det, 2)),
             hjust = 0, vjust = 1.5, size = 5, color = "blue") +
    labs(title = paste("NDVI vs TDiff (Detrended) -", sp),
         x = "Average TDiff", y = ndvi_label) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  combined_corr_sp_det <- plot_psi_sp_det + plot_tdiff_sp_det + plot_layout(ncol = 2, guides = "collect")
  ggsave(filename = file.path(output_dir, paste0("correlation_two_panel_detrended_", sp, ".png")),
         plot = combined_corr_sp_det, width = 15, height = 7, dpi = 300)
  
  # --- Correlation Matrix for Species (Detrended)
  png(filename = file.path(output_dir, paste0("correlation_matrix_detrended_", sp, ".png")), 
      width = 800, height = 600)
  cor_mat_sp_det <- cor(df_sp_det[, c("avg_quantile", "avg_psi", "avg_tdiff")], use = "complete.obs")
  corrplot(cor_mat_sp_det, 
           method = "color", 
           type = "upper", 
           order = "hclust",
           addCoef.col = "white",   # White coefficient labels
           number.cex = 2,          # Coefficient text size set to 2
           tl.col = "black", 
           tl.srt = 45, 
           diag = FALSE,
           title = paste("Correlation Matrix (Detrended) -", sp),
           mar = c(0,0,1,0))
  dev.off()
  
  # --- ACF Plots for Species (Detrended)
  png(filename = file.path(output_dir, paste0("acf_plots_detrended_", sp, ".png")), width = 1200, height = 400)
  par(mfrow = c(1, 3))
  acf(df_sp_det$avg_quantile, main = paste("ACF:", sp, "avg_quantile (Detrended)"))
  acf(df_sp_det$avg_psi, main = paste("ACF:", sp, "avg_psi (Detrended)"))
  acf(df_sp_det$avg_tdiff, main = paste("ACF:", sp, "avg_tdiff (Detrended)"))
  par(mfrow = c(1, 1))
  dev.off()
  
  # --- Power Spectrum for Species (Detrended)
  ps_quantile_sp_det <- compute_power_spectrum_gg(df_sp_det$avg_quantile, paste("avg_quantile (Detrended) -", sp))
  ps_psi_sp_det      <- compute_power_spectrum_gg(df_sp_det$avg_psi, paste("avg_psi (Detrended) -", sp))
  ps_tdiff_sp_det    <- compute_power_spectrum_gg(df_sp_det$avg_tdiff, paste("avg_tdiff (Detrended) -", sp))
  combined_power_sp_det <- ps_quantile_sp_det$plot / ps_psi_sp_det$plot / ps_tdiff_sp_det$plot
  ggsave(filename = file.path(output_dir, paste0("power_spectrum_detrended_", sp, ".png")),
         plot = combined_power_sp_det, width = 8, height = 12, dpi = 300)
}

## 2E. Combined Species-Level Detrended Correlation Matrix
# Save a combined figure with four species (ordered as oak, beech, spruce, pine)
species_order <- c("oak", "beech", "spruce", "pine")
png(filename = file.path(output_dir, "species_correlation_matrix_detrended.png"), 
    width = 1600, height = 1200)
par(mfrow = c(2, 2), bg = "white")
for (sp in species_order) {
  sp_data_det <- subset(df_species_det, tolower(species) == sp)
  if(nrow(sp_data_det) == 0) next
  cor_matrix_sp_det <- cor(sp_data_det[, c("avg_quantile", "avg_psi", "avg_tdiff")], use = "complete.obs")
  corrplot(cor_matrix_sp_det, 
           method = "color", 
           type = "upper", 
           order = "hclust",
           addCoef.col = "white",
           number.cex = 2,
           tl.col = "black",
           tl.srt = 0,
           diag = FALSE,
           main = paste("Correlation Matrix (Detrended) -", sp),
           mar = c(0, 0, 1, 0))
}
par(mfrow = c(1, 1))
dev.off()

# ------------------------------------------------------------------------------
# End of Script
# ------------------------------------------------------------------------------
