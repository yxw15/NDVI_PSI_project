# =============================================================================
# Script Name: Multi-Species Pixel-wise ACF Analysis for Soil Water Potential 
#              (Dataframe Version)
#
# Script Description:
#   This script analyzes time series data for several parameters across four species 
#   (Oak, Beech, Spruce, Pine) using an aggregated dataframe. Each pixel (defined by 
#   its x and y coordinates) has yearly measurements that are first averaged and then 
#   used to compute the autocorrelation function (ACF) for that pixel's time series.
#
#   Analysis steps:
#   1. Data Preparation:
#      - Filter data for July and August and pivot it so that each pixel forms a 
#        complete time series (columns = years).
#
#   2. Computing the ACF:
#      - Compute the ACF for each pixel’s time series up to a specified maximum lag.
#
#   3. Bootstrapping Confidence Intervals:
#      - Use bootstrap resampling on each pixel’s time series to generate confidence 
#        intervals for the ACF.
#      - Print progress messages for each pixel.
#
#   4. Reshaping Data for Visualization:
#      - Reshape the ACF values and bootstrap confidence intervals into a long format.
#
#   5. Plotting:
#      - Generate a boxplot of pixel-wise ACF values for each lag with a ribbon showing 
#        the 25th and 75th percentiles.
#      - The x-axis shows only numeric lag values (excluding lag 0).
#
#   6. Saving the Output:
#      - Save each plot in "results/key_displays_July_August/".
#
# Required Libraries:
#   - dplyr, ggplot2, reshape2, tidyr
#
# Usage:
#   - Adjust working directory and file paths as necessary.
#   - Ensure the aggregated dataframe (final_df) contains columns:
#         x, y, year, mean_soil_water_potential, mean_transpiration_deficit, 
#         mean_Quantiles, species, month, etc.
# =============================================================================

# Clear workspace (optional)
rm(list = ls())

# Set working directory (adjust if needed)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

# Load required libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

# Define color palette for species
cb_palette <- c("Oak"   = "#E69F00",  # Orange
                "Beech" = "#0072B2",  # Deep blue
                "Spruce"= "#009E73",  # Bluish-green
                "Pine"  = "#F0E442")  # Yellow

# =============================================================================
# Functions for Dataframe Analysis
# =============================================================================

# Function: Pivot DataFrame to Create Time Series Matrix
pivotTimeSeries <- function(df, parameter_col) {
  message("Pivoting dataframe to time series matrix ...")
  ts_df <- df %>% 
    select(x, y, year, !!sym(parameter_col)) %>%
    arrange(x, y, year)
  
  # Pivot so that each pixel (x, y) becomes a row; columns are years.
  ts_wide <- pivot_wider(ts_df, names_from = year, values_from = !!sym(parameter_col))
  
  # Extract numeric matrix (dropping the x and y columns) and set row names.
  ts_matrix <- as.matrix(ts_wide[, -c(1,2)])
  rownames(ts_matrix) <- paste(ts_wide$x, ts_wide$y, sep = "_")
  
  # Keep only complete time series (rows with no NA).
  valid_ts_matrix <- ts_matrix[apply(ts_matrix, 1, function(x) all(!is.na(x))), ]
  message("Found ", nrow(valid_ts_matrix), " valid pixel time series out of ", nrow(ts_matrix))
  return(valid_ts_matrix)
}

# Function 2: Compute the ACF Matrix for Each Pixel's Time Series
computeACFMatrix <- function(valid_vals, lag.max) {
  message("Computing ACF for each pixel time series ...")
  acf_matrix <- t(apply(valid_vals, 1, function(ts) {
    acf(ts, plot = FALSE, lag.max = lag.max)$acf
  }))
  message("ACF matrix computed.")
  return(acf_matrix)
}

# Function 4: Bootstrap Confidence Intervals for Each Pixel's ACF
bootstrapCI <- function(valid_vals, lag.max, n_boot = 1000, alpha = 0.05) {
  message("Bootstrapping confidence intervals ...")
  set.seed(123)
  n_pixels <- nrow(valid_vals)
  boot_ci_matrix_lower <- matrix(NA, nrow = n_pixels, ncol = lag.max + 1)
  boot_ci_matrix_upper <- matrix(NA, nrow = n_pixels, ncol = lag.max + 1)
  
  for (i in 1:n_pixels) {
    message(sprintf("Bootstrapping: Processing pixel %d of %d", i, n_pixels))
    ts <- valid_vals[i, ]
    acf_boot <- replicate(n_boot, {
      ts_perm <- sample(ts)
      as.numeric(acf(ts_perm, plot = FALSE, lag.max = lag.max)$acf)
    })
    boot_ci_matrix_lower[i, ] <- apply(acf_boot, 1, quantile, probs = alpha / 2)
    boot_ci_matrix_upper[i, ] <- apply(acf_boot, 1, quantile, probs = 1 - alpha / 2)
  }
  message("Bootstrapping complete.")
  
  global_lower_bound <- apply(boot_ci_matrix_lower, 2, function(x) quantile(x, probs = alpha / 2, na.rm = TRUE))
  global_upper_bound <- apply(boot_ci_matrix_upper, 2, function(x) quantile(x, probs = 1 - alpha / 2, na.rm = TRUE))
  
  return(list(lower = boot_ci_matrix_lower, 
              upper = boot_ci_matrix_upper,
              global_lower = global_lower_bound,
              global_upper = global_upper_bound))
}

# Function 5: Reshape the ACF and Bootstrap Data for Plotting
reshapeACFData <- function(acf_matrix, boot_ci_lower, boot_ci_upper) {
  message("Reshaping data for plotting ...")
  acf_df <- melt(acf_matrix)
  lower_df <- melt(boot_ci_lower)
  upper_df <- melt(boot_ci_upper)
  
  colnames(acf_df) <- c("Pixel", "Lag", "ACF")
  colnames(lower_df) <- c("Pixel", "Lag", "Lower_CI")
  colnames(upper_df) <- c("Pixel", "Lag", "Upper_CI")
  
  acf_combined <- cbind(acf_df, Lower_CI = lower_df$Lower_CI, Upper_CI = upper_df$Upper_CI)
  # Adjust Lag index: subtract 1 so that lag 0 corresponds to LagNum 0.
  acf_combined$LagNum <- acf_combined$Lag - 1
  
  ci_summary <- acf_combined %>%
    group_by(LagNum) %>%
    summarise(
      ci_lower = quantile(Lower_CI, 0.25, na.rm = TRUE),
      ci_upper = quantile(Upper_CI, 0.75, na.rm = TRUE)
    )
  
  message("Data reshaped for plotting.")
  return(list(acf_combined = acf_combined, ci_summary = ci_summary))
}

# Function 6: Generate the ACF Plot with Bootstrapped Confidence Intervals
generateACFPlot <- function(acf_combined, ci_summary, lag.max, species, var_short, cb_palette) {
  message("Generating plot for ", species, " - ", var_short, " ...")
  # Exclude lag 0 from the plot.
  acf_combined_filtered <- acf_combined %>% dplyr::filter(LagNum > 0)
  ci_summary_filtered <- ci_summary %>% dplyr::filter(LagNum > 0)
  
  p <- ggplot(acf_combined_filtered, aes(x = LagNum, y = ACF)) +
    geom_ribbon(data = ci_summary_filtered,
                aes(x = LagNum, ymin = ci_lower, ymax = ci_upper),
                fill = cb_palette[species], alpha = 0.3, inherit.aes = FALSE) +
    geom_boxplot(aes(group = LagNum), fill = cb_palette[species], alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = sort(unique(acf_combined_filtered$LagNum)),
                       labels = sort(unique(acf_combined_filtered$LagNum))) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(title = paste("Pixel-wise ACF with Bootstrapped CI -", species, var_short),
         x = "Lag", y = "ACF") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  message("Plot generated for ", species, " - ", var_short)
  return(p)
}

# Function 7: Save the Generated Plot to File
saveACFPlot <- function(plot, species, var_short, out_dir = "results/key_displays_July_August", width = 12, height = 6) {
  message("Saving plot for ", species, " - ", var_short, " ...")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  output_file <- file.path(out_dir, paste0(species, "_", var_short, "_ACF_CI_Boxplot.png"))
  ggsave(filename = output_file, plot = plot, width = width, height = height)
  message("Plot saved to ", output_file)
}

# =============================================================================
# Main Function: Process the DataFrame for One Species and Parameter
# =============================================================================
processDataFrame <- function(df, parameter_col, lag.max = 21, n_boot = 1000, alpha = 0.05,
                             var_short, species, cb_palette, out_dir = "results/key_displays_July_August") {
  message("=== Processing for ", species, " - ", var_short, " ===")
  
  # Pivot the dataframe to create a time series matrix.
  valid_ts_matrix <- pivotTimeSeries(df, parameter_col)
  
  # Determine the number of years available and adjust lag.max accordingly.
  n_years <- ncol(valid_ts_matrix)
  lag_actual <- min(lag.max, n_years - 1)
  
  # Compute the ACF matrix using the adjusted lag.
  acf_matrix <- computeACFMatrix(valid_ts_matrix, lag_actual)
  
  # Perform bootstrap resampling to compute confidence intervals for the ACF.
  boot_ci <- bootstrapCI(valid_ts_matrix, lag_actual, n_boot, alpha)
  
  # Reshape the ACF and bootstrap CI data for plotting.
  reshaped <- reshapeACFData(acf_matrix, boot_ci$lower, boot_ci$upper)
  
  # Generate the final plot using the adjusted lag.
  p <- generateACFPlot(reshaped$acf_combined, reshaped$ci_summary, lag_actual, species, var_short, cb_palette)
  
  # Save the final plot.
  saveACFPlot(p, species, var_short, out_dir)
  
  message("=== Completed processing for ", species, " - ", var_short, " ===")
  return(list(plot = p, acf_data = reshaped$acf_combined, ci_summary = reshaped$ci_summary))
}

# =============================================================================
# Execution: Loop Over Species and Parameters
# =============================================================================

# List of species to process
species_list <- c("Oak", "Beech", "Pine", "Spruce")

# List of parameters and their short names
parameters <- c("mean_soil_water_potential", "mean_transpiration_deficit", "mean_Quantiles")
var_names <- c("mean_soil_water_potential" = "PSI",
               "mean_transpiration_deficit" = "TDiff",
               "mean_Quantiles" = "Q")





# Loop over each species
for (sp in species_list) {
  message("========== Processing species: ", sp, " ==========")
  # Filter data for current species and for July and August.
  df_sp <- final_df %>% filter(month %in% c("July", "August"), species == sp)
  
  # Aggregate data by pixel (x, y) and year.
  df_mean_sp <- df_sp %>%
    group_by(x, y, year) %>%
    summarise(
      mean_Quantiles = mean(Quantiles, na.rm = TRUE),
      mean_soil_water_potential = mean(soil_water_potential, na.rm = TRUE),
      mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE),
      species = first(species)
    ) %>% 
    ungroup()
  
  # Loop over each parameter for the current species.
  for (param in parameters) {
    message("----- Processing parameter: ", param, " for species: ", sp, " -----")
    processDataFrame(
      df = df_mean_sp,
      parameter_col = param,
      lag.max = 22,
      n_boot = 1000,
      alpha = 0.05,
      var_short = var_names[param],
      species = sp,
      cb_palette = cb_palette,
      out_dir = "results/key_displays_July_August"
    )
  }
}
