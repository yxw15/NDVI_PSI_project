# =============================================================================
# Script Description:
#
# This script is designed to process and analyze raster data for different tree species
# based on various parameters. The workflow consists of several key steps:
#
# 1. Loading Raster Data:
#    - The script loads a stack of raster mask files for a given species and variable 
#      (e.g., soil water potential, transpiration deficit, or NDVI quantiles) from a specified 
#      directory structure.
#
# 2. Extracting Valid Pixel Values:
#    - It extracts pixel values from the raster stack and filters out those with any NA 
#      values across the time series, ensuring that only complete time series are analyzed.
#
# 3. Computing the Autocorrelation Function (ACF):
#    - For each pixel's time series, the script computes the ACF up to a specified lag.
#
# 4. Bootstrapping Confidence Intervals:
#    - It performs bootstrap resampling on the time series to generate confidence intervals 
#      for the ACF estimates, providing a measure of uncertainty.
#
# 5. Reshaping Data for Visualization:
#    - The ACF values and their corresponding bootstrap confidence intervals are reshaped into 
#      a format suitable for plotting.
#
# 6. Generating ACF Plots:
#    - The script creates detailed ACF plots using ggplot2, overlaying boxplots of the ACF values 
#      with a ribbon representing the bootstrapped confidence intervals.
#
# 7. Saving the Output:
#    - Finally, the generated plots are saved to a specified output directory.
#
# The script is structured to loop through multiple species (e.g., Oak, Beech, Spruce, Pine) and
# various parameters (e.g., soil water potential, transpiration deficit, NDVI quantiles), making it 
# a flexible tool for spatial and temporal analysis of raster datasets.
#
# Required Libraries:
#    - terra: for handling raster data
#    - ggplot2: for creating visualizations
#    - reshape2: for data reshaping operations
#    - dplyr: for data manipulation and summarization
#
# Usage:
#    - Adjust the working directory and paths as needed.
#    - Define the list of species and parameters.
#    - Run the main loop to process each combination and generate the corresponding ACF plots.
#
# =============================================================================

# Load required packages
library(terra)
library(ggplot2)
library(reshape2)
library(dplyr)

# Function 1: Load the raster stack for a given species and variable
loadRasterData <- function(species, var_short, raster_base_path = "results") {
  message("Loading raster files for ", species, " - ", var_short, " ...")
  raster_folder <- file.path(raster_base_path, species, paste0(var_short, "_mask"))
  tif_files <- list.files(raster_folder, pattern = "^.*_\\d{4}_mask\\.tif$", full.names = TRUE)
  if (length(tif_files) == 0) stop("No files found in ", raster_folder)
  
  years <- as.numeric(gsub("^.*_(\\d{4})_mask\\.tif$", "\\1", tif_files))
  tif_files <- tif_files[order(years)]
  raster_stack <- rast(tif_files)
  names(raster_stack) <- as.character(sort(years))
  message("Loaded ", length(tif_files), " raster files for ", species, " - ", var_short)
  return(raster_stack)
}

# Function 2: Extract pixel values and filter for valid pixels (non-NA across time)
extractValidPixelValues <- function(raster_stack) {
  message("Extracting pixel values and filtering valid pixels ...")
  vals <- values(raster_stack)
  valid_pixels <- apply(vals, 1, function(x) all(!is.na(x)))
  valid_vals <- vals[valid_pixels, ]
  message("Found ", sum(valid_pixels), " valid pixels out of ", nrow(vals))
  return(valid_vals)
}

# Function 3: Compute the ACF matrix for each pixel time series
computeACFMatrix <- function(valid_vals, lag.max = 12) {
  message("Computing ACF for each pixel time series ...")
  acf_matrix <- t(apply(valid_vals, 1, function(ts) {
    acf(ts, plot = FALSE, lag.max = lag.max)$acf
  }))
  message("ACF matrix computed.")
  return(acf_matrix)
}

# Function 4: Bootstrap confidence intervals for each pixel's ACF
bootstrapCI <- function(valid_vals, lag.max = 12, n_boot = 1000, alpha = 0.05) {
  message("Bootstrapping confidence intervals ...")
  set.seed(123)
  n_pixels <- nrow(valid_vals)
  boot_ci_matrix_lower <- matrix(NA, nrow = n_pixels, ncol = lag.max + 1)
  boot_ci_matrix_upper <- matrix(NA, nrow = n_pixels, ncol = lag.max + 1)
  
  for (i in 1:n_pixels) {
    ts <- valid_vals[i, ]
    acf_boot <- replicate(n_boot, {
      ts_perm <- sample(ts)
      acf(ts_perm, plot = FALSE, lag.max = lag.max)$acf
    })
    boot_ci_matrix_lower[i, ] <- apply(acf_boot, 1, quantile, probs = alpha / 2)
    boot_ci_matrix_upper[i, ] <- apply(acf_boot, 1, quantile, probs = 1 - alpha / 2)
  }
  message("Bootstrapping complete.")
  return(list(lower = boot_ci_matrix_lower, upper = boot_ci_matrix_upper))
}

# Function 5: Reshape the ACF and bootstrap data for plotting
reshapeACFData <- function(acf_matrix, boot_ci_lower, boot_ci_upper) {
  message("Reshaping data for plotting ...")
  acf_df <- melt(acf_matrix)
  lower_df <- melt(boot_ci_lower)
  upper_df <- melt(boot_ci_upper)
  
  colnames(acf_df) <- c("Pixel", "Lag", "ACF")
  colnames(lower_df) <- c("Pixel", "Lag", "Lower_CI")
  colnames(upper_df) <- c("Pixel", "Lag", "Upper_CI")
  
  acf_combined <- cbind(acf_df, Lower_CI = lower_df$Lower_CI, Upper_CI = upper_df$Upper_CI)
  acf_combined$LagNum <- acf_combined$Lag - 1
  
  # Summarize the bootstrap confidence intervals for the ribbon
  ci_summary <- acf_combined %>%
    group_by(LagNum) %>%
    summarise(
      ci_lower = quantile(Lower_CI, 0.25, na.rm = TRUE),
      ci_upper = quantile(Upper_CI, 0.75, na.rm = TRUE)
    )
  
  message("Data reshaped for plotting.")
  return(list(acf_combined = acf_combined, ci_summary = ci_summary))
}

# Function 6: Generate the ACF plot with bootstrapped confidence intervals
generateACFPlot <- function(acf_combined, ci_summary, lag.max, species, var_short, cb_palette) {
  message("Generating plot for ", species, " - ", var_short, " ...")
  p <- ggplot(acf_combined, aes(x = LagNum, y = ACF)) +
    geom_ribbon(data = ci_summary,
                aes(x = LagNum, ymin = ci_lower, ymax = ci_upper),
                fill = cb_palette[species], alpha = 0.3, inherit.aes = FALSE) +
    geom_boxplot(aes(group = LagNum), fill = cb_palette[species], alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = 0:lag.max, labels = paste0("Lag_", 0:lag.max)) +
    labs(title = paste("Pixel-wise ACF with Bootstrapped CI -", species, var_short),
         x = "Lag", y = "ACF") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  message("Plot generated for ", species, " - ", var_short)
  return(p)
}

# Function 7: Save the generated plot to file
saveACFPlot <- function(plot, species, var_short, out_dir = "results_acf", width = 12, height = 6) {
  message("Saving plot for ", species, " - ", var_short, " ...")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  output_file <- file.path(out_dir, paste0(species, "_", var_short, "_ACF_CI_Boxplot.png"))
  ggsave(filename = output_file, plot = plot, width = width, height = height)
  message("Plot saved to ", output_file)
}

# Main function: Process the raster data for one species and parameter
processRasterData <- function(species, var, var_names, cb_palette,
                              raster_base_path = "results",
                              out_dir = "results_acf",
                              lag.max = 12,
                              n_boot = 1000,
                              alpha = 0.05) {
  var_short <- var_names[var]
  message("=== Processing for ", species, " - ", var_short, " ===")
  
  # Step 1: Load Raster Data
  raster_stack <- loadRasterData(species, var_short, raster_base_path)
  
  # Step 2: Extract Valid Pixel Values
  valid_vals <- extractValidPixelValues(raster_stack)
  
  # Step 3: Compute ACF Matrix
  acf_matrix <- computeACFMatrix(valid_vals, lag.max)
  
  # Step 4: Bootstrap Confidence Intervals
  boot_ci <- bootstrapCI(valid_vals, lag.max, n_boot, alpha)
  
  # Step 5: Reshape Data for Plotting
  reshaped <- reshapeACFData(acf_matrix, boot_ci$lower, boot_ci$upper)
  
  # Step 6: Generate the Plot
  p <- generateACFPlot(reshaped$acf_combined, reshaped$ci_summary, lag.max, species, var_short, cb_palette)
  
  # Step 7: Save the Plot
  saveACFPlot(p, species, var_short, out_dir)
  
  message("=== Completed processing for ", species, " - ", var_short, " ===")
  
  return(list(plot = p, acf_data = reshaped$acf_combined, ci_summary = reshaped$ci_summary))
}

# ---------------------------
# Example: Looping through all species and parameters

# Set working directory (adjust if necessary)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Define species and parameters
species_all <- c("Pine")
parameter_all <- c("transpiration_deficit")

# Variable names for parameters
var_names <- c("soil_water_potential" = "PSI",
               "transpiration_deficit" = "TDiff",
               "Quantiles" = "NDVI")

# Color palette for species
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# Loop through each species and parameter combination
results_list <- list()
for (sp in species_all) {
  for (par in parameter_all) {
    message("===== Starting processing for ", sp, " with parameter: ", par, " =====")
    res <- tryCatch({
      processRasterData(species = sp, var = par, var_names = var_names, cb_palette = cb_palette,
                        raster_base_path = "results", out_dir = "results_acf",
                        lag.max = 12, n_boot = 1000, alpha = 0.05)
    }, error = function(e) {
      message("Error processing ", sp, " with parameter ", par, ": ", e$message)
      return(NULL)
    })
    results_list[[paste0(sp, "_", par)]] <- res
    message("===== Completed processing for ", sp, " with parameter: ", par, " =====")
  }
}

# Optionally, inspect results_list for further analysis
