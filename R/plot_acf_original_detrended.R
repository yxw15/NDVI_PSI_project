#########################################################
# NDVI_PSI Project - Temporal Auto-correlation
# Plot ACF box plots for original and detrended dataframe
# Author: Yixuan Wang (adapted)
# Date: Sys.Date()
#########################################################

#### 1. SETUP & DATA LOADING ####

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Define output folders
output_folder_raster <- "results_temporal/FFT/Raster"
output_folder_figures <- "results_temporal/FFT/Figures"
output_folder_data   <- "results/Data"

# Create output folders if they don't exist
dirs <- c(output_folder_raster, output_folder_figures, output_folder_data)
for(d in dirs){
  if (!dir.exists(d)) { dir.create(d, recursive = TRUE) }
}

# Load required packages
library(dplyr)
library(tidyr)
library(terra)
library(rlang)
library(stats)
library(ggplot2)
library(viridis)
library(grid)  # for unit()


#### 2. FUNCTION DEFINITIONS ####

# Function to generate ACF plots for a given dataset
generate_acf_plots <- function(data, data_type, species_order, vars, output_dir, cb_palette) {
  # Ensure species factor ordering
  data$species <- factor(data$species, levels = species_order)
  
  # Loop through each variable
  for (var in vars) {
    rast_list_species <- list()
    
    # Build a raster stack for each species
    for (sp in species_order) {
      sp_df <- data %>% filter(species == sp)
      sp_df$year <- as.character(sp_df$year)
      yrs <- sort(unique(sp_df$year))
      
      rast_list <- lapply(yrs, function(yr) {
        sub_df <- sp_df %>% filter(year == yr)
        rast(sub_df[, c("x", "y", var)], type = "xyz")
      })
      names(rast_list) <- yrs
      sp_rast <- rast(rast_list)
      rast_list_species[[sp]] <- sp_rast
    }
    
    # A. Random Pixel ACF Plot
    output_file_rand <- file.path(output_dir, paste0(var, "_ACF_one_pixel_subpanel_", data_type, ".png"))
    png(filename = output_file_rand, width = 1200, height = 800)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    for (sp in species_order) {
      sp_rast <- rast_list_species[[sp]]
      raster_values <- values(sp_rast)
      valid_cells <- which(apply(raster_values, 1, function(x) all(!is.na(x))))
      if (length(valid_cells) == 0) {
        plot.new()
        title(main = paste(sp, var, "\nNo valid pixel data"), col.main = cb_palette[sp])
        next
      }
      cell_idx <- sample(valid_cells, 1)
      ts_values <- raster_values[cell_idx, ]
      acf(ts_values, main = paste(sp, var, "ACF for Cell", cell_idx), col.main = cb_palette[sp])
    }
    dev.off()
    
    # B. Boxplot of ACF for All Pixels
    output_file_box <- file.path(output_dir, paste0(var, "_ACF_all_pixels_boxplot_subpanel_", data_type, ".png"))
    png(filename = output_file_box, width = 1200, height = 800)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    for (sp in species_order) {
      sp_rast <- rast_list_species[[sp]]
      vals <- values(sp_rast)
      valid_pixels <- apply(vals, 1, function(x) all(!is.na(x)))
      valid_vals <- vals[valid_pixels, ]
      n_years <- ncol(valid_vals)
      lag_max <- n_years - 1
      acf_matrix <- t(apply(valid_vals, 1, function(ts) {
        acf_res <- acf(ts, plot = FALSE, lag.max = lag_max)
        as.vector(acf_res$acf)
      }))
      lag_names <- paste0("lag_", 0:lag_max)
      acf_df <- as.data.frame(acf_matrix)
      colnames(acf_df) <- lag_names
      ci_bound <- 1.96 / sqrt(n_years)
      boxplot(acf_df,
              main = paste("Distribution of", sp, var, "ACF"),
              xlab = "Lag",
              ylab = "ACF",
              col = cb_palette[sp])
      abline(h = ci_bound, col = "blue", lty = 2)
      abline(h = -ci_bound, col = "blue", lty = 2)
      abline(h = 0, col = "red", lty = 2)
    }
    dev.off()
  }
}


#### 3. CONFIGURATION ####

# Define species order and analysis years
species_order <- c("Oak", "Beech", "Spruce", "Pine")
years <- as.character(2003:2024)  # adjust as needed

# Define color palette
cb_palette <- c("Oak"    = "#E69F00",  # Orange
                "Beech"  = "#0072B2",  # Deep blue
                "Spruce" = "#009E73",  # Bluish-green
                "Pine"   = "#F0E442")  # Yellowish

# Variables to process
vars <- c("Quantiles", "soil_water_potential", "transpiration_deficit")

vars <- c("Proportions", "soil_water_potential", "transpiration_deficit")

# Define common output directory for results
output_dir <- "results_temporal/FFT/results"
if (!dir.exists(output_dir)) { 
  dir.create(output_dir, recursive = TRUE)
}


#### 4. PROCESS ORIGINAL DATA ####

# Load the combined data for original analysis
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")  # loads all_results_df

load("results/Data/All_Species_Proportions_PSI_TDiff.RData") 

# Generate ACF plots for original data
generate_acf_plots(data = all_results_df,
                   data_type = "original",
                   species_order = species_order,
                   vars = vars,
                   output_dir = output_dir,
                   cb_palette = cb_palette)


#### 5. PROCESS DETRENDED DATA ####

# Load the detrended data for FFT analysis
load("results/Data/All_Quantiles_PSI_TDiff_species_year_fft_detrended.RData")  # loads final_df_clean

# Generate ACF plots for detrended data
generate_acf_plots(data = final_df_clean,
                   data_type = "fft_detrended",
                   species_order = species_order,
                   vars = vars,
                   output_dir = output_dir,
                   cb_palette = cb_palette)
