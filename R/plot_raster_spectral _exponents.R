###############################################################################
# Description:
# This script processes raster data of spectral exponents for four species
# (Oak, Beech, Spruce, Pine) and three parameters (Quantiles, soil_water_potential,
# transpiration_deficit). For each species and parameter combination, it loads both
# original and detrended spectral exponent rasters, computes their difference, and
# reshapes the data for plotting.
#
# The script bins the spectral exponent values into four noise categories based on:
# - "Red noise": slope < -2
# - "Pink noise": slope between -2 and -1
# - "White noise": slope between -1 and 0
# - "High-frequency": slope >= 0
#
# A composite plot is then created (faceted by parameter and method) where the points
# are colored by these noise categories, and the Germany boundary is overlaid on each
# panel. Both the boundary and the point data are converted to sf objects so that they
# can be used with coord_sf().
#
# NOTE: A custom panel border has been removed from the theme to avoid the error:
# "no applicable method for 'depth' applied to an object of class 'NULL'".
#
# The resulting plots are saved in the specified output directory.
###############################################################################

library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(grid)       # for unit()

# Additional libraries for spatial boundaries and conversion
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Define species, parameters, and output directory
species_list <- c("Oak", "Beech", "Spruce", "Pine")
params <- c("Quantiles", "soil_water_potential", "transpiration_deficit")
output_dir <- "results_temporal/FFT/results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load Germany boundary as an sf object (ensure CRS compatibility, typically WGS84)
germany <- ne_countries(country = "Germany", scale = "large", returnclass = "sf")

# Function to create a long-format data frame for one species and one parameter
get_species_param_df <- function(species, param) {
  # Construct file paths for original and detrended spectral exponent rasters
  orig_file <- file.path("results_temporal/FFT/Raster", 
                         paste0(species, "_", param, "_original_spectral_exponent.tif"))
  det_file  <- file.path("results_temporal/FFT/Raster", 
                         paste0(species, "_", param, "_fft_detrended_spectral_exponent.tif"))
  
  # Load the rasters
  rast_orig <- rast(orig_file)
  rast_det  <- rast(det_file)
  
  # Convert to data frames with x, y coordinates
  df_orig <- as.data.frame(rast_orig, xy = TRUE)
  df_det  <- as.data.frame(rast_det, xy = TRUE)
  
  # Rename value columns to "original" and "detrended"
  colnames(df_orig)[3] <- "original"
  colnames(df_det)[3]  <- "detrended"
  
  # Merge by x and y coordinates
  df_merge <- inner_join(df_orig, df_det, by = c("x", "y"))
  
  # Compute difference (Original - Detrended)
  df_merge <- df_merge %>% mutate(diff = original - detrended)
  
  # Add columns for the parameter and species name
  df_merge <- df_merge %>% mutate(param = param, species = species)
  
  # Pivot longer so that one column contains method (Original, Detrended, Difference)
  df_long <- df_merge %>%
    pivot_longer(cols = c("original", "detrended", "diff"),
                 names_to = "method",
                 values_to = "value")
  
  df_long
}

# Loop over species and parameters to combine all data into one data frame
all_data <- do.call(rbind,
                    lapply(species_list, function(sp) {
                      do.call(rbind, lapply(params, function(p) {
                        get_species_param_df(sp, p)
                      }))
                    }))

# Set factors for proper ordering
all_data$species <- factor(all_data$species, levels = species_list)
all_data$param <- factor(all_data$param, levels = params)
all_data$method <- factor(all_data$method, levels = c("original", "detrended", "diff"),
                          labels = c("Original", "Detrended", "Difference"))

# Bin the continuous values based on theory:
# - "Red noise": very strong low-frequency dominance (slope < -2)
# - "Pink noise": intermediate low-frequency (slope between -2 and -1)
# - "White noise": near white noise (slope between -1 and 0)
# - "High-frequency": positive slopes (>= 0)
all_data <- all_data %>%
  mutate(value_group = cut(value,
                           breaks = c(-Inf, -2, -1, 0, Inf),
                           labels = c("Red noise", "Pink noise", "White noise", "High-frequency")))

# Generate four colors from the magma palette (magma is color-blind friendly)
magma_colors <- magma(4)

# Loop over species to create and save one composite plot per species.
for(sp in species_list) {
  # Subset data for the species
  sp_data <- filter(all_data, species == sp)
  # Convert point data to an sf object (using x, y as coordinates, assuming WGS84)
  sp_data_sf <- st_as_sf(sp_data, coords = c("x", "y"), crs = 4326, remove = FALSE)
  
  # Create composite plot: rows = parameters, columns = methods.
  # Both the Germany boundary and the point data are plotted using geom_sf().
  p <- ggplot() +
    geom_sf(data = germany, fill = NA, color = "black", size = 0.5) +
    geom_sf(data = sp_data_sf, aes(color = value_group), size = 0.5) +
    facet_grid(param ~ method) +
    scale_color_manual(values = setNames(magma_colors,
                                         c("Red noise", "Pink noise", "White noise", "High-frequency"))) +
    coord_sf() +
    labs(title = paste(sp, "Spectral Exponent"),
         x = "Longitude", y = "Latitude", color = "Spectral Exponent Category") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.key.size = unit(1.5, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
      # Removed panel.border to avoid "depth" error with coord_sf()
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Print and save the plot for the species
  print(p)
  output_file <- file.path(output_dir, paste0(sp, "_Composite_Spectral_Exponent_Binned.png"))
  ggsave(output_file, p, width = 15, height = 12, dpi = 300)
}
