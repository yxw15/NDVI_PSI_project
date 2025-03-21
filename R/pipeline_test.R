# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

source("R/functions_NDVI.R")
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

# Define input files and parameters
input_file <- "../TreeSpeciesGermany/TreeSpeciesMODIS.tif"
output_folder <- "results/Species_Maps"
NDVI <- "../WZMAllDOYs/Proportions_209.nc"
start_date <- "2003-01-01"
month_day <- "07-28"
depth <- 50
years <- 2003:2024

# Helper function to create directories if missing
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

# Ensure output directories exist
ensure_directory(output_folder)

# Define species file paths
species_list <- list(
  Beech = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Beech_mask.tif",
    output_dir = "results/Beech"
  ),
  Oak = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Oak_mask.tif",
    output_dir = "results/Oak"
  ),
  Spruce = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Spruce_mask.tif",
    output_dir = "results/Spruce"
  ),
  Pine = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Pine_mask.tif",
    output_dir = "results/Pine"
  )
)

# Function to process a single species
process_species <- function(species_name, paths) {
  cat(paste("Processing", species_name, "...\n"))
  
  ensure_directory(paths$output_dir)
  
  # Convert PSI and TDiff NetCDFs to rasters
  psi_df <- transfer_psi_to_df(paths$psi_nc, start_date)
  save_psi_raster(psi_df, month_day, depth, paths$output_dir)

  tdiff_df <- transfer_tdiff_to_df(paths$tdiff_nc, start_date)
  save_tdiff_raster(tdiff_df, month_day, paths$output_dir)

  # Process NDVI, PSI, TDiff for each year
  species_results <- lapply(years, function(year) {
    process_NDVI_PSI_TDiff_species_year(NDVI, species_name, paths$mask, year, paths$output_dir)
  })
  
  # Combine yearly results
  species_df <- bind_rows(species_results)
  
  # Save data
  save_path <- file.path(paths$output_dir, "Proportions_PSI_TDiff_df.RData")
  save(species_df, file = save_path)
  cat(paste("Saved data for", species_name, "to", save_path, "\n"))
  
  return(species_df)
}

# Process all species
all_results_df <- bind_rows(lapply(names(species_list), function(species) {
  process_species(species, species_list[[species]])
}))

all_results_df <- rbind(Oak, Beech, Spruce, Pine)

# Save combined dataset
output_data_dir <- "results/Data"
ensure_directory(output_data_dir)

final_results_path <- file.path(output_data_dir, "All_Species_Proportions_PSI_TDiff.RData")
save(all_results_df, file = final_results_path)
cat("Final dataset saved.\n")

# Generate plots if data exists
if (nrow(all_results_df) > 0) {
  filtered_data <- all_results_df %>% filter(species %in% c("Beech", "Oak", "Spruce", "Pine"))
  # filtered_drought_years <- all_results_df %>% filter(year %in% c(2003, 2015, 2018, 2019, 2020, 2022))
  
  # Save plots
  plot_functions <- list(
    # plot_series_NDVI_PSI_same_month,
    # plot_species_NDVI_PSI,
    # plot_correlation_NDVI_PSI_avg,
    # plot_correlation_NDVI_PSI_species,
    # plot_mean_box_NDVI_PSI_with_slope,
    # plot_mean_box_NDVI_PSI_with_slope_bin,
    # plot_mean_box_TDiff_PSI_with_slope,
    # plot_mean_box_TDiff_PSI_with_slope_bin,
    # plot_box_NDVI_PSI,
    # plot_box_TDiff_PSI, 
    # plot_box_merge_TDiff_PSI,
    # plot_mean_box_TDiff_PSI,
    # plot_NDVI_PDM_PSIbin_slope,
    # plot_TDiff_PSIbin_slope, 
    plot_NDVI_PDM_TDiff_slope
  )
  
  plot_paths <- c(
    # "Proportions_PSI_same_month.png",
    # "Proportions_PSI_time_series_species.png",
    # "Proportions_PSI_correlation_avg.png",
    # "Proportions_PSI_correlation_species.png",
    # "mean_box_Proportions_PSI_with_slope.png",
    # "mean_box_Proportions_PSI_with_slope_bin.png",
    # "box_Proportion_PSI.png",
    # "mean_box_TDiff_PSI_with_slope.png",
    # "mean_box_TDiff_PSI_with_slope_bin.png",
    # "mean_box_TDiff_PSI.png",
    # "mean_box_merge_TDiff_PSI.png",
    # "mean_TDiff_PSI.png",
    # "Proportions_PSI_non_linear.png",
    # TDiff_PSIbin_non_linear.png",
    "Proportions_TDiff_non_linear.png"
  )
  
  for (i in seq_along(plot_functions)) {
    plot_functions[[i]](all_results_df, paste0("results/Figures/", plot_paths[i]))
  }
  
  cat("All plots generated successfully.\n")
} else {
  cat("No data available for plotting.\n")
}

cat("All analyses completed successfully.\n")
