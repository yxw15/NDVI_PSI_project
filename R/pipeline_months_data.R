# -----------------------------------------------------------------------------
# Title: NDVI, PSI, and TDiff Multi-Depth Processing Script
# Description:
#   This script processes NDVI quantiles, soil water potential (PSI), and temperature
#   difference (TDiff) data for multiple tree species across defined months and depths.
#   Depths processed: 50 cm, 100 cm, and 150 cm.
#   Outputs include:
#     - Raster exports for PSI and TDiff per species, month, and depth (in results_monthly_<depth>/)
#     - Combined monthly RData files per depth (in results/Data/)
#     - Final standardized dataframes per depth containing columns:
#       year, month, species, soil_water_potential, transpiration_deficit,
#       quantiles, x, y.
# Author: Your Name
# Date: YYYY-MM-DD
# -----------------------------------------------------------------------------

# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Source required R functions
source("R/functions_NDVI.R")
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

# Define global parameters
start_date <- "2003-01-01"
years <- 2003:2024

# Month-specific settings
months_config <- list(
  April  = list(month_day = "04-23", NDVI = "../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(month_day = "05-25", NDVI = "../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(month_day = "06-26", NDVI = "../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(month_day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"), 
  August = list(month_day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

# Species-specific settings
species_config <- list(
  Beech = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Beech_mask.tif"
  ),
  Oak   = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Oak_mask.tif"
  ),
  Spruce= list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Spruce_mask.tif"
  ),
  Pine  = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Pine_mask.tif"
  )
)

# Utility: ensure directory exists
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Function: process_species
#   Processes PSI and TDiff rasters and combines with NDVI per species/month/depth
# -----------------------------------------------------------------------------
process_species <- function(species_name, species_paths, NDVI_file, month_day, depth_val, month_name) {
  message(sprintf("[%s | Depth=%d] -> %s: Start", month_name, depth_val, species_name))
  
  # Prepare output directory: results_monthly_<depth>/<month>/<species>/
  species_paths$output_dir <- file.path(
    sprintf("results_monthly_%d", depth_val), month_name, species_name
  )
  ensure_directory(species_paths$output_dir)
  
  # 1) PSI: NetCDF -> dataframe -> save raster
  psi_df <- tryCatch({
    transfer_psi_to_df(species_paths$psi_nc, start_date)
  }, error = function(e) stop(sprintf("PSI conversion error (%s): %s", species_name, e$message)))
  save_psi_raster(psi_df, month_day, depth_val, species_paths$output_dir)
  message("   - PSI raster saved")
  
  # 2) TDiff: NetCDF -> dataframe -> save raster
  tdiff_df <- tryCatch({
    transfer_tdiff_to_df(species_paths$tdiff_nc, start_date)
  }, error = function(e) stop(sprintf("TDiff conversion error (%s): %s", species_name, e$message)))
  save_tdiff_raster(tdiff_df, month_day, depth_val, species_paths$output_dir)
  message("   - TDiff raster saved")
  
  # 3) Combine yearly NDVI, PSI, and TDiff into a dataframe
  species_results <- lapply(years, function(yr) {
    process_NDVI_PSI_TDiff_species_year(
      NDVI_file, species_name, species_paths$mask,
      yr, species_paths$output_dir
    )
  })
  species_df <- bind_rows(species_results) %>%
    mutate(depth = depth_val)
  
  # Save combined species dataframe
  save_path <- file.path(
    species_paths$output_dir,
    sprintf("NDVI_PSI_TDiff_%s_depth%d.RData", gsub("-", "", month_day), depth_val)
  )
  save(species_df, file = save_path)
  message(sprintf("   - Combined data saved: %s", basename(save_path)))
  
  return(species_df)
}

# -----------------------------------------------------------------------------
# Function: process_month
#   Processes all species for one month at a given depth
# -----------------------------------------------------------------------------
process_month <- function(month_name, depth_val) {
  message(sprintf("--- [%s | Depth=%d] Month start ---", month_name, depth_val))
  cfg <- months_config[[month_name]]
  
  # Process each species
  month_list <- lapply(names(species_config), function(species) {
    sp <- species_config[[species]]
    process_species(
      species, sp, cfg$NDVI, cfg$month_day,
      depth_val, month_name
    )
  })
  
  # Combine all species into one month dataframe
  month_df <- bind_rows(month_list) %>% mutate(month = month_name)
  
  # Save monthly combined data
  out_dir <- file.path("results", "Data", sprintf("depth_%d", depth_val))
  ensure_directory(out_dir)
  out_file <- file.path(
    out_dir,
    sprintf("AllSpecies_%s_depth%d.RData", month_name, depth_val)
  )
  save(month_df, file = out_file)
  message(sprintf("--- [%s | Depth=%d] Month data saved: %s ---", month_name, depth_val, basename(out_file)))
  
  return(month_df)
}

# -----------------------------------------------------------------------------
# Main loop: iterate over depths and months
# -----------------------------------------------------------------------------
depths <- c(100)
for (d in depths) {
  message(sprintf("===== Full run for depth = %d cm: START =====", d))
  
  # Process each month
  all_months <- lapply(names(months_config), process_month, depth_val = d)
  
  # Combine all months for this depth
  combined <- bind_rows(all_months) %>% mutate(depth = d)
  
  # Save raw combined data
  raw_dir <- file.path("results", "Data")
  ensure_directory(raw_dir)
  raw_file <- file.path(raw_dir, sprintf("AllSpecies_AllMonths_depth%d.RData", d))
  save(combined, file = raw_file)
  message(sprintf("   -> Raw combined for depth %d saved: %s", d, basename(raw_file)))
  
  # Standardize & select final columns
  final_df <- combined %>%
    # rename(
    #   soil_water_potential = psi,        # adjust 'psi' to actual PSI column name
    #   transpiration_deficit  = tdiff,     # adjust 'tdiff' to actual TDiff column name
    #   quantiles              = ndvi_quantile  # adjust to your NDVI quantile column name
    # ) %>%
    select(year, month, species,
           soil_water_potential,
           transpiration_deficit,
           Quantiles, x, y)
  
  # Assign and save the standardized dataframe
  assign(sprintf("final_df_depth%d", d), final_df)
  final_file <- file.path(raw_dir, sprintf("final_df_depth%d.RData", d))
  save(list = sprintf("final_df_depth%d", d), file = final_file)
  message(sprintf("   -> Standardized final_df_depth%d saved: %s", d, basename(final_file)))
  message(sprintf("===== Full run for depth = %d cm: COMPLETE =====\n", d))
}

message("=== All depths processed successfully! 🎉 ===")
