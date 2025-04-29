# -----------------------------------------------------------------------------
# Title: NDVI, PSI, and TDiff Multi-Depth Processing Script (Parallelized) 🚀
# Description: Adds parallel processing with mclapply and fun stickers in messages 🎉
# -----------------------------------------------------------------------------

# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)  # for mclapply

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Source required R functions
source("R/functions_NDVI.R")
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

# Define global parameters
start_date <- "2003-01-01"
years <- 2003:2024

# Detect number of cores and reserve one for the OS
n_cores <- max(1, detectCores() - 1)
message(sprintf("🖥️ Detected %d cores, using %d for parallel tasks 😎", detectCores(), n_cores))

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
  Oak = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Oak_mask.tif"
  ),
  Spruce = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Spruce_mask.tif"
  ),
  Pine = list(
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
#   Processes PSI & TDiff rasters, combines with NDVI, with skip logic and stickers 🌿
# -----------------------------------------------------------------------------
process_species <- function(species_name, species_paths, NDVI_file, month_day, depth_val, month_name) {
  message(sprintf("🌱 [%s | Depth=%dcm] -> %s: Starting work...", month_name, depth_val, species_name))
  
  # Prepare output
  output_dir <- file.path(sprintf("results_monthly_%d", depth_val), month_name, species_name)
  ensure_directory(output_dir)
  save_path <- file.path(output_dir, sprintf("NDVI_PSI_TDiff_%s_depth%d.RData", gsub("-", "", month_day), depth_val))
  
  # Skip if exists
  if (file.exists(save_path)) {
    message(sprintf("⏭️ Skipping %s (already done)! 👍", species_name))
    load(save_path)  # loads species_df
    return(species_df)
  }
  
  # PSI conversion
  psi_df <- tryCatch({
    transfer_psi_to_df(species_paths$psi_nc, start_date)
  }, error = function(e) stop(sprintf("🔴 PSI conversion error for %s: %s", species_name, e$message)))
  save_psi_raster(psi_df, month_day, depth_val, output_dir)
  message("✅ PSI raster saved 🎯")
  
  # TDiff conversion
  tdiff_df <- tryCatch({
    transfer_tdiff_to_df(species_paths$tdiff_nc, start_date)
  }, error = function(e) stop(sprintf("🔴 TDiff conversion error for %s: %s", species_name, e$message)))
  save_tdiff_raster(tdiff_df, month_day, depth_val, output_dir)
  message("✅ TDiff raster saved 🎯")
  
  # Combine yearly data in parallel
  message(sprintf("🔄 Combining NDVI, PSI & TDiff for %s across %d years...", species_name, length(years)))
  species_results <- mclapply(years, function(yr) {
    process_NDVI_PSI_TDiff_species_year(
      NDVI_file, species_name, species_paths$mask,
      yr, output_dir
    )
  }, mc.cores = n_cores)
  species_df <- bind_rows(species_results) %>% mutate(depth = depth_val)
  
  # Save result
  save(species_df, file = save_path)
  message(sprintf("🎉 Combined data saved for %s: %s", species_name, basename(save_path)))
  return(species_df)
}

# -----------------------------------------------------------------------------
# Function: process_month
#   Processes all species for a given month in parallel 🗓️
# -----------------------------------------------------------------------------
process_month <- function(month_name, depth_val) {
  cfg <- months_config[[month_name]]
  message(sprintf("🌸 [%s | Depth=%dcm] Month start 🌸", month_name, depth_val))
  
  # Prepare monthly output
  out_dir <- file.path("results", "Data", sprintf("depth_%d", depth_val))
  ensure_directory(out_dir)
  out_file <- file.path(out_dir, sprintf("AllSpecies_%s_depth%d.RData", month_name, depth_val))
  
  # Skip if exists
  if (file.exists(out_file)) {
    message(sprintf("⏭️ Skipping month %s (already done)! 👍", month_name))
    load(out_file)  # loads month_df
    return(month_df)
  }
  
  # Parallel species processing
  message("🚀 Launching species processing in parallel...")
  month_list <- mclapply(names(species_config), function(species) {
    sp <- species_config[[species]]
    process_species(species, sp, cfg$NDVI, cfg$month_day, depth_val, month_name)
  }, mc.cores = n_cores)
  
  month_df <- bind_rows(month_list) %>% mutate(month = month_name)
  save(month_df, file = out_file)
  message(sprintf("🎉 Month data saved: %s", basename(out_file)))
  return(month_df)
}

# -----------------------------------------------------------------------------
# Main loop: iterate over depths and months in parallel 🌎
# -----------------------------------------------------------------------------
depths <- c(150)

for (d in depths) {
  message(sprintf("🧭 ===== Depth = %d cm RUN START ===== 🧭", d))
  
  # Raw combined filename
  raw_dir <- file.path("results", "Data")
  ensure_directory(raw_dir)
  raw_file <- file.path(raw_dir, sprintf("AllSpecies_AllMonths_depth%d.RData", d))
  
  if (file.exists(raw_file)) {
    message(sprintf("⏭️ Skipping full combined for depth %d (already exists)! 👍", d))
    load(raw_file)  # loads 'combined'
  } else {
    # Parallel month processing
    message("🚀 Starting month-level parallel processing...")
    all_months <- mclapply(names(months_config), process_month, depth_val = d, mc.cores = n_cores)
    combined <- bind_rows(all_months) %>% mutate(depth = d)
    save(combined, file = raw_file)
    message(sprintf("🎉 Raw combined saved: %s", basename(raw_file)))
  }
  
  # Standardize final dataframe
  final_df <- combined %>%
    rename(
      soil_water_potential = psi,
      transpiration_deficit = tdiff,
      quantiles = ndvi_quantile
    ) %>%
    select(year, month, species,
           soil_water_potential,
           transpiration_deficit,
           quantiles, x, y)
  
  # Assign and save
  assign(sprintf("final_df_depth%d", d), final_df)
  final_file <- file.path(raw_dir, sprintf("final_df_depth%d.RData", d))
  save(list = sprintf("final_df_depth%d", d), file = final_file)
  message(sprintf("🥳 Standardized final_df_depth%d saved: %s", d, basename(final_file)))
  
  message(sprintf("🧭 ===== Depth = %d cm RUN COMPLETE ===== 🧭\n", d))
}

message("🏁 ALL DEPTHS PROCESSED SUCCESSFULLY! 🏁")
