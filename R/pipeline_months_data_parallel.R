# -----------------------------------------------------------------------------
# Title: NDVI, PSI, and TDiff Multi-Depth Processing Script (Parallelized & Robust) 🚀
# Description: Adds parallel processing with mclapply, robust missing-data handling, always rewrite final results, and fun stickers 🎉
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
years      <- 2003:2024

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
  # Explicit running message
  message(sprintf("🔍 Processing species '%s' for month '%s' at depth %dcm...", 
                  species_name, month_name, depth_val))
  # Detailed start message
  message(sprintf("🌱 [%s | Depth=%dcm] -> %s: Starting work...", 
                  month_name, depth_val, species_name))
  
  # Prepare output
  output_dir <- file.path(sprintf("results_monthly_%d", depth_val), month_name, species_name)
  ensure_directory(output_dir)
  save_path <- file.path(output_dir, sprintf("NDVI_PSI_TDiff_%s_depth%d.RData", gsub("-", "", month_day), depth_val))
  
  # Skip if exists for efficiency
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
#   Processes all species for a given month in parallel, handling missing files 🗓️
# -----------------------------------------------------------------------------
process_month <- function(month_name, depth_val) {
  # Month start message
  message(sprintf("🔖 Starting month '%s' at depth %dcm...", month_name, depth_val))
  cfg <- months_config[[month_name]]
  
  species_dfs <- lapply(names(species_config), function(species) {
    sp <- species_config[[species]]
    file <- file.path(
      sprintf("results_monthly_%d", depth_val),
      month_name, species,
      sprintf("NDVI_PSI_TDiff_%s_depth%d.RData",
              gsub("-", "", cfg$month_day), depth_val)
    )
    if (file.exists(file)) {
      message(sprintf("⏭️ [%s | %s] already done – loading…", month_name, species))
      e <- new.env()
      load(file, envir = e)
      return(e$species_df)
    } else {
      return(process_species(species, sp, cfg$NDVI, cfg$month_day, depth_val, month_name))
    }
  })
  
  month_df <- bind_rows(species_dfs) %>% mutate(month = month_name, depth = depth_val)
  
  # Save month summary (always overwrite)
  out_dir <- file.path("results", "Data", sprintf("depth_%d", depth_val))
  ensure_directory(out_dir)
  month_file <- file.path(out_dir, sprintf("AllSpecies_%s_depth%d.RData",
                                           month_name, depth_val))
  if (file.exists(month_file)) {
    file.remove(month_file)
    message(sprintf("🗑️ Removed existing month summary for %s to rewrite", month_name))
  }
  save(month_df, file = month_file)
  message(sprintf("🎉 Completed month '%s' at depth %dcm", month_name, depth_val))
  return(month_df)
}

# ----------------------------------------------------------------------------
# Main loop: iterate over depths and months (with robust missing-data handling) 🌎
# -----------------------------------------------------------------------------
depths <- c(100, 150)

for (d in depths) {
  message(sprintf("🧭 ===== Depth = %d cm RUN START ===== 🧭", d))
  
  all_months <- mclapply(names(months_config), function(month_name) {
    message(sprintf("🚀 Processing month '%s' at depth %dcm...", month_name, d))
    process_month(month_name, d)
  }, mc.cores = n_cores)
  
  # Combine all months into one table
  combined <- bind_rows(all_months)
  
  # Standardize final dataframe
  final_df <- combined %>%
    select(year, month, species,
           soil_water_potential,
           transpiration_deficit,
           Quantiles, x, y, depth)
  
  # Save final object (always overwrite)
  raw_dir <- file.path("results", "Data")
  ensure_directory(raw_dir)
  final_file <- file.path(raw_dir, sprintf("final_df_depth%d.RData", d))
  if (file.exists(final_file)) {
    file.remove(final_file)
    message(sprintf("🗑️ Removed existing final_df_depth%d to rewrite", d))
  }
  save(final_df, file = final_file)
  message(sprintf("🥳 final_df_depth%d saved: %s", d, basename(final_file)))
  
  message(sprintf("🧭 ===== Depth = %d cm RUN COMPLETE ===== 🧭\n", d))
}

message("🏁 ALL DEPTHS PROCESSED SUCCESSFULLY! 🏁")
