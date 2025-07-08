# 🎉 Start Script
message("🎉 Starting NDVI-PSI-TDiff multi-depth processing script...")

# 📚 Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ncdf4)

# 📂 Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 🔧 Global parameters
dates_start <- "2003-01-01"
years <- 2003:2024
depths <- c(50, 100, 150)

# 🗓️ Month config
months_config <- list(
  April  = list(day="04-23", NDVI="../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(day="05-25", NDVI="../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(day="06-26", NDVI="../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(day="07-28", NDVI="../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day="08-29", NDVI="../WZMAllDOYs/Quantiles_241.nc")
)

# 🌳 Species config
species_config <- list(
  Beech  = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
                mask="species_map_MODIS/Beech.tif"),
  Oak    = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
                mask="species_map_MODIS/Oak.tif"),
  Spruce = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
                mask="species_map_MODIS/Spruce.tif"),
  Pine   = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
                mask="species_map_MODIS/Pine.tif")
)

ensure_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive=TRUE)
}

# -----------------------------------------------------------------------------
# Function: transfer_psi_to_df
#   Read PSI NetCDF into long data.frame (x, y, depth, time, soil_water_potential)
# -----------------------------------------------------------------------------
transfer_psi_to_df <- function(nc_file) {
  message("🌊 Loading PSI NetCDF: ", nc_file)
  nc <- nc_open(nc_file)
  time <- ncvar_get(nc, "time")
  depth_vals <- ncvar_get(nc, "depth")
  x <- ncvar_get(nc, "x"); y <- ncvar_get(nc, "y")
  psi <- ncvar_get(nc, "psi")
  dates <- as.Date(dates_start) + time
  nc_close(nc)
  
  dimnames(psi) <- list(x=x, y=y, depth=depth_vals, time=as.character(dates))
  df <- reshape2::melt(psi,
                       varnames=c("x","y","depth","time"),
                       value.name="soil_water_potential") %>% na.omit()
  message("✅ PSI df ready: ", nrow(df), " rows.")
  return(df)
}

# -----------------------------------------------------------------------------
# Function: save_psi_raster
#   Save yearly PSI rasters filtered by month_day and depth
# -----------------------------------------------------------------------------
save_psi_raster <- function(df, month_day, depth_val, out_dir) {
  message("🌊 Saving PSI rasters for depth ", depth_val, " on ", month_day)
  ensure_directory(out_dir)
  sub <- df %>%
    filter(format(as.Date(time), "%m-%d")==month_day,
           depth==depth_val)
  sub$time <- as.Date(sub$time)
  
  for (yr in years) {
    yr_pts <- sub %>% filter(format(time, "%Y")==yr)
    if (nrow(yr_pts)==0) next
    r <- rast(yr_pts[,c("x","y","soil_water_potential")], type="xyz")
    crs(r) <- "epsg:31467"
    path <- file.path(out_dir, sprintf("psi_%d.tif", yr))
    writeRaster(r, path, overwrite=TRUE)
    message("   🔹 PSI year ", yr, " saved.")
  }
  message("✅ PSI rasters done for depth ", depth_val)
}

# -----------------------------------------------------------------------------
# Function: transfer_tdiff_to_df
#   Read TDiff NetCDF into long data.frame (x, y, time, transpiration_deficit)
# -----------------------------------------------------------------------------
transfer_tdiff_to_df <- function(nc_file) {
  message("🔥 Loading TDiff NetCDF: ", nc_file)
  nc <- nc_open(nc_file)
  time <- ncvar_get(nc, "time")
  x <- ncvar_get(nc, "x"); y <- ncvar_get(nc, "y")
  tdiff <- ncvar_get(nc, "tdiff")
  dates <- as.Date(dates_start) + time
  nc_close(nc)
  
  dimnames(tdiff) <- list(x=x, y=y, time=as.character(dates))
  df <- reshape2::melt(tdiff,
                       varnames=c("x","y","time"),
                       value.name="transpiration_deficit") %>% na.omit()
  message("✅ TDiff df ready: ", nrow(df), " rows.")
  return(df)
}

# -----------------------------------------------------------------------------
# Function: save_tdiff_raster
#   Save yearly TDiff rasters filtered by month_day (no depth dimension)
# -----------------------------------------------------------------------------
save_tdiff_raster <- function(df, month_day, depth_val, out_dir) {
  message("🔥 Saving TDiff rasters for depth ", depth_val, " on ", month_day)
  ensure_directory(out_dir)
  sub <- df %>%
    filter(format(as.Date(time), "%m-%d")==month_day)
  sub$time <- as.Date(sub$time)
  
  for (yr in years) {
    yr_pts <- sub %>% filter(format(time, "%Y")==yr)
    if (nrow(yr_pts)==0) next
    r <- rast(yr_pts[,c("x","y","transpiration_deficit")], type="xyz")
    crs(r) <- "epsg:31467"
    path <- file.path(out_dir, sprintf("tdiff_%d.tif", yr))
    writeRaster(r, path, overwrite=TRUE)
    message("   🔹 TDiff year ", yr, " saved.")
  }
  message("✅ TDiff rasters done for depth ", depth_val)
}

# -----------------------------------------------------------------------------
# Function: process_NDVI_PSI_TDiff_species_year
# -----------------------------------------------------------------------------
process_NDVI_PSI_TDiff_species_year <- function(NDVI_file, species, year, out_dir) {
  message("🌿 [", species, "] Combining year ", year)
  ndvi <- rast(NDVI_file)[[year - min(years) + 1]]
  mask_r <- rast(species_config[[species]]$mask)
  
  psi_r   <- project(rast(file.path(out_dir, sprintf("psi_%d.tif", year))), mask_r)
  tdiff_r <- project(rast(file.path(out_dir, sprintf("tdiff_%d.tif", year))), mask_r)
  ndvi_m  <- mask(ndvi, mask_r)
  
  writeRaster(ndvi_m,  file.path(out_dir, sprintf("NDVI_%d_mask.tif", year)), overwrite=TRUE)
  writeRaster(psi_r,   file.path(out_dir, sprintf("PSI_%d_mask.tif", year)), overwrite=TRUE)
  writeRaster(tdiff_r, file.path(out_dir, sprintf("TDiff_%d_mask.tif", year)), overwrite=TRUE)
  
  df <- as.data.frame(c(ndvi_m, psi_r, tdiff_r), xy=TRUE, na.rm=TRUE)
  df$year    <- year
  df$species <- species
  names(df)  <- sub("^Quantiles_\\d+$", "Quantiles", names(df))
  message("   ✅ Combined data rows: ", nrow(df))
  return(df)
}

# -----------------------------------------------------------------------------
# Function: process_species
# -----------------------------------------------------------------------------
process_species <- function(species, month_name, depth_val) {
  message(sprintf("🎬 [%s | %dcm] -> %s", month_name, depth_val, species))
  sp_cfg <- species_config[[species]]
  mo_cfg <- months_config[[month_name]]
  
  out_dir <- file.path(sprintf("results_monthly_%d", depth_val), month_name, species)
  ensure_directory(out_dir)
  rdata_fp <- file.path(out_dir, sprintf("NDVI_PSI_TDiff_%s_depth%d.RData", gsub("-","", mo_cfg$day), depth_val))
  
  if (file.exists(rdata_fp)) {
    message("🎯 Already processed: ", species, " @ ", depth_val, " - skipping all processing.")
    load(rdata_fp)
    return(species_df)
  }
  
  psi_df   <- transfer_psi_to_df(sp_cfg$psi)
  tdiff_df <- transfer_tdiff_to_df(sp_cfg$tdiff)
  save_psi_raster(psi_df, mo_cfg$day, depth_val, out_dir)
  save_tdiff_raster(tdiff_df, mo_cfg$day, depth_val, out_dir)
  
  lst <- lapply(years, process_NDVI_PSI_TDiff_species_year, NDVI_file=mo_cfg$NDVI, species=species, out_dir=out_dir)
  species_df <- bind_rows(lst) %>% mutate(depth=depth_val)
  save(species_df, file=rdata_fp)
  message("📦 Saved species RData: ", rdata_fp)
  return(species_df)
}

# -----------------------------------------------------------------------------
# Function: process_month
# -----------------------------------------------------------------------------
process_month <- function(month_name, depth_val) {
  message(sprintf("🗓️ Processing month %s at depth %dcm", month_name, depth_val))
  ensure_directory(file.path("results","Data",sprintf("depth_%d", depth_val)))
  dfs <- lapply(names(species_config), process_species,
                month_name=month_name,
                depth_val=depth_val)
  month_df <- bind_rows(dfs) %>% mutate(month=month_name)
  
  out_dir <- file.path("results","Data",sprintf("depth_%d", depth_val))
  save(month_df, file=file.path(out_dir,
                                sprintf("AllSpecies_%s_depth%d.RData",month_name,depth_val)))
  message("💾 Saved month data: ", month_name, " @ depth ", depth_val)
  return(month_df)
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
for (d in depths) {
  message(sprintf("🚀 Starting full run for depth %dcm...", d))
  all_months <- lapply(names(months_config), process_month, depth_val=d)
  combined <- bind_rows(all_months) %>% mutate(depth=d)
  base_dir <- file.path("results","Data")
  ensure_directory(base_dir)
  save(combined, file=file.path(base_dir,
                                sprintf("AllSpecies_AllMonths_depth%d.RData", d)))
  message("📊 Combined all months for depth ", d)
}
message("🎉 All depths processed successfully! 🎉🎉🎉")

