# -----------------------------------------------------------------------------
# Title: NDVI, PSI, and TDiff Full Multi-Depth Pipeline Script
# Description:
#   1) Processes NDVI quantiles, soil water potential (PSI), and temperature difference (TDiff)
#      for multiple tree species across defined months and soil depths.
#   2) Skips already existing outputs to resume partial runs.
#   3) Combines per-species per-month results into depth-level datasets and saves them.
# -----------------------------------------------------------------------------

# 🎉 Start Script
message("🎉 Starting full NDVI-PSI-TDiff multi-depth pipeline...")

# 📚 Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ncdf4)
message("📚 Libraries loaded: terra, ggplot2, dplyr, tidyr, reshape2, ncdf4")

# 📂 Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
message("📂 Working directory set to: ", getwd())

# 🔧 Global parameters
dates_start <- "2003-01-01"
years       <- 2003:2024
# depths      <- c(50, 100, 150)
depths <- c(50)
message("🔧 Global settings: years ", min(years), "-", max(years), "; depths: ", paste(depths, collapse=", "))

# 🗓️ Month configuration
months_config <- list(
  April  = list(day="04-23", NDVI="../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(day="05-25", NDVI="../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(day="06-26", NDVI="../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(day="07-28", NDVI="../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day="08-29", NDVI="../WZMAllDOYs/Quantiles_241.nc")
)
message("🗓️ Months configured: ", paste(names(months_config), collapse=", "))

# 🌳 Species configuration
species_config <- list(
  Beech  = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
                mask="results/Species_Maps/Beech_mask.tif"),
  Oak    = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
                mask="results/Species_Maps/Oak_mask.tif"),
  Spruce = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
                mask="results/Species_Maps/Spruce_mask.tif"),
  Pine   = list(psi="../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
                tdiff="../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
                mask="results/Species_Maps/Pine_mask.tif")
)
message("🌳 Species configured: ", paste(names(species_config), collapse=", "))

# 📁 Utility: ensure directory exists
ensure_directory <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive=TRUE)
    message("📁 Created directory: ", path)
  }
}

# ----------------------------------------------------------------------------
# Data-transfer & save functions
# ----------------------------------------------------------------------------
transfer_psi_to_df <- function(nc_file) {
  nc <- nc_open(nc_file)
  time       <- ncvar_get(nc, "time")
  depth_vals <- ncvar_get(nc, "depth")
  x          <- ncvar_get(nc, "x")
  y          <- ncvar_get(nc, "y")
  psi        <- ncvar_get(nc, "psi")
  dates      <- as.Date(dates_start) + time
  nc_close(nc)
  dimnames(psi) <- list(x=x, y=y, depth=depth_vals, time=as.character(dates))
  reshape2::melt(psi,
                 varnames=c("x","y","depth","time"),
                 value.name="soil_water_potential") %>% na.omit()
}

save_psi_raster <- function(df, month_day, depth_val, out_dir) {
  sub <- df %>%
    filter(format(as.Date(time), "%m-%d")==month_day,
           depth==depth_val)
  sub$time <- as.Date(sub$time)
  ensure_directory(out_dir)
  for (yr in years) {
    pts <- sub %>% filter(format(time, "%Y")==yr)
    if (nrow(pts)==0) next
    r <- rast(pts[,c("x","y","soil_water_potential")], type="xyz")
    crs(r) <- "epsg:31467"
    writeRaster(r,
                file.path(out_dir, sprintf("psi_%d.tif", yr)),
                overwrite=TRUE)
  }
}

transfer_tdiff_to_df <- function(nc_file) {
  nc <- nc_open(nc_file)
  time  <- ncvar_get(nc, "time")
  x     <- ncvar_get(nc, "x")
  y     <- ncvar_get(nc, "y")
  tdiff <- ncvar_get(nc, "tdiff")
  dates <- as.Date(dates_start) + time
  nc_close(nc)
  dimnames(tdiff) <- list(x=x, y=y, time=as.character(dates))
  reshape2::melt(tdiff,
                 varnames=c("x","y","time"),
                 value.name="transpiration_deficit") %>% na.omit()
}

save_tdiff_raster <- function(df, month_day, depth_val, out_dir) {
  sub <- df %>% filter(format(as.Date(time), "%m-%d")==month_day)
  sub$time <- as.Date(sub$time)
  ensure_directory(out_dir)
  for (yr in years) {
    pts <- sub %>% filter(format(time, "%Y")==yr)
    if (nrow(pts)==0) next
    r <- rast(pts[,c("x","y","transpiration_deficit")], type="xyz")
    crs(r) <- "epsg:31467"
    writeRaster(r,
                file.path(out_dir, sprintf("tdiff_%d.tif", yr)),
                overwrite=TRUE)
  }
}

# ----------------------------------------------------------------------------
# Processing functions
# ----------------------------------------------------------------------------
process_NDVI_PSI_TDiff_species_year <- function(NDVI_file, species, year, out_dir, mask_r) {
  ndvi_r   <- rast(NDVI_file)[[year - min(years) + 1]]
  psi_r    <- project(rast(file.path(out_dir, sprintf("psi_%d.tif", year))), mask_r)
  tdiff_r  <- project(rast(file.path(out_dir, sprintf("tdiff_%d.tif", year))), mask_r)
  ndvi_m   <- mask(ndvi_r, mask_r)
  writeRaster(ndvi_m,
              file.path(out_dir, sprintf("NDVI_%d_mask.tif", year)),
              overwrite=TRUE)
  writeRaster(psi_r,
              file.path(out_dir, sprintf("PSI_%d_mask.tif", year)),
              overwrite=TRUE)
  writeRaster(tdiff_r,
              file.path(out_dir, sprintf("TDiff_%d_mask.tif", year)),
              overwrite=TRUE)
  df <- as.data.frame(c(ndvi_m, psi_r, tdiff_r), xy=TRUE, na.rm=TRUE)
  df$year    <- year
  df$species <- species
  df$depth   <- depth_val
  df
}

process_species <- function(species, month_name, depth_val) {
  sp_cfg <- species_config[[species]]
  mo_cfg <- months_config[[month_name]]
  out_dir <- file.path(sprintf("results_monthly_%d", depth_val), month_name, species)
  ensure_directory(out_dir)
  rdata_fp <- file.path(out_dir,
                        sprintf("NDVI_PSI_TDiff_%s_depth%d.RData",
                                gsub("-","", mo_cfg$day),
                                depth_val))
  if (file.exists(rdata_fp) && length(list.files(out_dir, "\\.tif$")) == length(years)*3) {
    load(rdata_fp)
    return(species_df)
  }
  psi_df   <- transfer_psi_to_df(sp_cfg$psi)
  tdiff_df <- transfer_tdiff_to_df(sp_cfg$tdiff)
  save_psi_raster(psi_df, mo_cfg$day, depth_val, out_dir)
  save_tdiff_raster(tdiff_df, mo_cfg$day, depth_val, out_dir)
  mask_r <- rast(sp_cfg$mask)
  lst    <- lapply(years, process_NDVI_PSI_TDiff_species_year,
                   NDVI_file=mo_cfg$NDVI,
                   species=species,
                   out_dir=out_dir,
                   mask_r=mask_r)
  species_df <- bind_rows(lst)
  save(species_df, file=rdata_fp)
  species_df
}

process_month <- function(month_name, depth_val) {
  dfs <- lapply(names(species_config), process_species,
                month_name=month_name,
                depth_val=depth_val)
  month_df <- bind_rows(dfs) %>% mutate(month=month_name)
  month_out <- file.path("results","Data", sprintf("depth_%d", depth_val))
  ensure_directory(month_out)
  save(month_df, file=file.path(month_out,
                                sprintf("AllSpecies_%s_depth%d.RData",
                                        month_name,
                                        depth_val)))
  month_df
}

# ----------------------------------------------------------------------------
# Main: run per-depth processing and save
# ----------------------------------------------------------------------------
for (d in depths) {
  message(sprintf("🚀 Processing depth %dcm...", d))
  month_dfs <- lapply(names(months_config), process_month, depth_val=d)
  depth_df  <- bind_rows(month_dfs)
  depth_df$depth <- d
  ensure_directory("results/Data")
  save(depth_df,
       file=file.path("results","Data",
                      sprintf("all_species_months_depth%d.RData", d)))
  message(sprintf("💾 Saved depth %d combined data.", d))
}
message("🎉 Full pipeline complete. Depth-level datasets saved under results/Data.")
