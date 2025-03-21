library(ncdf4)
library(terra)
library(reshape2)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Create output directory if it doesn't exist
if (!dir.exists("results/minPSI")) { 
  dir.create("results/minPSI", recursive = TRUE) 
}

transfer_psi_to_df <- function(nc_file, start_date) {
  nc <- nc_open(nc_file)
  time <- ncvar_get(nc, "time")
  depth <- ncvar_get(nc, "depth")
  x <- ncvar_get(nc, "x")
  y <- ncvar_get(nc, "y")
  psi <- ncvar_get(nc, "psi")
  dates <- as.Date(start_date) + time
  nc_close(nc)
  dimnames(psi) <- list(x = x, y = y, depth = depth, time = as.character(dates))
  psi_melted <- melt(psi, varnames = c("x", "y", "depth", "time"), value.name = "soil_water_potential")
  psi_melted <- na.omit(psi_melted)
  return(psi_melted)
}

species_list <- list(
  Beech  = list(psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
                mask   = "results/Species_Maps/Beech_mask.tif"),
  Oak    = list(psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
                mask   = "results/Species_Maps/Oak_mask.tif"),
  Spruce = list(psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
                mask   = "results/Species_Maps/Spruce_mask.tif"),
  Pine   = list(psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
                mask   = "results/Species_Maps/Pine_mask.tif")
)

start_date <- "2003-01-01"

#############################################
## PART 1: Unmasked Outputs
#############################################
for (species in names(species_list)) {
  cat("Processing unmasked outputs for", species, "\n")
  nc_file <- species_list[[species]]$psi_nc
  psi_df <- transfer_psi_to_df(nc_file, start_date)
  
  psi_agg <- psi_df %>%
    group_by(x, y, depth) %>%
    summarise(min_soil_water_potential = min(soil_water_potential, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(x = as.numeric(x),
           y = as.numeric(y),
           depth = as.numeric(depth))
  
  psi_agg_plot <- psi_agg %>%
    mutate(depth = factor(depth, levels = sort(unique(depth))))
  
  p <- ggplot(psi_agg_plot, aes(x = x, y = y, fill = min_soil_water_potential)) +
    geom_tile() +
    facet_wrap(~ depth, ncol = 3) +
    coord_equal() +
    scale_fill_viridis_c(option = "magma") +
    labs(title = paste(species, "Minimum Soil Water Potential (Unmasked)"),
         x = "X coordinate (m)",
         y = "Y coordinate (m)",
         fill = "Min PSI (kPa)") +
    theme_minimal()
  
  out_plot_file <- file.path("results", "minPSI", paste0(species, "_minPSI_threePanel_unmasked.png"))
  ggsave(out_plot_file, p, width = 12, height = 6)
  cat("Unmasked three-panel plot for", species, "saved to", out_plot_file, "\n")
  
  unique_depths <- sort(unique(psi_agg$depth))
  for (d in unique_depths) {
    psi_depth_df <- psi_agg %>% filter(depth == d) %>% select(x, y, min_soil_water_potential)
    r <- rast(psi_depth_df, type = "xyz")
    crs(r) <- "EPSG:31467"
    out_tif_file <- file.path("results", "minPSI", paste0(species, "_", d, "_minPSI_unmasked.tif"))
    writeRaster(r, out_tif_file, overwrite = TRUE)
    cat("Unmasked raster for", species, "depth", d, "saved to", out_tif_file, "\n")
  }
}

#############################################
## PART 2: Masked Outputs
#############################################
for (species in names(species_list)) {
  cat("Processing masked outputs for", species, "\n")
  nc_file <- species_list[[species]]$psi_nc
  psi_df <- transfer_psi_to_df(nc_file, start_date)
  
  psi_agg <- psi_df %>%
    group_by(x, y, depth) %>%
    summarise(min_soil_water_potential = min(soil_water_potential, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(x = as.numeric(x),
           y = as.numeric(y),
           depth = as.numeric(depth))
  
  mask_file <- species_list[[species]]$mask
  unique_depths <- sort(unique(psi_agg$depth))
  masked_list <- list()
  
  for (d in unique_depths) {
    cat("  Processing depth", d, "\n")
    
    psi_depth_df <- psi_agg %>% filter(depth == d) %>% select(x, y, min_soil_water_potential)
    
    r <- rast(psi_depth_df, type = "xyz")
    crs(r) <- "EPSG:31467"
    
    mask_r <- rast(mask_file)
    r_proj <- project(r, mask_r)
    r_masked <- mask(r_proj, mask_r)
    
    masked_tif_file <- file.path("results", "minPSI", paste0(species, "_", d, "_minPSI_masked.tif"))
    writeRaster(r_masked, masked_tif_file, overwrite = TRUE)
    cat("    Masked raster for", species, "depth", d, "saved to", masked_tif_file, "\n")
    
    masked_df <- as.data.frame(r_masked, xy = TRUE)
    masked_df <- na.omit(masked_df)
    colnames(masked_df)[3] <- "min_soil_water_potential"
    masked_df$depth <- d
    masked_list[[as.character(d)]] <- masked_df
  }
  
  masked_combined <- do.call(rbind, masked_list)
  masked_combined$depth <- factor(masked_combined$depth, levels = unique_depths)
  
  p_masked <- ggplot(masked_combined, aes(x = x, y = y, fill = min_soil_water_potential)) +
    geom_point(size = 0.5, shape = 21, stroke = 0) +
    facet_wrap(~ depth, ncol = 3) +
    coord_quickmap() +
    scale_fill_viridis_c(option = "magma") +
    labs(title = paste(species, "Minimum Soil Water Potential (Masked)"),
         x = "X",
         y = "Y",
         fill = "Min PSI (kPa)") +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(hjust = 0.5))
  
  masked_png_file <- file.path("results", "minPSI", paste0(species, "_minPSI_threePanel_masked.png"))
  ggsave(masked_png_file, p_masked, width = 12, height = 6)
  cat("Masked three-panel plot for", species, "saved to", masked_png_file, "\n")
}
