setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)

# 1. NDVI export from netcdf
months_config <- list(
  July   = list(day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

out_dir <- "results_whole_Germany"
if (!dir.exists(out_dir)) dir.create(out_dir)

for (month in names(months_config)) {
  ndvi <- rast(months_config[[month]]$NDVI)
  years <- format(time(ndvi), "%Y")
  names(ndvi) <- paste0("NDVI_", years, "_", month)
  out_path <- file.path(out_dir, paste0("NDVI_", month, ".tif"))
  writeRaster(ndvi, filename = out_path, overwrite = TRUE)
}

# 2. Export PSI and TDiff for each species, month
species <- c("Beech", "Oak", "Spruce", "Pine")
months <- c("July", "August")
years <- 2003:2024
types <- c("psi", "tdiff")
output_names <- c(psi = "PSI", tdiff = "TDiff")  # All-lowercase names

for (sp in species) {
  for (mo in months) {
    for (ty in types) {
      input_dir <- file.path("results_monthly", mo, sp)
      tif_files <- file.path(input_dir, paste0(ty, "_", years, ".tif"))
      existing_idx <- file.exists(tif_files)
      tif_files_exist <- tif_files[existing_idx]
      if (length(tif_files_exist) > 0) {
        raster_stack <- rast(tif_files_exist)
        layer_years <- years[existing_idx]
        names(raster_stack) <- paste0(ty, "_", layer_years, "_", mo)
        out_prefix <- output_names[ty]
        output_file <- file.path("results_whole_Germany", paste0(sp, "_", out_prefix, "_", mo, ".tif"))
        writeRaster(raster_stack, filename = output_file, overwrite = TRUE)
        cat("Saved:", output_file, "\n")
      } else {
        cat("No files found for", sp, mo, ty, "\n")
      }
    }
  }
}

# 3. Species masks
Oak    <- rast("species_map_MODIS/Oak.tif")
Beech  <- rast("species_map_MODIS/Beech.tif")
Spruce <- rast("species_map_MODIS/Spruce.tif")
Pine   <- rast("species_map_MODIS/Pine.tif")
species_masks <- list(Oak=Oak, Beech=Beech, Spruce=Spruce, Pine=Pine)

# 4. NDVI and PSI filepaths for all species and months
get_raster_path <- function(sp, var, month) {
  file.path("results_whole_Germany", paste0(sp, "_", var, "_", month, ".tif"))
}

# 5. All cross-species pairs (NDVI/PSI, NDVI_species != PSI_species)
species_pairs <- expand.grid(
  NDVI_species = species,
  PSI_species  = species,
  stringsAsFactors = FALSE
) %>% filter(NDVI_species != PSI_species)

months <- c("July", "August")
years  <- 2003:2024

results_list <- list()

for (i in seq_len(nrow(species_pairs))) {
  ndvi_sp <- species_pairs$NDVI_species[i]
  psi_sp  <- species_pairs$PSI_species[i]
  cat("Processing NDVI:", ndvi_sp, "PSI:", psi_sp, "\n")
  
  for (month in months) {
    # NDVI/PSI rasters for this month/species
    ndvi_r <- rast(get_raster_path(ndvi_sp, "NDVI", month))
    psi_r  <- rast(get_raster_path(psi_sp, "PSI", month))
    mask_r <- species_masks[[ndvi_sp]]
    
    ndvi_r <- mask(ndvi_r, mask_r)
    psi_r  <- mask(psi_r, mask_r)
    
    # project PSI if necessary
    if (!compareGeom(ndvi_r, psi_r, stopOnError = FALSE)) {
      psi_r <- project(psi_r, ndvi_r)
    }
    
    # For all years: stack matching layers (assume same number of layers per month)
    n_layers <- min(nlyr(ndvi_r), nlyr(psi_r))
    for (lyr in seq_len(n_layers)) {
      ndvi_lyr <- ndvi_r[[lyr]]
      psi_lyr  <- psi_r[[lyr]]
      year     <- years[lyr]
      
      # Layer names for tracking
      names(ndvi_lyr) <- "NDVI"
      names(psi_lyr)  <- "PSI"
      comb_r <- c(ndvi_lyr, psi_lyr)
      df <- as.data.frame(comb_r, xy=TRUE, na.rm=TRUE)
      df$NDVI_species <- ndvi_sp
      df$PSI_species  <- psi_sp
      df$Year         <- year
      df$Month        <- month
      results_list[[paste(ndvi_sp, psi_sp, month, year, sep = "_")]] <- df
    }
  }
}

# 6. Combine to big dataframe
big_df <- bind_rows(results_list)

# 7. Bin and summarize NDVI by PSI for each NDVI_species, PSI_species, Year, Month
NDVI_PSIbin_multi <- function(df, bin_width = 50) {
  df <- na.omit(df)
  # Bin PSI
  psi_min <- floor(min(df$PSI, na.rm = TRUE))
  psi_max <- ceiling(max(df$PSI, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  df <- df %>%
    mutate(PSI_bin = cut(PSI, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  species_totals <- df %>%
    group_by(NDVI_species, PSI_species, Year, Month) %>%
    summarise(total_pixels = n(), .groups = "drop")
  meanNDVI_PSIbin_species <- df %>%
    group_by(NDVI_species, PSI_species, Year, Month, PSI_bin) %>%
    summarise(
      avg_value = mean(NDVI, na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = c("NDVI_species", "PSI_species", "Year", "Month")) %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.001) %>%
    select(NDVI_species, PSI_species, Year, Month, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_PSIbin_species)
}

big_bin_df <- NDVI_PSIbin_multi(big_df)

# 8. Colorblind-friendly palette
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# 9. Plot: NDVI (y) vs. PSI bin median (x), color = NDVI_species, facet = NDVI_species vs PSI_species
p <- ggplot(big_bin_df, aes(x = bin_median, y = avg_value, color = NDVI_species)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = cb_palette) +
  facet_grid(NDVI_species ~ PSI_species, labeller = label_both) +
  labs(
    x = "Soil water potential (kPa)",
    y = "NDVI",
    color = "NDVI species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 8)
  )
print(p)

# Optionally, save
ggsave("results_whole_Germany/NDVI_vs_PSI_crossspecies.png", p, width = 12, height = 8, dpi = 300)

# 10. Save summary table if desired
write.csv(big_bin_df, "results_whole_Germany/NDVI_PSI_bin_summary_crossspecies.csv", row.names = FALSE)

cat("All done!\n")
