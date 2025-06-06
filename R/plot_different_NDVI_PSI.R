setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

months_config <- list(
  July   = list(day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

# Ensure output directory exists
out_dir <- "results_whole_Germany"
if (!dir.exists(out_dir)) dir.create(out_dir)

for (month in names(months_config)) {
  ndvi <- rast(months_config[[month]]$NDVI)
  years <- format(time(ndvi), "%Y")
  # Set names as "YYYY_Month" (e.g. 2003_July)
  names(ndvi) <- paste0("NDVI_", years, "_", month)
  out_path <- file.path(out_dir, paste0("NDVI_", month, ".tif"))
  writeRaster(ndvi, filename = out_path, overwrite = TRUE)
}

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
        # Set names to "YYYY_Month"
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

Oak <- rast("species_map_MODIS/Oak.tif")
Beech <- rast("species_map_MODIS/Beech.tif")
Spruce <- rast("species_map_MODIS/Spruce.tif")
Pine <- rast("species_map_MODIS/Pine.tif")

NDVI.July <- rast("results_whole_Germany/NDVI_July.tif")
NDVI.July.Oak <- mask(NDVI.July, Oak)
NDVI.August <- rast("results_whole_Germany/NDVI_August.tif")
NDVI.August.Oak <- mask(NDVI.August, Oak)
NDVI.Oak <- c(NDVI.July.Oak, NDVI.August.Oak)
writeRaster(NDVI.Oak, filename = "results_whole_Germany/NDVI_Oak_July_August.tif", overwrite = TRUE)

Beech.PSI.July <- rast("results_whole_Germany/Beech_PSI_July.tif")
Beech.PSI.July.prjOak <- project(Beech.PSI.July, Oak)
Beech.PSI.July.prjOak <- mask(Beech.PSI.July.prjOak, Oak)

Beech.PSI.August <- rast("results_whole_Germany/Beech_PSI_August.tif")
Beech.PSI.August.prjOak <- project(Beech.PSI.August, Oak)
Beech.PSI.August.prjOak <- mask(Beech.PSI.August.prjOak, Oak)

Beech.PSI.prjOak <- c(Beech.PSI.July.prjOak, Beech.PSI.August.prjOak)
writeRaster(Beech.PSI.prjOak, filename = "results_whole_Germany/Beech_PSI_prjOak_July_August.tif", overwrite = TRUE)

Oak_NDVI.Beech_PSI <- c(NDVI.Oak, Beech.PSI.prjOak)

Oak_NDVI.Beech_PSI.df <- as.data.frame(Oak_NDVI.Beech_PSI, xy=T, na.rm=T)
Oak_NDVI.Beech_PSI.df.df_long <- Oak_NDVI.Beech_PSI.df %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "var_time",
    values_to = "value"
  )
Oak_NDVI.Beech_PSI.df.df_long <- Oak_NDVI.Beech_PSI.df.df_long %>%
  mutate(
    Variable = str_extract(var_time, "NDVI|psi"),
    Year = str_extract(var_time, "\\d{4}"),
    Month = str_extract(var_time, "July|August")
  )
Oak_NDVI.Beech_PSI.df.df_wide <- Oak_NDVI.Beech_PSI.df.df_long %>%
  select(-var_time) %>%
  pivot_wider(names_from = Variable, values_from = value)
Oak_NDVI.Beech_PSI.df.df_wide$species = "Beech"

NDVI_PSIbin <- function(df, bin_width = 50) {
  
  library(dplyr)
  
  df <- na.omit(df)
  
  # Identify correct value column
  value_column <- "NDVI"
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks
  psi_min <- floor(min(df$psi, na.rm = TRUE))
  psi_max <- ceiling(max(df$psi, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  # Bin the soil water potential values
  df <- df %>%
    mutate(PSI_bin = cut(psi, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute statistics for each species and bin, including filtering out bins with count < 2000
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.001) %>%
    select(species, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_PSIbin_species)
}

Oak_NDVI.Beech_PSIbin.df <- NDVI_PSIbin(Oak_NDVI.Beech_PSI.df.df_wide)

cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

ggplot(Oak_NDVI.Beech_PSIbin.df, aes(x = bin_median, y = avg_value, color = species)) +
  geom_point(size = 2) +
  scale_color_manual(values = cb_palette) +
  labs(
    x = "soil water potential (kPa)",
    y = "Oak NDVI",
    color = ""
  ) +
  theme_minimal()

