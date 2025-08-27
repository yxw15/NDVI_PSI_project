# --- Setup --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)              # unit()
  library(sf)                # spatial features
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(terra)
  library(purrr)
})

# Project root (adjust if needed)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# --- Config -------------------------------------------------------------------
target_date <- as.Date("2024-07-28")
target_year <- 2024

species_config <- list(
  Beech  = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc"
  ),
  Oak    = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc"
  ),
  Spruce = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc"
  ),
  Pine   = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc"
  )
)

# Color-blind–friendly species palette (used in legends/titles if needed)
cb_palette <- c(
  Oak    = "#E69F00",
  Beech  = "#0072B2",
  Spruce = "#009E73",
  Pine   = "#F0E442"
)

# Output dirs
dir.create("results_rootzone/Project_mean", recursive = TRUE, showWarnings = FALSE)
dir.create("results_rootzone/Figures",      recursive = TRUE, showWarnings = FALSE)

# Template & borders
NDVI_template <- rast("../WZMAllDOYs/Quantiles_241.nc")[[1]]
germany_border_gk <- ne_countries(country = "Germany", scale = "medium", returnclass = "sf") |>
  st_transform(crs(NDVI_template))

# --- Helpers ------------------------------------------------------------------
# Variable-specific color palette
palette_for <- function(var) {
  if (tolower(var) == "psi") {
    rev(c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"))
  } else {
    c("blue", "dodgerblue", "cyan", "yellow", "orange", "red")
  }
}

# Build file names
points_csv_path <- function(var) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_points_%s.csv", toupper(var), format(target_date, "%Y_%m_%d")))
}
projected_tif_path <- function(var, sp) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_%s_%s.tif", toupper(var), sp, format(target_date, "%Y_%m_%d")))
}

# Safe layer extractor for a given date
extract_date_layer <- function(nc_path, date) {
  r_all <- rast(nc_path)
  idx   <- which(time(r_all) == date)
  if (length(idx) != 1) return(NULL)
  r_all[[idx]]
}

# --- 1) Process data to GeoTIFFs & CSV ---------------------------------------
process_variable_data <- function(var) {
  message(sprintf("--- Processing %s ---", toupper(var)))
  
  dfs <- imap(species_config, function(paths, sp) {
    nc <- paths[[tolower(var)]]
    if (is.null(nc) || !file.exists(nc)) {
      warning(sprintf("[%s] Missing file for '%s': %s", sp, var, nc))
      return(NULL)
    }
    
    lyr <- extract_date_layer(nc, target_date)
    if (is.null(lyr)) {
      warning(sprintf("[%s] Date %s not found in %s", sp, target_date, basename(nc)))
      return(NULL)
    }
    
    # Reproject to template and save GeoTIFF
    proj <- project(lyr, NDVI_template)
    writeRaster(proj, projected_tif_path(var, sp), overwrite = TRUE)
    
    # Raster -> points
    df <- as.data.frame(proj, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- var
    df$species <- sp
    df$year    <- target_year
    df
  })
  
  dfs <- compact(dfs)
  if (length(dfs) == 0) {
    warning(sprintf("No data processed for %s", var)); return(invisible(FALSE))
  }
  
  all_df <- bind_rows(dfs)
  write.csv(all_df, points_csv_path(var), row.names = FALSE)
  message(sprintf("✓ Done: %s", points_csv_path(var)))
  invisible(TRUE)
}

# --- 2) Plot via ggplot2 (point cloud over Germany outline) -------------------
plot_distribution <- function(var, legend_title) {
  message(sprintf("--- ggplot2 map: %s ---", toupper(var)))
  
  csv_file <- points_csv_path(var)
  if (!file.exists(csv_file)) {
    warning(sprintf("CSV not found: %s", csv_file))
    return(invisible(FALSE))
  }
  
  all_df <- read.csv(csv_file)
  
  dfm <- all_df |>
    group_by(x, y, species, year) |>
    summarise(mean_val = mean(.data[[var]], na.rm = TRUE), .groups = "drop") |>
    mutate(species = factor(species, levels = c("Oak","Beech","Spruce","Pine")))
  
  pal <- palette_for(var)
  
  p <- ggplot(dfm, aes(x = x, y = y, color = mean_val)) +
    geom_point(size = 0.5) +
    geom_sf(data = germany_border_gk, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(colours = pal, name = legend_title) +
    facet_wrap(~ species, nrow = 1) +
    coord_sf(crs = crs(NDVI_template), expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle     = element_text(hjust = 0.5, size = 14),
      axis.title        = element_text(face = "bold", size = 14),
      axis.text.y       = element_text(color = "black", size = 12),
      axis.text.x       = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "top",
      legend.text       = element_text(size = 12),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = 12)
    )
  
  out <- file.path("results_rootzone/Figures",
                   sprintf("mean_%s_all_species_%s.png", var, format(target_date, "%Y_%m_%d")))
  ggsave(out, p, width = 10, height = 4, dpi = 300)
  message(sprintf("✓ Saved: %s", out))
  invisible(TRUE)
}

# --- 3) Plot via terra/base with shared legend (1x4 panels) -------------------
# Reads 4 GeoTIFFs, computes shared z-limits and draws a shared legend.
plot_four_panel <- function(var, legend_title = "") {
  message(sprintf("--- terra/base 1x4 panels: %s ---", toupper(var)))
  
  # Collect rasters per species in display order
  species_order <- c("Oak","Beech","Spruce","Pine")
  files <- setNames(vapply(species_order, \(sp) projected_tif_path(var, sp), character(1)), species_order)
  
  missing_files <- names(files)[!file.exists(files)]
  if (length(missing_files)) {
    warning(sprintf("Missing GeoTIFFs for %s: %s", var, paste(missing_files, collapse = ", ")))
    return(invisible(FALSE))
  }
  
  rasters <- lapply(files, rast)
  
  pal <- palette_for(var)
  all_vals <- do.call(c, lapply(rasters, values))
  all_vals <- all_vals[!is.na(all_vals)]
  zlim <- range(all_vals)
  
  legend_ticks  <- pretty(zlim, n = 6)
  legend_labels <- format(legend_ticks, trim = TRUE)
  
  out <- file.path("results_rootzone/Figures",
                   sprintf("%s_all_species_1x4_shared_legend.png", toupper(var)))
  png(out, width = 16, height = 5, units = "in", res = 300)
  
  # 4 panels + 1 legend column
  layout(matrix(c(1,2,3,4,5), nrow = 1), widths = c(1,1,1,1,0.35))
  par(mar = c(5, 5, 2, 1), oma = c(5, 5, 0, 0))
  
  plot_one <- function(r, title) {
    plot(r, main = "", col = pal, axes = TRUE, frame = FALSE,
         zlim = zlim, legend = FALSE, pax = list(cex.axis = 1.6))
    mtext(title, side = 3, line = 0.6, font = 2)
  }
  
  imap(rasters, ~plot_one(.x, .y))
  
  # Shared legend
  if (!requireNamespace("fields", quietly = TRUE)) {
    install.packages("fields", quiet = TRUE)
  }
  par(mar = c(5, 1.5, 2, 3))
  fields::image.plot(
    legend.only   = TRUE,
    zlim          = zlim,
    col           = pal,
    legend.width  = 2,
    legend.shrink = 0.8,
    axis.args     = list(at = legend_ticks, labels = legend_labels, cex.axis = 1.6),
    legend.args   = list(text = legend_title, side = 3, line = 1, font = 2, cex = 1.2)
  )
  
  mtext("longitude", side = 1, outer = TRUE, line = 0.8, font = 2)
  mtext("latitude",  side = 2, outer = TRUE, line = 0.8, font = 2)
  
  dev.off()
  message(sprintf("✓ Saved: %s", out))
  invisible(TRUE)
}

# --- Run pipeline -------------------------------------------------------------
# 1) Generate GeoTIFFs + CSV for both variables
process_variable_data("psi")
process_variable_data("tdiff")

# 2) ggplot2 method
plot_distribution("psi",   "soil water potential (kPa)")
plot_distribution("tdiff", "transpiration deficit (mm)")

# 3) terra/base method (shared legend 1x4 panels)
plot_four_panel("psi",   legend_title = "")
plot_four_panel("tdiff", legend_title = "")
