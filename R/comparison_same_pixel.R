# ======================= PSI & TDiff cross-species analysis ===================
# Process -> GeoTIFFs/CSVs -> Maps -> Same-pixel correlations (PSI & TDIFF)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(terra)
  library(purrr)
  library(GGally)
  library(readr)
  library(patchwork)
})

# ---------------------------- Root & Config -----------------------------------
root_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project"
setwd(root_dir)

# will be set in the loop:
sample_cap  <- 10000L
target_year <- NA_integer_
target_date <- as.Date("2003-07-28")  # placeholder updated in loop

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

dir.create("results_rootzone/Project_mean", recursive = TRUE, showWarnings = FALSE)
dir.create("results_rootzone/Figures",      recursive = TRUE, showWarnings = FALSE)

# Template & borders
NDVI_template <- rast("../WZMAllDOYs/Quantiles_241.nc")[[1]]
germany_border_gk <- ne_countries(country = "Germany", scale = "medium", returnclass = "sf") |>
  st_transform(crs(NDVI_template))

# ------------------------------ Visual Style ----------------------------------
theme_consistent <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle     = element_text(hjust = 0.5, size = 14),
      axis.title        = element_text(face = "bold", size = 14),
      axis.text.y       = element_text(color = "black", size = 12),
      axis.text.x       = element_text(color = "black", size = 12),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "top",
      legend.text       = element_text(size = 12),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = 12)
    )
}

palette_for <- function(var) {
  if (tolower(var) == "psi") {
    rev(c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"))
  } else {
    c("blue", "dodgerblue", "cyan", "yellow", "orange", "red")
  }
}

save_plot <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi)
  message("✅ Saved: ", normalizePath(path, winslash = "/"))
}

# ------------------------------- File Helpers ---------------------------------
# Per-date files (rasters kept per-date to avoid collisions)
points_csv_path <- function(var) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_points_%s.csv", toupper(var), format(target_date, "%Y_%m_%d")))
}

projected_tif_path <- function(var, sp) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_%s_%s.tif", toupper(var), sp, format(target_date, "%Y_%m_%d")))
}

# NEW: per-year CSV (both dates for the year appended here)
points_csv_year_path <- function(var) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_points_%d.csv", toupper(var), target_year))
}

# Figure directory for combined Jul/Aug figures
figure_dir_julaug <- function() {
  file.path("results_rootzone/Figures",
            sprintf("%04d", target_year),
            "JulAug")
}

# --------------------------- NetCDF Date Extractor ----------------------------
extract_date_layer <- function(nc_path, date) {
  r_all <- rast(nc_path)
  idx   <- which(time(r_all) == date)
  if (length(idx) != 1) return(NULL)
  r_all[[idx]]
}

# ------------------------- 1) Process PSI & TDIFF -----------------------------
# Writes per-date rasters, and appends both dates into one per-year CSV per var.
process_variable_data <- function(var) {
  message(sprintf("--- Processing %s for %s ---", toupper(var), format(target_date)))
  dfs <- imap(species_config, function(paths, sp) {
    nc <- paths[[tolower(var)]]
    if (is.null(nc) || !file.exists(nc)) return(NULL)
    lyr <- extract_date_layer(nc, target_date)
    if (is.null(lyr)) return(NULL)
    proj <- project(lyr, NDVI_template)
    writeRaster(proj, projected_tif_path(var, sp), overwrite = TRUE)
    df <- as.data.frame(proj, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- var
    df$species  <- sp
    df$year     <- target_year
    df$monthday <- format(target_date, "%m-%d")  # keep which date
    df
  })
  dfs <- compact(dfs)
  if (!length(dfs)) return(invisible(FALSE))
  all_df <- bind_rows(dfs)
  
  # OPTIONAL: also keep per-date CSVs if you want them for debugging
  # write.csv(all_df, points_csv_path(var), row.names = FALSE)
  
  # Append into the per-year CSV (holds both 07-28 and 08-29)
  csv_year <- points_csv_year_path(var)
  if (file.exists(csv_year)) {
    old <- read.csv(csv_year)
    # Avoid duplicate appends if re-running; de-dup by keys if present
    all_df <- bind_rows(old, all_df) |>
      distinct(x, y, species, year, monthday, .keep_all = TRUE)
  }
  write.csv(all_df, csv_year, row.names = FALSE)
  
  invisible(TRUE)
}

# -------------------- 2) Combined 1x4 figure (overlay Jul & Aug) --------------
# Columns = species; inside each panel, ● Jul 28 and ▲ Aug 29 are overplotted.
plot_distribution_jul_aug_1x4 <- function(var, legend_title) {
  csv_year <- points_csv_year_path(var)
  if (!file.exists(csv_year)) {
    message("No per-year CSV found for ", toupper(var), " in ", target_year)
    return(invisible(FALSE))
  }
  all_df <- read.csv(csv_year)
  
  md_levels <- c("07-28","08-29")
  all_df <- all_df |>
    dplyr::filter(monthday %in% md_levels) |>
    dplyr::distinct(x, y, species, year, monthday, .keep_all = TRUE) |>
    dplyr::mutate(
      species  = factor(species, levels = c("Oak","Beech","Spruce","Pine")),
      monthday = factor(monthday, levels = md_levels, labels = c("Jul 28", "Aug 29"))
    )
  
  # use the actual column ("psi" or "tdiff") for color
  p <- ggplot(all_df, aes(x = x, y = y)) +
    geom_point(aes(color = .data[[var]], shape = monthday), size = 0.45, alpha = 0.7) +
    geom_sf(data = germany_border_gk, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(colours = palette_for(var), name = legend_title) +
    scale_shape_manual(values = c(16, 17), name = NULL) +  # ● Jul 28, ▲ Aug 29
    facet_wrap(~ species, nrow = 1) +
    coord_sf(crs = crs(NDVI_template), expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    theme_consistent(12) +
    guides(
      color = guide_colorbar(order = 1),
      shape = guide_legend(override.aes = list(size = 2), order = 0)
    )
  
  out_dir <- figure_dir_julaug()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir,
                   sprintf("%s_all_species_1x4_JulAug_%d.png", toupper(var), target_year))
  save_plot(p, out, width = 12, height = 4.5, dpi = 300)
}


# --------- (Optional) Correlations — keep available but not used in run -------
cross_species_correlations <- function(var, sample_cap = 10000L, seed = 1234) {
  date_tag <- format(target_date, "%Y_%m_%d")
  r_dir <- "results_rootzone/Project_mean"
  r_oak    <- rast(file.path(r_dir, sprintf("%s_Oak_%s.tif",   toupper(var), date_tag)))
  r_beech  <- rast(file.path(r_dir, sprintf("%s_Beech_%s.tif", toupper(var), date_tag)))
  r_spruce <- rast(file.path(r_dir, sprintf("%s_Spruce_%s.tif",toupper(var), date_tag)))
  r_pine   <- rast(file.path(r_dir, sprintf("%s_Pine_%s.tif",  toupper(var), date_tag)))
  
  names(r_oak) <- "Oak"; names(r_beech) <- "Beech"
  names(r_spruce) <- "Spruce"; names(r_pine) <- "Pine"
  stk <- c(r_oak, r_beech, r_spruce, r_pine)
  
  set.seed(seed)
  complete_mask <- terra::app(stk, function(v) as.integer(all(!is.na(v))))
  available <- as.integer(terra::global(complete_mask, "sum", na.rm = TRUE)[1,1])
  n_samp <- min(sample_cap, available)
  pts <- terra::spatSample(complete_mask, size = n_samp, method = "random",
                           as.points = TRUE, na.rm = TRUE, values = FALSE)
  samp_df <- terra::extract(stk, pts, xy = TRUE) |> dplyr::select(-ID)
  
  num_df <- dplyr::select(samp_df, Oak, Beech, Spruce, Pine)
  corr_pearson  <- cor(num_df, use = "complete.obs", method = "pearson")
  
  p_pairs <- GGally::ggpairs(num_df,
                             progress = FALSE,
                             upper = list(continuous = GGally::wrap("cor", size = 3)),
                             lower = list(continuous = GGally::wrap("points", alpha = 0.4, size = 0.6)),
                             diag  = list(continuous = GGally::wrap("densityDiag"))) +
    theme_consistent(12)
  pairs_ggally <- file.path(figure_dir_julaug(),
                            sprintf("%s_pairs_sample%d_%d.png", toupper(var), n_samp, target_year))
  save_plot(p_pairs, pairs_ggally, width = 10, height = 10, dpi = 300)
  
  corr_long <- as.data.frame(as.table(corr_pearson))
  names(corr_long) <- c("Var1","Var2","corr")
  p_heat <- ggplot(corr_long, aes(Var1, Var2, fill = corr)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", corr)), size = 4) +
    scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Pearson r") +
    coord_fixed() +
    labs(title = paste(toupper(var), "correlations across species"),
         x = NULL, y = NULL) +
    theme_consistent(12)
  heat_path <- file.path(figure_dir_julaug(),
                         sprintf("%s_correlation_heatmap_sample%d_%d.png",
                                 toupper(var), n_samp, target_year))
  save_plot(p_heat, heat_path, width = 6.5, height = 5.2, dpi = 300)
}

# ------------------------------- Run Pipeline ---------------------------------
years <- 2003:2024
monthdays <- c("07-28", "08-29")

for (yr in years) {
  target_year <- yr
  
  message("\n==============================")
  message(" Year ", target_year, " — processing both dates")
  message("==============================\n")
  
  # 1) Process both dates (fills the per-year CSV with monthday = 07-28 & 08-29)
  for (md in monthdays) {
    target_date <- as.Date(sprintf("%04d-%s", yr, md))
    message("  • Date: ", target_date)
    process_variable_data("psi")
    process_variable_data("tdiff")
  }
  
  # 2) After BOTH dates are processed, make the combined 1x4 overlay figures
  plot_distribution_jul_aug_1x4("psi",   "soil water potential (kPa)")
  plot_distribution_jul_aug_1x4("tdiff", "transpiration deficit (mm)")
  
  message("\nFigures written to: ", normalizePath(figure_dir_julaug(), winslash = "/"))
  print(list.files(figure_dir_julaug(), full.names = TRUE))
}

message("\nAll processing finished ✅")
# ==============================================================================
