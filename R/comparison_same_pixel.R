# ======================= PSI & TDiff cross-species analysis ===================
# Process -> GeoTIFFs/CSVs -> 2x4 Maps -> Same-pixel correlations (PSI & TDIFF)
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
# (unchanged theme)
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

# Per-year CSV (both dates appended here)
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
  
  # OPTIONAL: per-date CSV
  # write.csv(all_df, points_csv_path(var), row.names = FALSE)
  
  # Append into per-year CSV (holds both 07-28 and 08-29)
  csv_year <- points_csv_year_path(var)
  if (file.exists(csv_year)) {
    old <- read.csv(csv_year)
    # Avoid duplicates if re-running
    all_df <- bind_rows(old, all_df) |>
      distinct(x, y, species, year, monthday, .keep_all = TRUE)
  }
  write.csv(all_df, csv_year, row.names = FALSE)
  
  invisible(TRUE)
}

# -------------------- 2) Combined 2x4 figure (rows=months, cols=species) ------
# Row 1 = Jul 28, Row 2 = Aug 29; no shape legend needed.
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
      monthday = factor(monthday, levels = md_levels, labels = c("DOY 209 (July)", "DOY 241 (August)"))
    )
  
  p <- ggplot(all_df, aes(x = x, y = y)) +
    geom_point(aes(color = .data[[var]]), size = 0.45, alpha = 0.7) +
    geom_sf(data = germany_border_gk, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(colours = palette_for(var), name = legend_title) +
    facet_grid(monthday ~ species) +
    coord_sf(crs = crs(NDVI_template), expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    theme_consistent(12) +
    guides(
      color = guide_colorbar(
        order = 1,
        barwidth = unit(6, "cm"),   # << make bar longer
        barheight = unit(0.5, "cm") # << keep it slim so text fits
      )
    )
  
  out_dir <- figure_dir_julaug()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(out_dir,
                   sprintf("%s_all_species_2x4_JulAug_%d.png", toupper(var), target_year))
  save_plot(p, out, width = 12, height = 8.5, dpi = 300)
}

# ------------------ 3) Correlations — pool Jul 28 + Aug 29 together ----------
# Pairs plot with: corr text same size as species name (12),
# blue bold-ish 1:1 line, and identical X/Y limits per panel.
cross_species_correlations <- function(var, sample_cap = 10000L, seed = 1234) {
  monthdays <- c("07-28", "08-29")
  r_dir <- "results_rootzone/Project_mean"
  
  sample_one_month <- function(md, n_samp, seed_off = 0) {
    date_tag <- format(as.Date(sprintf("%04d-%s", target_year, md)), "%Y_%m_%d")
    
    r_oak    <- rast(file.path(r_dir, sprintf("%s_Oak_%s.tif",    toupper(var), date_tag)))
    r_beech  <- rast(file.path(r_dir, sprintf("%s_Beech_%s.tif",  toupper(var), date_tag)))
    r_spruce <- rast(file.path(r_dir, sprintf("%s_Spruce_%s.tif", toupper(var), date_tag)))
    r_pine   <- rast(file.path(r_dir, sprintf("%s_Pine_%s.tif",   toupper(var), date_tag)))
    
    names(r_oak) <- "Oak"; names(r_beech) <- "Beech"
    names(r_spruce) <- "Spruce"; names(r_pine) <- "Pine"
    stk <- c(r_oak, r_beech, r_spruce, r_pine)
    
    complete_mask <- terra::app(stk, function(v) as.integer(all(!is.na(v))))
    available <- as.integer(terra::global(complete_mask, "sum", na.rm = TRUE)[1,1])
    if (is.na(available) || available == 0) return(NULL)
    
    n_take <- min(n_samp, available)
    set.seed(seed + seed_off)
    pts <- terra::spatSample(complete_mask, size = n_take, method = "random",
                             as.points = TRUE, na.rm = TRUE, values = FALSE)
    terra::extract(stk, pts, xy = TRUE) |>
      dplyr::select(-ID) |>
      dplyr::mutate(monthday = md)
  }
  
  n_per_month <- max(1L, floor(sample_cap / length(monthdays)))
  samp_list <- list(
    sample_one_month("07-28", n_per_month, seed_off = 0),
    sample_one_month("08-29", n_per_month, seed_off = 1)
  )
  samp_list <- purrr::compact(samp_list)
  if (!length(samp_list)) {
    message("No data available for correlations in year ", target_year)
    return(invisible(FALSE))
  }
  
  samp_df <- dplyr::bind_rows(samp_list)
  num_df  <- dplyr::select(samp_df, Oak, Beech, Spruce, Pine)
  
  # Common axis limits across all species so X/Y are identical in every panel
  axis_limits <- range(as.matrix(num_df), finite = TRUE, na.rm = TRUE)
  
  # Custom lower panel: points + bold(ish) blue 1:1 line + equal aspect + shared limits
  panel_scatter_11 <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.4, size = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 0.9, lineend = "round") +
      scale_x_continuous(limits = axis_limits) +
      scale_y_continuous(limits = axis_limits) +
      coord_equal() +
      theme_consistent(12)
  }
  
  # Correlation matrix (Pearson)
  corr_pearson <- cor(num_df, use = "complete.obs", method = "pearson")
  
  # Pairs plot: corr text sized like species strip labels (12)
  p_pairs <- GGally::ggpairs(
    num_df,
    progress = FALSE,
    upper = list(continuous = GGally::wrap("cor", size = 12)),
    lower = list(continuous = panel_scatter_11),
    diag  = list(continuous = GGally::wrap("densityDiag"))
  ) +
    theme_consistent(12)
  
  out_dir <- figure_dir_julaug()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pairs_ggally <- file.path(
    out_dir,
    sprintf("%s_pairs_pooledJulAug_sample%d_%d.png",
            toupper(var), nrow(num_df), target_year)
  )
  save_plot(p_pairs, pairs_ggally, width = 10, height = 10, dpi = 300)
  
  invisible(TRUE)
}

# ------------------------------- Run Pipeline ---------------------------------
years <- c(2005, 2022)
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
  
  # 2) After BOTH dates are processed, make the combined 2x4 figures
  plot_distribution_jul_aug_1x4("psi",   "soil water potential (kPa)")
  plot_distribution_jul_aug_1x4("tdiff", "transpiration deficit (mm)")
  
  # (Optional) enable if you want correlations each year:
  # cross_species_correlations("psi",   sample_cap = sample_cap)
  # cross_species_correlations("tdiff", sample_cap = sample_cap)
  
  message("\nFigures written to: ", normalizePath(figure_dir_julaug(), winslash = "/"))
  print(list.files(figure_dir_julaug(), full.names = TRUE))
}

message("\nAll processing finished ✅")
# ==============================================================================
