# ======================= PSI & TDiff cross-species analysis ===================
# Process -> GeoTIFFs/CSVs -> Same-pixel correlations (PSI & TDIFF)
# Output -> ONLY the pairs (scatter-matrix) figure, for 2018 Jul 28 & Aug 29
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(purrr)
  library(GGally)
  library(rlang)
})

# ---------------------------- Root & Config -----------------------------------
root_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project"
setwd(root_dir)

# Controls
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

# Template used for projection
NDVI_template <- rast("../WZMAllDOYs/Quantiles_241.nc")[[1]]

# ------------------------------ Visual Style ----------------------------------
theme_consistent <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title        = element_text(hjust = 0.5, size = base_size + 4, face = "bold"),
      plot.subtitle     = element_text(hjust = 0.5, size = base_size + 2),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text.y       = element_text(color = "black", size = base_size - 2),
      axis.text.x       = element_text(color = "black", size = base_size - 2),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "top",
      legend.text       = element_text(size = base_size - 2),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = base_size - 2)
    )
}

save_plot <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi)
  message("✅ Saved: ", normalizePath(path, winslash = "/"))
}

# ------------------------------- File Helpers ---------------------------------
# Per-year CSV (both dates appended here)
points_csv_year_path <- function(var) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_points_%d.csv", toupper(var), target_year))
}

# Per-date projected GeoTIFF
projected_tif_path <- function(var, sp) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_%s_%s.tif", toupper(var), sp, format(target_date, "%Y_%m_%d")))
}

# Figure directory
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
  dfs <- purrr::imap(species_config, function(paths, sp) {
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
    df$monthday <- format(target_date, "%m-%d")
    df
  })
  dfs <- purrr::compact(dfs)
  if (!length(dfs)) return(invisible(FALSE))
  all_df <- dplyr::bind_rows(dfs)
  
  # Append into per-year CSV (holds both 07-28 and 08-29)
  csv_year <- points_csv_year_path(var)
  if (file.exists(csv_year)) {
    old <- read.csv(csv_year)
    all_df <- dplyr::bind_rows(old, all_df) |>
      dplyr::distinct(x, y, species, year, monthday, .keep_all = TRUE)
  }
  write.csv(all_df, csv_year, row.names = FALSE)
  invisible(TRUE)
}

# ------------------ 2) Correlations — pool Jul 28 + Aug 29 together ----------
# Only writes the PAIRS figure (no heatmap). Show a blue 1:1 line; let axes free.
cross_species_correlations <- function(var,
                                       sample_cap = 10000L,
                                       seed = 1234,
                                       axis_title_size = 14,      # axis title size
                                       corr_text_size = NULL) {   # corr text size ties to species label size
  monthdays <- c("07-28", "08-29")
  r_dir <- "results_rootzone/Project_mean"
  
  # If corr_text_size not provided, match the species/variable label size (strip.text in theme_consistent)
  if (is.null(corr_text_size)) corr_text_size <- axis_title_size - 2
  
  # custom lower panel: points + blue 1:1 line; DO NOT force equal limits or coord_equal
  panel_scatter_11 <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.4, size = 0.6) +
      geom_abline(slope = 1, intercept = 0,
                  linewidth = 0.8, linetype = "solid", color = "blue") +
      scale_x_continuous(n.breaks = 3, guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(n.breaks = 3) +
      theme_consistent(axis_title_size) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(angle = 0)
      )
  }
  
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
  
  # build the pairs plot: corr text size == species label size
  p_pairs <- GGally::ggpairs(
    num_df,
    progress = FALSE,
    upper = list(continuous = GGally::wrap("cor", size = corr_text_size)),
    lower = list(continuous = panel_scatter_11),
    diag  = list(continuous = GGally::wrap("densityDiag"))
  ) +
    theme_consistent(axis_title_size) +
    theme(
      strip.text    = element_text(face = "bold", size = axis_title_size - 2),
      axis.title    = element_text(size = axis_title_size),
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.spacing = unit(0.8, "lines")
    )
  
  out_dir <- figure_dir_julaug()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pairs_ggally <- file.path(
    out_dir,
    sprintf("%s_pairs_pooledJulAug_sample%d_%d.png",
            toupper(var), nrow(num_df), target_year)
  )
  save_plot(p_pairs, pairs_ggally, width = 12, height = 12, dpi = 300)
  
  invisible(TRUE)
}

# ------------------------------- Run Pipeline ---------------------------------
# years <- c(2003, 2005, 2018, 2022)    
years <- c(2018, 2005)
monthdays <- c("07-28", "08-29")   # Jul 28 & Aug 29

for (yr in years) {
  target_year <- yr

  message("\n==============================")
  message(" Year ", target_year, " — processing both dates")
  message("==============================\n")

  # 1) Process both dates (creates per-date TIFFs used by the pairs plot)
  for (md in monthdays) {
    target_date <- as.Date(sprintf("%04d-%s", yr, md))
    message("  • Date: ", target_date)
    process_variable_data("psi")
    process_variable_data("tdiff")
  }
}

for (yr in years) {
  # 2) Pairs plot only (no other figures)
  cross_species_correlations("psi",   sample_cap = sample_cap,
                             axis_title_size = 14, corr_text_size = 6)
  cross_species_correlations("tdiff", sample_cap = sample_cap,
                             axis_title_size = 14, corr_text_size = 6)
  
  message("\nPairs figure(s) written to: ",
          normalizePath(figure_dir_julaug(), winslash = "/"))
  print(list.files(figure_dir_julaug(), full.names = TRUE, pattern = "pairs_"))
}

message("\nAll processing finished ✅")
# ==============================================================================
