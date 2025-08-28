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
points_csv_path <- function(var) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_points_%s.csv", toupper(var), format(target_date, "%Y_%m_%d")))
}
projected_tif_path <- function(var, sp) {
  file.path("results_rootzone/Project_mean",
            sprintf("%s_%s_%s.tif", toupper(var), sp, format(target_date, "%Y_%m_%d")))
}

# Figure directory is per-year and per-mm-dd; filenames end with YEAR only.
figure_dir <- function() {
  file.path("results_rootzone/Figures",
            sprintf("%04d", target_year),
            format(target_date, "%m-%d"))
}
figure_path <- function(filename) {
  file.path(figure_dir(), filename)
}

# --------------------------- NetCDF Date Extractor ----------------------------
extract_date_layer <- function(nc_path, date) {
  r_all <- rast(nc_path)
  idx   <- which(time(r_all) == date)
  if (length(idx) != 1) return(NULL)
  r_all[[idx]]
}

# ------------------------- 1) Process PSI & TDIFF -----------------------------
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
    df$species <- sp
    df$year    <- target_year
    df
  })
  dfs <- compact(dfs)
  if (!length(dfs)) return(invisible(FALSE))
  all_df <- bind_rows(dfs)
  write.csv(all_df, points_csv_path(var), row.names = FALSE)
  invisible(TRUE)
}

# -------------------- 2) Faceted ggplot maps (per variable) -------------------
plot_distribution <- function(var, legend_title) {
  csv_file <- points_csv_path(var)
  if (!file.exists(csv_file)) return(invisible(FALSE))
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
    theme_consistent(12)
  
  # filename ends with year only, and lives in <Figures>/<YEAR>/<MM-DD>/
  out <- figure_path(sprintf("%s_all_species_1x4_clean_%d.png", toupper(var), target_year))
  save_plot(p, out, width = 10, height = 4, dpi = 300)
}

# --------- 3) Same-pixel correlations (pairs, heatmap, selected pairs) -------
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
  
  # Correlations
  num_df <- dplyr::select(samp_df, Oak, Beech, Spruce, Pine)
  corr_pearson  <- cor(num_df, use = "complete.obs", method = "pearson")
  
  # A) GGally pairs
  p_pairs <- GGally::ggpairs(num_df,
                             progress = FALSE,
                             upper = list(continuous = GGally::wrap("cor", size = 3)),
                             lower = list(continuous = GGally::wrap("points", alpha = 0.4, size = 0.6)),
                             diag  = list(continuous = GGally::wrap("densityDiag"))) +
    theme_consistent(12)
  pairs_ggally <- figure_path(sprintf("%s_pairs_sample%d_%d.png", toupper(var), n_samp, target_year))
  save_plot(p_pairs, pairs_ggally, width = 10, height = 10, dpi = 300)
  
  # B) Heatmap
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
  heat_path <- figure_path(sprintf("%s_correlation_heatmap_sample%d_%d.png",
                                   toupper(var), n_samp, target_year))
  save_plot(p_heat, heat_path, width = 6.5, height = 5.2, dpi = 300)
  
  # C) Selected pairwise scatter panels
  plot_pair <- function(df, xvar, yvar) {
    r <- suppressWarnings(cor(df[[xvar]], df[[yvar]], use = "complete.obs"))
    ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
      geom_point(alpha = 0.4, size = 0.8) +
      suppressWarnings(geom_smooth(method = "lm", se = FALSE)) +
      annotate("text",
               x = min(df[[xvar]], na.rm = TRUE),
               y = max(df[[yvar]], na.rm = TRUE),
               label = paste0("r = ", round(r, 2)),
               hjust = 0, vjust = 1, size = 5, fontface = "bold") +
      labs(x = xvar, y = yvar) +
      theme_consistent(12)
  }
  pairs_to_plot <- list(
    c("Oak", "Beech"), c("Oak", "Spruce"), c("Oak", "Pine"),
    c("Beech", "Spruce"), c("Beech", "Pine"), c("Spruce", "Pine")
  )
  plots <- lapply(pairs_to_plot, function(p) plot_pair(num_df, p[1], p[2]))
  combined <- (plots[[1]] | plots[[2]] | plots[[3]]) /
    (plots[[4]] | plots[[5]] | plots[[6]])
  pairs_custom <- figure_path(sprintf("%s_pairs_corr_sample%d_%d.png", toupper(var), n_samp, target_year))
  save_plot(combined, pairs_custom, width = 12, height = 8, dpi = 300)
}

# ------------------------------- Run Pipeline ---------------------------------
run_for_date <- function(yr, month_day, sample_cap = 10000L) {
  target_year <<- yr
  target_date <<- as.Date(sprintf("%04d-%s", yr, month_day))
  
  message("\n==============================")
  message(" Processing ", target_date, " (Year: ", target_year, ")")
  message("==============================\n")
  
  # Ensure per-year/per-date figure dir exists
  dir.create(figure_dir(), recursive = TRUE, showWarnings = FALSE)
  
  # 1) Process variables
  process_variable_data("psi")
  process_variable_data("tdiff")
  
  # 2) Maps
  plot_distribution("psi",   "soil water potential (kPa)")
  plot_distribution("tdiff", "transpiration deficit (mm)")
  
  # 3) Cross-species correlations
  cross_species_correlations("psi",   sample_cap = sample_cap)
  cross_species_correlations("tdiff", sample_cap = sample_cap)
  
  message("\nContents of Figures dir for this run:")
  print(list.files(figure_dir(), full.names = TRUE))
}

# Years and dates to process
years <- 2003:2024
monthdays <- c("07-28", "08-29")

for (yr in years) {
  for (md in monthdays) {
    run_for_date(yr, md, sample_cap = sample_cap)
  }
}

message("\nAll processing finished ✅")
# ==============================================================================
