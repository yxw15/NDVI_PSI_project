# ─── 0. Setup ────────────────────────────────────────────────────────────────
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(tibble)
library(purrr)

# Create output directory
results_dir <- "results_rootzone/Diff_match/NDVI_PSI"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Define species, DOYs, years
species <- c("Oak", "Beech", "Spruce", "Pine")
doys    <- c("07-28", "08-29")  # Month-day strings
years   <- 2003:2024

# Color palette used later
cb_palette <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")

# ─── 1. Load rasters ─────────────────────────────────────────────────────────
PSI.rasters <- list(
  Oak    = rast("../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ei_bfv_20032024_compressed.nc"),
  Beech  = rast("../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Bu_bfv_20032024_compressed.nc"),
  Spruce = rast("../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Fi_bfv_20032024_compressed.nc"),
  Pine   = rast("../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ki_bfv_20032024_compressed.nc")
)

NDVI.rasters <- list(
  "07-28" = rast("../WZMAllDOYs/Quantiles_209.nc"),
  "08-29" = rast("../WZMAllDOYs/Quantiles_241.nc")
)

species_masks <- list(
  Oak    = rast("species_map_MODIS/Oak.tif"),
  Beech  = rast("species_map_MODIS/Beech.tif"),
  Spruce = rast("species_map_MODIS/Spruce.tif"),
  Pine   = rast("species_map_MODIS/Pine.tif")
)

# ─── 2. Helper: subset PSI by DOY ────────────────────────────────────────────
subset_by_doy <- function(psi_stack, mmdd) {
  dates <- time(psi_stack)
  idx   <- which(format(dates, "%m-%d") == mmdd)
  psi_stack[[idx]]
}

# ─── 3. Process one NDVI_species × PSI_species × DOY ─────────────────────────
process_pair <- function(ndvi_sp, psi_sp, mmdd) {
  psi_sub  <- subset_by_doy(PSI.rasters[[psi_sp]], mmdd)
  ndvi_sub <- NDVI.rasters[[mmdd]]
  
  # Mask NDVI to species footprint
  mask_nd <- species_masks[[ndvi_sp]]
  if (!compareGeom(ndvi_sub, mask_nd, stopOnError=FALSE)) {
    mask_nd <- project(mask_nd, ndvi_sub, method="near") # categorical mask → nearest
  }
  ndvi_m <- mask(ndvi_sub, mask_nd)
  
  # Project & mask PSI onto same grid
  if (!compareGeom(ndvi_m, psi_sub, stopOnError=FALSE)) {
    psi_sub <- project(psi_sub, ndvi_m)  # continuous → default (bilinear) OK
  }
  psi_m <- mask(psi_sub, mask_nd)
  
  # Extract pixel pairs per year (assumes NDVI & PSI layers are aligned by year)
  yrs <- format(time(ndvi_m), "%Y")
  out <- vector("list", length(yrs))
  for (i in seq_along(yrs)) {
    combo <- c(ndvi_m[[i]], psi_m[[i]])
    df <- as.data.frame(combo, xy=TRUE, na.rm=TRUE)
    names(df)[3:4] <- c("NDVI","PSI")
    df$NDVI_species <- ndvi_sp
    df$PSI_species  <- psi_sp
    df$DOY          <- mmdd
    df$Year         <- yrs[i]
    out[[i]] <- df
  }
  bind_rows(out)
}

# ─── 4. Build all combinations & extract matched pixels ──────────────────────
pairs <- expand.grid(
  NDVI_sp = species,
  PSI_sp  = species,
  DOY     = doys,
  stringsAsFactors = FALSE
)

res_list <- vector("list", nrow(pairs))
for (j in seq_len(nrow(pairs))) {
  nd_sp <- pairs$NDVI_sp[j]
  ps_sp <- pairs$PSI_sp[j]
  md    <- pairs$DOY[j]
  message("Processing NDVI=", nd_sp, " | PSI=", ps_sp, " @ DOY=", md)
  res_list[[j]] <- process_pair(nd_sp, ps_sp, md)
}
big_df <- bind_rows(res_list)

write.csv(big_df,
          file.path(results_dir, "pixel_pairs_rootzone.csv"),
          row.names = FALSE)

# ─── 5. Bin PSI and summarise NDVI ───────────────────────────────────────────
bin_summary <- big_df %>%
  filter(!is.na(NDVI), !is.na(PSI)) %>%
  mutate(
    PSI_bin = cut(
      PSI,
      breaks = seq(floor(min(PSI)), ceiling(max(PSI)), by = 50),
      include.lowest = TRUE, right = FALSE
    )
  ) %>%
  group_by(NDVI_species, PSI_species, DOY, PSI_bin) %>%
  summarise(
    avg_NDVI = mean(NDVI, na.rm = TRUE),
    count    = n(),
    .groups  = 'drop'
  ) %>%
  group_by(NDVI_species, PSI_species, DOY) %>%
  mutate(
    total_pixels     = sum(count),
    pixel_percentage = count / total_pixels,
    bin_median       = sapply(as.character(PSI_bin), function(lbl) {
      nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", lbl), ",")[[1]])
      mean(nums)
    })
  ) %>%
  filter(pixel_percentage >= 0.001) %>%   # keep bins with >= 0.1% of pixels
  ungroup()

write.csv(bin_summary,
          file.path(results_dir, "bin_summary_rootzone.csv"),
          row.names = FALSE)

# ─── 6. Prep data for plotting ───────────────────────────────────────────────
plot_df <- bin_summary %>%
  mutate(
    NDVI_species = factor(NDVI_species, levels = species),
    PSI_species  = factor(PSI_species,  levels = species),
    pixel_pct    = pixel_percentage * 100
  )

plot_df_clean <- plot_df %>%
  filter(!is.na(avg_NDVI), is.finite(bin_median))

# ─── 7b. Fit curves + points, faceted by NDVI_species only ───────────────────

# AIC-based selection between linear and exponential; force exponential for Spruce|Spruce
# AIC-only model selection (no forcing)
fit_curves <- function(df) {
  vals <- df$bin_median[is.finite(df$bin_median)]
  if (length(vals) < 2) {
    return(tibble(bin_median = numeric(0), fit = numeric(0), model = character(0)))
  }
  
  # Linear (always available & used as fallback)
  lm_fit <- lm(avg_NDVI ~ bin_median, data = df)
  aic_lm <- AIC(lm_fit)
  
  # Exponential (try)
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  nls_fit <- tryCatch(
    nls(avg_NDVI ~ a + b * exp(c * bin_median),
        data = df, start = start_list, control = ctrl),
    error = function(e) NULL
  )
  aic_exp <- if (!is.null(nls_fit)) AIC(nls_fit) else Inf
  
  # Choose the lower AIC (if nls converged)
  use_exp <- !is.null(nls_fit) && is.finite(aic_exp) && (aic_exp < aic_lm)
  
  grid <- tibble(bin_median = seq(min(vals), max(vals), length.out = 100))
  if (use_exp) {
    grid$fit   <- predict(nls_fit, newdata = grid)
    grid$model <- "exponential"
  } else {
    grid$fit   <- predict(lm_fit, newdata = grid)
    grid$model <- "linear"
  }
  grid
}

# Compute fitted curves (grouped by NDVI & PSI species)
fitted_df <- plot_df_clean %>%
  dplyr::group_by(NDVI_species, PSI_species) %>%
  dplyr::group_modify(~ fit_curves(.x)) %>%
  dplyr::ungroup()

# Save fitted-data for reference
write.csv(
  fitted_df,
  file.path(results_dir, "fitted_data_rootzone.csv"),
  row.names = FALSE
)

# 7b plot
p_combined_fitted <- ggplot(plot_df_clean, aes(
  x     = bin_median,
  y     = avg_NDVI,
  color = PSI_species,
  shape = PSI_species,
  size  = pixel_pct
)) +
  geom_point(alpha = 0.85) +
  geom_line(
    data = fitted_df,
    aes(
      x        = bin_median,
      y        = fit,
      color    = PSI_species,
      linetype = model,
      group    = interaction(PSI_species, model)
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = cb_palette, name = expression(Psi[~~soil])) +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = expression(Psi[~~soil])) +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_wrap(~ NDVI_species, ncol = 2) +
  labs(
    x = "soil water potential (kPa)",
    y = "NDVI quantiles (rank)",
    title = ""
  ) +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title        = element_text(face = "bold", size = 14),
    axis.text         = element_text(color = "black", size = 12),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12),
    legend.position   = "right",
    legend.text       = element_text(size = 12)
  )

ggsave(
  file.path(results_dir, "NDVI_PSI_fitted_combined_rootzone.png"),
  plot  = p_combined_fitted,
  width = 10,
  height= 8,
  dpi   = 300
)
