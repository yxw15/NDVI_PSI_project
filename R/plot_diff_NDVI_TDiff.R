# ─── 0. Setup ────────────────────────────────────────────────────────────────
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)

# Create output directory
results_dir <- "results_rootzone/Data"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Define species, DOYs, years
species <- c("Oak", "Beech", "Spruce", "Pine")
doys    <- c("07-28", "08-29")  # Month-day strings
years   <- 2003:2024

#─── 1. Load rasters ─────────────────────────────────────────────────────────
TDiff.rasters <- list(
  Oak    = rast("../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc"),
  Beech  = rast("../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc"),
  Spruce = rast("../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc"),
  Pine   = rast("../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc")
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

#─── 2. Helper: subset TDIFF by DOY ─────────────────────────────────────────────
subset_by_doy <- function(tdiff_stack, mmdd) {
  dates <- time(tdiff_stack)
  idx   <- which(format(dates, "%m-%d") == mmdd)
  tdiff_stack[[idx]]
}

#─── 3. Process one NDVI_species × tdiff_species × DOY ──────────────────────────
process_pair <- function(ndvi_sp, tdiff_sp, mmdd) {
  tdiff_sub <- subset_by_doy(TDiff.rasters[[tdiff_sp]], mmdd)
  ndvi_sub <- NDVI.rasters[[mmdd]]
  
  # Mask NDVI to species footprint
  mask_nd <- species_masks[[ndvi_sp]]
  if (!compareGeom(ndvi_sub, mask_nd, stopOnError=FALSE)) {
    mask_nd <- project(mask_nd, ndvi_sub, method="near")
  }
  ndvi_m <- mask(ndvi_sub, mask_nd)
  
  # Project & mask TDIFF onto same grid
  if (!compareGeom(ndvi_m, tdiff_sub, stopOnError=FALSE)) {
    tdiff_sub <- project(tdiff_sub, ndvi_m)
  }
  tdiff_m <- mask(tdiff_sub, mask_nd)
  
  # Extract pixel pairs per year
  yrs <- format(time(ndvi_m), "%Y")
  out <- vector("list", length(yrs))
  for (i in seq_along(yrs)) {
    combo <- c(ndvi_m[[i]], tdiff_m[[i]])
    df <- as.data.frame(combo, xy=TRUE, na.rm=TRUE)
    names(df)[3:4] <- c("NDVI","TDIFF")
    df$NDVI_species <- ndvi_sp
    df$tdiff_species  <- tdiff_sp
    df$DOY          <- mmdd
    df$Year         <- yrs[i]
    out[[i]] <- df
  }
  bind_rows(out)
}

#─── 4. Loop across all NDVI×TDIFF species combinations and DOYs ────────────────
pairs <- expand.grid(
  NDVI_sp = species,
  tdiff_sp  = species,
  DOY     = doys,
  stringsAsFactors = FALSE
)
res_list <- vector("list", nrow(pairs))
for (j in seq_len(nrow(pairs))) {
  nd_sp <- pairs$NDVI_sp[j]
  tdiff_sp <- pairs$tdiff_sp[j]
  md    <- pairs$DOY[j]
  message("Processing NDVI=", nd_sp, " | TDIFF=", tdiff_sp, " @ DOY=", md)
  res_list[[j]] <- process_pair(nd_sp, tdiff_sp, md)
}
big_df <- bind_rows(res_list)
write.csv(big_df,
          file.path(results_dir, "pixel_pairs_rootzone.csv"),
          row.names = FALSE)

#─── 5. Bin TDIFF and summarize NDVI ────────────────────────────────────────────
bin_summary <- big_df %>%
  filter(!is.na(NDVI), !is.na(TDIFF)) %>%
  mutate(TDIFF_bin = cut(TDIFF,
                         breaks = seq(floor(min(TDIFF)), ceiling(max(TDIFF)), by = 3),
                         include.lowest = TRUE, right = FALSE)) %>%
  group_by(NDVI_species, tdiff_species, DOY, TDIFF_bin) %>%
  summarise(
    avg_NDVI = mean(NDVI, na.rm = TRUE),
    count    = n(),
    .groups  = 'drop'
  ) %>%
  group_by(NDVI_species, tdiff_species, DOY) %>%
  mutate(
    total_pixels     = sum(count),
    pixel_percentage = count / total_pixels,
    bin_median       = sapply(as.character(TDIFF_bin), function(lbl) {
      nums <- as.numeric(strsplit(
        gsub("\\[|\\]|\\(|\\)", "", lbl),
        ","
      )[[1]])
      mean(nums)
    })
  ) %>%
  filter(pixel_percentage >= 0.0001) %>%
  ungroup()

# write out
write.csv(bin_summary,
          file.path(results_dir, "bin_summary_rootzone.csv"),
          row.names = FALSE)

#─── 6. Plot: scatter faceted by NDVI_species and DOY ─────────────────────────
cb_palette <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")
plot_df <- bin_summary %>%
  mutate(
    NDVI_species = factor(NDVI_species, levels = species),
    tdiff_species  = factor(tdiff_species,  levels = species),
    pixel_pct    = pixel_percentage * 100
  )

plot_df_clean <- plot_df %>%
  filter(!is.na(avg_NDVI), !is.na(bin_median))

p1 <- ggplot(plot_df, aes(
  x = bin_median,
  y = avg_NDVI,
  color = tdiff_species,
  shape = tdiff_species,
  size  = pixel_pct
)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = cb_palette, name = "TDIFF species") +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = "TDIFF species") +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_grid(rows = vars(DOY), cols = vars(NDVI_species)) +
  labs(x = "transpiration deficit (mm)", y = "NDVI quantiles") +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 14),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    shape = guide_legend(override.aes = list(size = 6)),
    size  = guide_legend()
  )

# Save scatter plot
ggsave(
  file.path(results_dir, "NDVI_TDIFF_scatter_rootzone.png"),
  plot = p1,
  width = 12,
  height = 8,
  dpi = 300
)

#─── 7. Robust model fitting (AICs, coefficients, exponential+linear) ─────────

fit_curves <- function(df) {
  vals <- df$bin_median[is.finite(df$bin_median)]
  if (length(vals) < 2) return(tibble())
  
  # Fit linear
  lm_fit <- tryCatch(lm(avg_NDVI ~ bin_median, data = df), error = function(e) NULL)
  aic_lm <- if (!is.null(lm_fit)) AIC(lm_fit) else NA
  
  # Fit exponential: NDVI = a + b * exp(-c * bin_median)
  start_nls <- list(a = 5, b = 7, c = 0.04)
  ctrl      <- nls.control(maxiter = 1200, minFactor = 1e-9)
  nls_fit <- tryCatch(
    nls(avg_NDVI ~ a + b * exp(-c * bin_median), data = df, start = start_nls, control = ctrl),
    error = function(e) NULL
  )
  aic_exp <- if (!is.null(nls_fit)) AIC(nls_fit) else NA
  
  # Get p-value for "c"
  p_val <- NA
  if (!is.null(nls_fit)) {
    s <- summary(nls_fit)$coefficients
    if ("c" %in% rownames(s)) p_val <- s["c", "Pr(>|t|)"]
  }
  
  # Use exponential if AIC lower and p-value significant
  use_exp <- !is.null(nls_fit) && !is.na(aic_exp) && (aic_exp < aic_lm) && !is.na(p_val) && (p_val < 0.05)
  
  grid <- tibble(bin_median = seq(min(vals), max(vals), length.out = 100))
  if (use_exp) {
    grid <- grid %>% mutate(fit = predict(nls_fit, newdata = grid), model = "exponential")
  } else if (!is.null(lm_fit)) {
    grid <- grid %>% mutate(fit = predict(lm_fit, newdata = grid), model = "linear")
  } else {
    return(tibble())
  }
  grid
}

# For AICs extraction
extract_aics <- function(df) {
  vals <- df$bin_median[is.finite(df$bin_median)]
  if (length(vals) < 2) return(tibble(AIC_linear=NA, AIC_exponential=NA))
  lm_fit <- tryCatch(lm(avg_NDVI ~ bin_median, data = df), error = function(e) NULL)
  aic_lm <- if (!is.null(lm_fit)) AIC(lm_fit) else NA
  start_nls <- list(a = 5, b = 7, c = 0.04)
  ctrl      <- nls.control(maxiter = 1200, minFactor = 1e-9)
  nls_fit   <- tryCatch(
    nls(avg_NDVI ~ a + b * exp(-c * bin_median), data = df, start = start_nls, control = ctrl),
    error = function(e) NULL
  )
  aic_exp <- if (!is.null(nls_fit)) AIC(nls_fit) else NA
  tibble(AIC_linear = aic_lm, AIC_exponential = aic_exp)
}

# For coefficient extraction (exponential only)
extract_exp_coefs <- function(df) {
  start_nls <- list(a = 5, b = 7, c = 0.04)
  ctrl      <- nls.control(maxiter = 1200, minFactor = 1e-9)
  nls_fit   <- tryCatch(
    nls(avg_NDVI ~ a + b * exp(-c * bin_median), data = df, start = start_nls, control = ctrl),
    error = function(e) NULL
  )
  if (is.null(nls_fit)) return(tibble())
  summ <- summary(nls_fit)$coefficients
  tibble(
    coef = rownames(summ),
    estimate = summ[, "Estimate"],
    pvalue = summ[, "Pr(>|t|)"]
  )
}

# Fit models
fitted_df <- plot_df %>%
  group_by(NDVI_species, tdiff_species) %>%
  group_modify(~ fit_curves(.x)) %>%
  ungroup()

aic_table <- plot_df %>%
  group_by(NDVI_species, tdiff_species) %>%
  group_modify(~ extract_aics(.x)) %>%
  ungroup()
write.csv(aic_table, file.path(results_dir, "aic_table_rootzone.csv"), row.names = FALSE)

coefs_table <- plot_df %>%
  group_by(NDVI_species, tdiff_species) %>%
  group_modify(~ extract_exp_coefs(.x)) %>%
  ungroup()
write.csv(coefs_table, file.path(results_dir, "exp_coeffs_rootzone.csv"), row.names = FALSE)

# Save fitted-data table
write.csv(
  fitted_df,
  file.path(results_dir, "fitted_data_rootzone.csv"),
  row.names = FALSE
)

#─── 8. Plots: fitted curves + points (faceted and combined) ────────────────

p2 <- ggplot(plot_df, aes(x = bin_median, y = avg_NDVI,
                          color = tdiff_species, shape = tdiff_species, size = pixel_pct)) +
  geom_point(alpha = 0.85) +
  geom_line(data = fitted_df,
            aes(bin_median, fit, linetype = model,
                group = interaction(tdiff_species, model)),
            linewidth = 1, inherit.aes = FALSE) +
  scale_color_manual(values = cb_palette, name = expression(TDIFF[~~soil])) +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = expression(TDIFF[~~soil])) +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_grid(rows = vars(DOY), cols = vars(NDVI_species)) +
  labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    size  = guide_legend()
  )

ggsave(
  file.path(results_dir, "NDVI_TDiff_fitted_rootzone.png"),
  plot = p2,
  width = 12,
  height = 8,
  dpi = 300
)

# Scatter, faceted by NDVI_species only
p_combined_scatter <- ggplot(plot_df_clean, aes(
  x     = bin_median,
  y     = avg_NDVI,
  color = tdiff_species,
  shape = tdiff_species,
  size  = pixel_pct
)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = cb_palette, name = expression(TDIFF[~~soil])) +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = expression(TDIFF[~~soil])) +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_wrap(~ NDVI_species, ncol = 2) +
  labs(
    x = "transpiration deficit (mm)",
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
  file.path(results_dir, "NDVI_TDiff_scatter_combined_rootzone.png"),
  plot = p_combined_scatter,
  width = 10,
  height = 8,
  dpi = 300
)

# Fitted curves + points, faceted by NDVI_species only
p_combined_fitted <- ggplot(plot_df_clean, aes(
  x     = bin_median,
  y     = avg_NDVI,
  color = tdiff_species,
  shape = tdiff_species,
  size  = pixel_pct
)) +
  geom_point(alpha = 0.85) +
  geom_line(
    data = fitted_df,
    aes(
      x        = bin_median,
      y        = fit,
      color    = tdiff_species,
      linetype = model,
      group    = interaction(tdiff_species, model)
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = cb_palette, name = expression(TDIFF[~~soil])) +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = expression(TDIFF[~~soil])) +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_wrap(~ NDVI_species, ncol = 2) +
  labs(
    x     = "transpiration deficit (mm)",
    y     = "NDVI quantiles (rank)",
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
  file.path(results_dir, "NDVI_TDiff_fitted_combined_rootzone.png"),
  plot  = p_combined_fitted,
  width = 10,
  height= 8,
  dpi   = 300
)
