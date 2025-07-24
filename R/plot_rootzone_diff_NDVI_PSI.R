# ─── 0. Setup ────────────────────────────────────────────────────────────────
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

# Create output directory
results_dir <- "results_rootzone"
if (!dir.exists(results_dir)) dir.create(results_dir)

# Define species, DOYs, years
species <- c("Oak", "Beech", "Spruce", "Pine")
doys    <- c("07-28", "08-29")  # Month-day strings
years   <- 2003:2024

#─── 1. Load rasters ─────────────────────────────────────────────────────────
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

#─── 2. Helper: subset PSI by DOY ─────────────────────────────────────────────
subset_by_doy <- function(psi_stack, mmdd) {
  dates <- time(psi_stack)
  idx   <- which(format(dates, "%m-%d") == mmdd)
  psi_stack[[idx]]
}

#─── 3. Process one NDVI_species × PSI_species × DOY ──────────────────────────
process_pair <- function(ndvi_sp, psi_sp, mmdd) {
  psi_sub <- subset_by_doy(PSI.rasters[[psi_sp]], mmdd)
  ndvi_sub <- NDVI.rasters[[mmdd]]
  
  # Mask NDVI to species footprint
  mask_nd <- species_masks[[ndvi_sp]]
  if (!compareGeom(ndvi_sub, mask_nd, stopOnError=FALSE)) {
    mask_nd <- project(mask_nd, ndvi_sub, method="near")
  }
  ndvi_m <- mask(ndvi_sub, mask_nd)
  
  # Project & mask PSI onto same grid
  if (!compareGeom(ndvi_m, psi_sub, stopOnError=FALSE)) {
    psi_sub <- project(psi_sub, ndvi_m)
  }
  psi_m <- mask(psi_sub, mask_nd)
  
  # Extract pixel pairs per year
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

#─── 4. Loop across all NDVI×PSI species combinations and DOYs ────────────────
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

#─── 5. Bin PSI and summarize NDVI ────────────────────────────────────────────
bin_summary <- big_df %>%
  filter(!is.na(NDVI), !is.na(PSI)) %>%
  mutate(PSI_bin = cut(PSI,
                       breaks = seq(floor(min(PSI)), ceiling(max(PSI)), by = 50),
                       include.lowest = TRUE, right = FALSE)) %>%
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
      # double-escaped backslashes so gsub() sees \[ \] \( \)
      nums <- as.numeric(strsplit(
        gsub("\\[|\\]|\\(|\\)", "", lbl),
        ","
      )[[1]])
      mean(nums)
    })
  ) %>%
  filter(pixel_percentage >= 0.001) %>%
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
    PSI_species  = factor(PSI_species,  levels = species),
    pixel_pct    = pixel_percentage * 100
  )

plot_df_clean <- plot_df %>%
  filter(!is.na(avg_NDVI), !is.na(bin_median))

p1 <- ggplot(plot_df, aes(
  x = bin_median,
  y = avg_NDVI,
  color = PSI_species,
  shape = PSI_species,
  size  = pixel_pct
)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = cb_palette, name = "PSI species") +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = "PSI species") +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_grid(rows = vars(DOY), cols = vars(NDVI_species)) +
  labs(x = "soil water potential (kPa)", y = "NDVI quantiles") +
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
  file.path(results_dir, "NDVI_PSI_scatter_rootzone.png"),
  plot = p1,
  width = 12,
  height = 8,
  dpi = 300
)

#─── 7. Fit curves and plot with linear & exponential ────────────────────────
library(broom)
library(dplyr)

fit_curves <- function(df) {
  # 1) guard against empty or all-NA bin_median
  vals <- df$bin_median[is.finite(df$bin_median)]
  if (length(vals) < 2) {
    # not enough data to fit a curve: return empty
    return(tibble(
      bin_median = numeric(0),
      fit        = numeric(0),
      model      = character(0)
    ))
  }
  
  # 2) fit linear
  lm_fit <- lm(avg_NDVI ~ bin_median, data = df)
  aic_lm <- AIC(lm_fit)
  
  # 3) try exponential
  nls_fit <- tryCatch({
    nls(avg_NDVI ~ a + b * exp(c * bin_median), data = df,
        start = list(
          a = min(df$avg_NDVI, na.rm = TRUE),
          b = (max(df$avg_NDVI, na.rm = TRUE) - min(df$avg_NDVI, na.rm = TRUE))/2,
          c = 0.001
        ),
        control = nls.control(maxiter = 1000))
  }, error = function(e) NULL)
  
  use_exp <- FALSE
  if (!is.null(nls_fit)) {
    aic_exp <- AIC(nls_fit)
    p_val   <- broom::tidy(nls_fit) %>%
      filter(term == "c") %>%
      pull(p.value)
    if (aic_exp < aic_lm && p_val < 0.05) use_exp <- TRUE
  }
  
  # 4) build grid over finite vals
  grid <- tibble(bin_median = seq(min(vals), max(vals), length.out = 100))
  
  if (use_exp) {
    grid <- grid %>%
      mutate(
        fit   = predict(nls_fit, newdata = grid),
        model = "exponential"
      )
  } else {
    grid <- grid %>%
      mutate(
        fit   = predict(lm_fit, newdata = grid),
        model = "linear"
      )
  }
  
  grid
}

# Compute fitted curves
fitted_df <- plot_df %>%
  group_by(NDVI_species, PSI_species) %>%
  group_modify(~ fit_curves(.x)) %>%
  ungroup()

# Save fitted-data table
write.csv(
  fitted_df,
  file.path(results_dir, "fitted_data_rootzone.csv"),
  row.names = FALSE
)

p2 <- ggplot(plot_df, aes(x = bin_median, y = avg_NDVI,
                          color = PSI_species, shape = PSI_species, size = pixel_pct)) +
  geom_point(alpha = 0.85) +
  geom_line(data = fitted_df,
            aes(bin_median, fit, linetype = model,
                group = interaction(PSI_species, model)),
            linewidth = 1, inherit.aes = FALSE) +
  scale_color_manual(values = cb_palette, name = expression(Psi[~~soil])) +
  scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = expression(Psi[~~soil])) +
  scale_size_continuous(name = "pixel percentage (%)", range = c(1, 8)) +
  facet_grid(rows = vars(DOY), cols = vars(NDVI_species)) +
  labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
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

# Save fitted-curve plot
ggsave(
  file.path(results_dir, "NDVI_PSI_fitted_rootzone.png"),
  plot = p2,
  width = 12,
  height = 8,
  dpi = 300
)


#─── 6b. Plot: scatter faceted by NDVI_species only ───────────────────────────

# we’ll reuse plot_df_clean from your existing code (it already has all DOYs)
p_combined_scatter <- ggplot(plot_df_clean, aes(
  x     = bin_median,
  y     = avg_NDVI,
  color = PSI_species,
  shape = PSI_species,
  size  = pixel_pct
)) +
  geom_point(alpha = 0.85) +
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
  file.path(results_dir, "NDVI_PSI_scatter_combined_rootzone.png"),
  plot = p_combined_scatter,
  width = 10,
  height = 8,
  dpi = 300
)


#─── 7b. Plot: fitted curves + points, faceted by NDVI_species only ───────────
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
      color    = PSI_species,               # <— add this
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
    x     = "soil water potential (kPa)",
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
  file.path(results_dir, "NDVI_PSI_fitted_combined_rootzone.png"),
  plot  = p_combined_fitted,
  width = 10,
  height= 8,
  dpi   = 300
)
