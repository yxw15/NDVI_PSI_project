setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
results_dir <- "results_rootzone/Figures_till2022/supplementary"
big_df <- read.csv(file.path(results_dir, "pixel_pairs_diff_NDVI_PSI_rootzone.csv"))

library(tidyverse)

##### Read saved data #####
bin_summary <- read.csv(file.path(results_dir, "bin_summary_diff_NDVI_PSI_rootzone.csv"))

#─── Plot: scatter faceted by NDVI_species and DOY ─────────────────────────
cb_palette <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")
plot_df <- bin_summary %>%
  mutate(
    NDVI_species = factor(NDVI_species, levels = species),
    PSI_species  = factor(PSI_species,  levels = species),
    pixel_pct    = pixel_percentage * 100
  )

plot_df_clean <- plot_df %>%
  filter(!is.na(avg_NDVI), !is.na(bin_median))

# ─── 7b. Fit curves + points, faceted by NDVI_species only ───────────────────

# AIC-based selection between linear and exponential; force exponential for Spruce|Spruce
# AIC-only model selection (no forcing)
fit_curves <- function(df, key) {
  vals <- df$bin_median[is.finite(df$bin_median)]
  if (length(vals) < 2) {
    return(tibble(bin_median = numeric(0), fit = numeric(0), model = character(0)))
  }
  
  # Get NDVI/PSI species from group keys (correct!)
  spp_NDVI <- key$NDVI_species
  spp_PSI  <- key$PSI_species
  
  cat("\n==============================\n")
  print(paste("FITTING for:", spp_NDVI, "|", spp_PSI))
  
  # Linear
  lm_fit <- lm(avg_NDVI ~ bin_median, data = df)
  aic_lm <- AIC(lm_fit)
  print(paste("Linear AIC:", aic_lm))
  
  # Exponential
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  nls_fit <- tryCatch(
    nls(avg_NDVI ~ a + b * exp(c * bin_median),
        data = df, start = start_list, control = ctrl),
    error = function(e) NULL
  )
  
  if (!is.null(nls_fit)) {
    aic_exp <- AIC(nls_fit)
    print(paste("Exponential AIC:", aic_exp))
    print(summary(nls_fit)$coefficients)
  } else {
    aic_exp <- Inf
    print("Exponential model FAILED to converge.")
  }
  
  cat("Model chosen:", ifelse(aic_exp < aic_lm, "exponential", "linear"), "\n")
  
  grid <- tibble(bin_median = seq(min(vals), max(vals), length.out = 100))
  if (!is.null(nls_fit) && is.finite(aic_exp) && (aic_exp < aic_lm)) {
    grid$fit <- predict(nls_fit, newdata = grid)
    grid$model <- "exponential"
  } else {
    grid$fit <- predict(lm_fit, newdata = grid)
    grid$model <- "linear"
  }
  
  grid
}

# Compute fitted curves (grouped by NDVI & PSI species)
fitted_df <- plot_df_clean %>%
  group_by(NDVI_species, PSI_species) %>%
  group_modify(~ fit_curves(.x, .y)) %>%
  ungroup()

# Save fitted-data for reference
write.csv(
  fitted_df,
  file.path(results_dir, "fitted_data_diff_NDVI_PSI_rootzone.csv"),
  row.names = FALSE
)

#─── Plot: fitted curves + points, faceted by NDVI_species only ───────────
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

print(p_combined_fitted)
ggsave(
  file.path(results_dir, "NDVI_PSI_fitted_combined_rootzone_till2022.png"),
  plot  = p_combined_fitted,
  width = 10,
  height= 8,
  dpi   = 300
)
