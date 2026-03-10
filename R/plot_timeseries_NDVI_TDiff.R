setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)
library(patchwork)

data <- combined %>%
  filter(month %in% c("July", "August")) %>%
  filter(year < 2023) %>%
  filter(Quantiles > 0) %>%
  na.omit()

# ---- theme ----
base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 14),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title        = element_text(size = 16),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_blank(),
    strip.text        = element_text(size = 14)
  )

plot_SI_time_series_NDVI_TD <- function(
    df_all,
    output_path = NULL,
    years = 2003:2022,
    months = NULL,
    species_order = c("Oak", "Beech", "Spruce", "Pine"),
    cb_palette = c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442"),
    base_theme = theme_bw(),
    legend_position = "bottom",
    width = 12,
    height = 9,
    dpi = 300
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
  })
  
  # ---- aggregate to yearly means by species ----
  # replace transpiration_deficit with your actual column name if needed
  df_species <- df_all %>%
    filter(year %in% years) %>%
    { if (!is.null(months) && "month" %in% names(.)) dplyr::filter(., month %in% months) else . } %>%
    group_by(year, species) %>%
    summarise(
      avg_quantile = mean(Quantiles, na.rm = TRUE),
      avg_td       = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(species = factor(species, levels = species_order))
  
  x_scale <- scale_x_continuous(
    breaks = seq(min(years), max(years), by = 2),
    limits = c(min(years), max(years))
  )
  
  # ---- (a) TDiff time series ----
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_td, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(
      x = "year",
      y = "transpiration deficit (mm)",
      color = ""
    ) +
    scale_color_manual(values = cb_palette) +
    x_scale +
    base_theme +
    ggtitle("(a)") +
    theme(
      plot.title.position = "panel",
      plot.title = element_text(size = 18, hjust = 0.0, vjust = -9, margin = margin(b = 5))
    )
  
  # ---- (b) regression: NDVI ~ transpiration deficit ----
  p2 <- ggplot(df_species, aes(x = avg_td, y = avg_quantile, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "transpiration deficit (mm)",
      y = expression(atop("NDVI quantiles", "(rank)"))
    ) +
    scale_color_manual(values = cb_palette) +
    base_theme +
    guides(color = "none") +
    ggtitle("(b)") +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  final_plot <- (p1 / p2) +
    plot_layout(guides = "collect") &
    theme(legend.position = legend_position)
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = final_plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  return(final_plot)
}

p <- plot_SI_time_series_NDVI_TD(
  df_all = data,
  output_path = "results_rootzone/Figures_till2022/SI_PSI/SI_time_series_NDVI_TD_species_till2022.png",
  base_theme = base_theme,
  legend_position = "top"
)

print(p)