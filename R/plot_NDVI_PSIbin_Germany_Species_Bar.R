load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")

# ---- theme & axes ----
# 1. Visual Style only
base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# 2. Legend Ordering only
base_guides <- guides(
  size  = guide_legend(order = 1, override.aes = list(alpha = 1)),
  color = guide_legend(order = 2, override.aes = list(linewidth = 1.2, alpha = 1)),
  shape = guide_legend(order = 2)
)


# ---- NDVI_PSIbin Germany vs Species Functions  ----

NDVI_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  # totals per (area, species) for percentage
  species_totals <- df %>%
    group_by(area, species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # global bin breaks so bins consistent across everything
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df %>%
    mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    group_by(area, species, PSI_bin) %>%
    summarise(
      bin_mean  = mean(soil_water_potential, na.rm = TRUE),
      avg_value = mean(NDVI, na.rm = TRUE),
      count     = n(),
      .groups   = "drop"
    ) %>%
    left_join(species_totals, by = c("area", "species")) %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    select(area, species, PSI_bin, bin_mean, avg_value, count, total_pixels, percentage)
}

NDVI_diff_species_minus_Germany_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df_binned <- df %>%
    mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  totals <- df_binned %>%
    group_by(area, species) %>%
    summarise(total_n = n(), .groups = "drop")
  
  sp_tbl <- df_binned %>%
    filter(area == "species") %>%
    group_by(species, PSI_bin) %>%
    summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      sp_mean  = mean(NDVI, na.rm = TRUE),
      sp_n     = n(),
      .groups  = "drop"
    ) %>%
    filter(sp_n >= min_count) %>%
    left_join(
      totals %>% filter(area == "species") %>% select(species, sp_total = total_n),
      by = "species"
    ) %>%
    mutate(sp_pct = sp_n / sp_total)
  
  de_tbl <- df_binned %>%
    filter(area == "Germany") %>%
    group_by(species, PSI_bin) %>%
    summarise(
      de_mean = mean(NDVI, na.rm = TRUE),
      de_n    = n(),
      .groups = "drop"
    ) %>%
    filter(de_n >= min_count) %>%
    left_join(
      totals %>% filter(area == "Germany") %>% select(species, de_total = total_n),
      by = "species"
    ) %>%
    mutate(de_pct = de_n / de_total)
  
  sp_tbl %>%
    inner_join(de_tbl, by = c("species", "PSI_bin")) %>%
    mutate(
      diff_NDVI = sp_mean - de_mean,
      pct_bin   = pmin(sp_pct, de_pct)
    ) %>%
    select(species, PSI_bin, bin_mean, sp_mean, de_mean, sp_n, de_n, sp_pct, de_pct, pct_bin, diff_NDVI)
}

# ---- NDVI_PSIbin Germany vs Species (Bar difference) ----

plot_F4_NDVI_PSIbin_area_Bar_Diff <- function(
    df,
    out_file = "results_rootzone/main_PSI/F4_NDVI_PSIbin_area_Bar_Diff.png",
    bin_width = 50,
    min_count = 1000,
    k_left = 3,                 # fixed default
    gam_method = "REML",         # fixed default
    bs_left = "tp",              # NEW: avoid undefined object in smooth
    width = 12,
    height = 13,
    dpi = 300,
    sp_levels = c("Oak", "Beech", "Spruce", "Pine"),
    sp_cols   = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    sp_shapes = c(Oak=16, Beech=17, Spruce=15, Pine=18),
    germany_grey = "grey50"
) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2); library(mgcv); library(scales); library(patchwork)
  })
  
  d_all <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  d_sp <- d_all %>% filter(area == "species")
  d_de <- d_all %>% filter(area == "Germany")
  
  legend_cols   <- c(sp_cols, "Germany" = germany_grey)
  legend_fills  <- c(sp_cols, "Germany" = "grey70")
  legend_shapes <- c(sp_shapes, "Germany" = 32)
  
  p_left <- ggplot(d_sp, aes(bin_mean, avg_value)) +
    # --- Germany: Points + Smooth with CI (REML, k=3) ---
    geom_point(
      data = d_de,
      aes(size = percentage),
      colour = germany_grey,
      alpha = 0.2,
      show.legend = FALSE
    ) +
    geom_smooth(
      data = d_de,
      aes(x = bin_mean, y = avg_value, colour = "Germany", fill = "Germany"),
      method = "gam",
      formula = y ~ s(x, k = k_left),
      method.args = list(method = gam_method),
      se = TRUE,
      alpha = 0.2,
      linewidth = 1
    ) +
    
    # --- Species: Points + Smooth with CI (REML, k=3) ---
    geom_point(aes(colour = species, shape = species, size = percentage), alpha = 0.4) +
    geom_smooth(
      aes(colour = species, fill = species),
      method = "gam",
      formula = y ~ s(x, k = k_left, bs = bs_left),
      method.args = list(method = gam_method),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1.2
    ) +
    
    facet_wrap(~ species, nrow = 2) +
    scale_colour_manual(values = legend_cols, name = "") +
    scale_fill_manual(values = legend_fills, name = "") +
    scale_shape_manual(values = legend_shapes, name = "") +
    scale_size_continuous(range = c(1, 8), labels = percent_format(accuracy = 1), name = "") +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)", title = "(a)") +
    base_theme +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          shape = 16,
          colour = "grey80",
          alpha = 1,
          fill = NA,
          linewidth = 0
        )
      ),
      colour = guide_legend(
        order = 2,
        override.aes = list(
          linewidth = 1,
          alpha = 0.4,
          fill = legend_fills
        )
      ),
      shape = "none",
      fill = "none"
    ) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  # --- Bar Logic ---
  df_diff <- NDVI_diff_species_minus_Germany_PSIbin(df, bin_width, min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  df_bar <- df_diff %>%
    group_by(species) %>%
    summarise(
      mean_unweighted = mean(diff_NDVI, na.rm = TRUE),
      mean_weighted   = weighted.mean(diff_NDVI, w = pct_bin, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(mean_weighted, mean_unweighted),
      names_to = "type",
      values_to = "value"
    ) %>%
    mutate(type = recode(type,
                         mean_weighted = "Weighted (by min area support)",
                         mean_unweighted = "Unweighted"
    ))
  
  p_right_top <- ggplot(
    df_bar %>% filter(type == "Weighted (by min area support)"),
    aes(x = species, y = value, fill = species)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) +
    scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = expression(Delta[paste("  ", NDVI, " ", quantiles, " ", (rank))]), title = "(b) weighted mean") +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  p_right_bottom <- ggplot(
    df_bar %>% filter(type == "Unweighted"),
    aes(x = species, y = value, fill = species)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) +
    scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = expression(Delta[paste("  ", NDVI, " ", quantiles, " ", (rank))]), title = "(c) unweighted mean") +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  p_bottom <- (p_right_top | p_right_bottom) + plot_layout(widths = c(1, 1))
  
  p <- (p_left / p_bottom) +
    plot_layout(heights = c(2.2, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  return(p)
}

plot_F4_NDVI_PSIbin_area_Bar_Diff(
  df = DT_combined,
  out_file = "results_rootzone/Figures_till2022/main_PSI/F4_NDVI_PSIbin_area_Bar_Diff.png",
  bin_width = 50,
  min_count = 1000)
