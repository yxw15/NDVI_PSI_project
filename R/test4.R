## F4 NDVI-PSIbin Germany vs Species-location Bar Difference (option 2) ###
plot_F4_NDVI_PSIbin_area_Bar_Diff <- function(
    df,
    out_file = "results_rootzone/main_PSI/F4_NDVI_PSIbin_area_Bar_Diff.png",
    bin_width = 50,
    min_count = 1000,
    k_left = 6,
    bs_left = "cs",
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
  
  # Define specific colors and shapes for the legend including Germany
  # We use a trick: named vectors for the legend keys
  legend_cols <- c(sp_cols, "Germany" = germany_grey)
  legend_fills <- c(sp_cols, "Germany" = "grey70")
  legend_shapes <- c(sp_shapes, "Germany" = 32) # 32 is empty/no shape for Germany line
  
  # --- Inside plot_F4_NDVI_PSIbin_area_Bar_Diff ---
  
  p_left <- ggplot(d_sp, aes(bin_mean, avg_value)) +
    # --- Germany: Points, Ribbon, Line ---
    geom_point(data = d_de, aes(size = percentage), colour = germany_grey, alpha = 0.2, show.legend = FALSE) +
    
    # 1. Germany Smooth: Explicitly map color and fill to "Germany"
    geom_smooth(
      data = d_de, 
      aes(x = bin_mean, y = avg_value, group = species, colour = "Germany", fill = "Germany"),
      method = "gam", 
      formula = y ~ s(x, k = k_left, bs = bs_left),
      se = TRUE,        # This ensures the Germany confidence band is drawn
      alpha = 0.2, 
      linewidth = 1
    ) +
    
    # --- Species: Points, Ribbon, Line ---
    geom_point(aes(colour = species, shape = species, size = percentage), alpha = 0.4) +
    
    # 2. Species Smooth: Mapping fill ensures ribbons are drawn for all 4 species
    geom_smooth(
      aes(colour = species, fill = species, group = species),
      method = "gam", 
      formula = y ~ s(x, k = k_left, bs = bs_left),
      se = TRUE,        # This ensures species confidence bands are drawn
      alpha = 0.25, 
      linewidth = 1.2
    ) +
    
    facet_wrap(~ species, nrow = 2) +
    
    # Scale configuration
    scale_colour_manual(values = legend_cols, name = "") +
    scale_fill_manual(values = legend_fills, name = "") +
    scale_shape_manual(values = legend_shapes, name = "") +
    scale_size_continuous(range = c(1, 8), labels = percent_format(accuracy = 1), name = "") +
    
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)", title = "(a)") +
    
    base_theme +
    guides(
      # Legend 1: Percentage (Grey dots only)
      size = guide_legend(
        order = 1,
        override.aes = list(
          shape = 16, 
          colour = "grey70", 
          alpha = 1,
          fill = NA, 
          linewidth = 0
        )
      ),
      # Legend 2: Species + Germany (Force BOTH Line and Ribbon Box)
      colour = guide_legend(
        order = 2,
        override.aes = list(
          linewidth = 1,
          alpha = 0.4,
          # This is the critical part: mapping the background fill for all 5 entries
          fill = legend_fills  
        )
      ),
      # Merging shape and fill into the color guide
      shape = "none",
      fill = "none"
    ) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  # --- Bar Logic (Identical to your provided snippet) ---
  df_diff <- NDVI_diff_species_minus_Germany_PSIbin(df, bin_width, min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  df_bar <- df_diff %>%
    group_by(species) %>%
    summarise(
      mean_unweighted = mean(diff_NDVI, na.rm = TRUE),
      mean_weighted   = weighted.mean(diff_NDVI, w = pct_bin, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean_weighted, mean_unweighted), names_to = "type", values_to = "value") %>%
    mutate(type = recode(type, mean_weighted = "Weighted (by min area support)", mean_unweighted = "Unweighted"))
  
  p_right_top <- ggplot(df_bar %>% filter(type == "Weighted (by min area support)"), aes(x = species, y = value, fill = species)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) + scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = expression(Delta[paste("  ", NDVI)]), title = "(b) weighted mean") + base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  p_right_bottom <- ggplot(df_bar %>% filter(type == "Unweighted"), aes(x = species, y = value, fill = species)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) + scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = expression(Delta[paste("  ", NDVI)]), title = "(c) unweighted mean") + base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  # Combine: (a) on top, (b) and (c) in a third row underneath
  p_bottom <- (p_right_top | p_right_bottom) +
    plot_layout(widths = c(1, 1))
  
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
