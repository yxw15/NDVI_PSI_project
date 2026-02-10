plot_PSI_PSI_pairwise_GAM <- function(
    years   = 2003:2022,
    months  = c("July", "August"),
    species = c("Oak", "Beech", "Spruce", "Pine"),
    base_dir = "results_monthly_rootzone",
    bin_width = 50,
    min_count = 1000,
    k = 6,                  # GAM basis dimension
    smooth_method = "REML", # mgcv smoothing parameter estimation
    point_size = 2,
    point_alpha = 0.4,
    smooth_alpha = 0.25,
    smooth_linewidth = 1.2,
    abline_lwd = 0.8,
    cb_palette = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    cb_shapes  = c(Oak=16, Beech=17, Spruce=15, Pine=18),
    save_path = NULL,
    save_width = 12,
    save_height = 8,
    save_dpi = 300
) {
  suppressPackageStartupMessages({
    library(terra)
    library(tidyverse)
    library(mgcv)
    library(scales)
    library(patchwork)
  })
  
  # ---- Theme + Guides (Figure-2 style) ----
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
  
  base_guides <- guides(
    size  = guide_legend(order = 1, override.aes = list(alpha = 1)),
    color = guide_legend(order = 2, override.aes = list(linewidth = 1.2, alpha = 1)),
    shape = guide_legend(order = 2)
  )
  
  # ---- 1) read rasters ----
  read_psi <- function(year, month, species_name) {
    file_path <- file.path(base_dir, month, species_name, paste0("psi_", year, ".tif"))
    r <- rast(file_path)
    
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    colnames(df)[3] <- "soil_water_potential"
    
    df %>%
      mutate(
        species = species_name,
        month   = month,
        year    = year
      )
  }
  
  psi_all <- expand.grid(
    year = years, month = months, species = species, stringsAsFactors = FALSE
  ) %>%
    pmap_dfr(~ read_psi(..1, ..2, ..3))
  
  # ---- 2) wide ----
  psi_wide <- psi_all %>%
    pivot_wider(names_from = species, values_from = soil_water_potential)
  
  # ---- 3) long + pairwise join ----
  psi_long <- psi_wide %>%
    pivot_longer(cols = all_of(species), names_to = "species", values_to = "psi_value")
  
  psi_pairs <- psi_long %>%
    inner_join(
      psi_long,
      by = c("x", "y", "month", "year"),
      suffix = c("_x", "_y"),
      relationship = "many-to-many"
    ) %>%
    filter(species_x != species_y) %>%
    mutate(
      species_x = factor(species_x, levels = species),
      species_y = factor(species_y, levels = species)
    )
  
  # ---- 4) bin x PSI and summarise ----
  psi_summary <- psi_pairs %>%
    group_by(species_x) %>%
    mutate(
      PSI_bin = cut(
        psi_value_x,
        breaks = seq(
          floor(min(psi_value_x, na.rm = TRUE)),
          ceiling(max(psi_value_x, na.rm = TRUE)),
          by = bin_width
        ),
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    group_by(species_x, species_y, PSI_bin) %>%
    summarise(
      bin_mean   = mean(psi_value_x, na.rm = TRUE),
      mean_PSI_y = mean(psi_value_y, na.rm = TRUE),
      count      = n(),
      .groups    = "drop"
    ) %>%
    filter(count >= min_count) %>%
    group_by(species_x, species_y) %>%
    mutate(percentage = count / sum(count, na.rm = TRUE)) %>%
    ungroup()
  
  # ---- 5) plotting (Figure-2 legend logic) ----
  # Key idea:
  #   - legend comes from POINTS (color+shape) + SIZE
  #   - ribbon exists but does NOT create its own legend
  #   - smooth line uses color mapping so it appears in the species legend
  p <- ggplot(
    psi_summary,
    aes(x = bin_mean, y = mean_PSI_y, color = species_y, shape = species_y)
  ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", alpha = 0.5, linewidth = abline_lwd
    ) +
    
    # points
    geom_point(aes(size = percentage), alpha = point_alpha) +
    
    # ribbon only (no legend)
    geom_smooth(
      aes(group = species_y, fill = species_y),
      method  = "gam",
      formula = y ~ s(x, k = k),
      method.args = list(method = smooth_method),
      se  = TRUE,
      alpha = smooth_alpha,
      linewidth = 0,
      show.legend = FALSE
    ) +
    
    # line only (drives "species" legend via color)
    geom_smooth(
      aes(group = species_y),
      method  = "gam",
      formula = y ~ s(x, k = k),
      method.args = list(method = smooth_method),
      se  = FALSE,
      linewidth = smooth_linewidth,
      show.legend = FALSE   # <- keep legend identical to Figure 2 (species legend from points)
    ) +
    
    facet_wrap(~ species_x, ncol = 2, scales = "free") +
    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette,  name = "") +
    scale_shape_manual(values = cb_shapes,  name = "") +
    scale_size_continuous(
      range  = c(1, 8),
      labels = percent_format(accuracy = 1),
      name   = "pixels per bin (%)"
    ) +
    
    labs(
      x = "soil water potential of panel species (kPa)",
      y = "soil water potential of line species (kPa)"
    ) +
    
    base_theme +
    base_guides +
    guides(fill = "none") +
    theme(legend.position = "bottom")
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(
      filename = save_path,
      plot = p,
      width = save_width,
      height = save_height,
      dpi = save_dpi,
      bg = "white"
    )
  }
  
  invisible(list(plot = p, psi_summary = psi_summary))
}

out <- plot_PSI_PSI_pairwise_GAM(
  years = 2003:2022,
  months = c("July","August"),
  base_dir = "results_monthly_rootzone",
  bin_width = 50,
  min_count = 1000,
  k = 6,
  save_path = "results_rootzone/Figures_till2022/main_PSI/F3_PSI_PSI_pairwise_4panel_GAM.png"
)
print(out$plot)