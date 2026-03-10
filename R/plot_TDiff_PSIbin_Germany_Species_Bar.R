load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")

library(tidyverse)

data <- combined %>% filter(month %in% c("July", "August")) %>% filter(year < 2023)
data <- data %>% filter(Quantiles > 0)
data <- na.omit(data)

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


TDiff_PSIbin_area <- function(df, bin_width = 50, min_count = 1000) {
  library(dplyr)
  
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
    ) %>%
    filter(!is.na(PSI_bin))  # <-- important
  
  totals <- df_binned %>%
    group_by(area, species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  df_binned %>%
    group_by(area, species, PSI_bin) %>%
    summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      avg_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE),
      count = n(),
      .groups = "drop"
    ) %>%
    left_join(totals, by = c("area", "species")) %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    tidyr::drop_na()
}


plot_TDiff_PSIbin_byArea <- function(data, figure_output = NULL,
                                     bin_width = 50, min_count = 1000,
                                     conf_level = 0.95) {
  
  library(tidyverse)
  library(scales)
  library(mgcv)
  
  df <- TDiff_PSIbin_area(data,
                          bin_width = bin_width,
                          min_count = min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      area    = factor(area, levels = c("Germany", "species"))
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  xcol <- "bin_mean"
  ycol <- "avg_transpiration_deficit"
  
  alpha <- 1 - conf_level
  
  # -------------------------
  # Fit GAM + predict per area × species
  # -------------------------
  pred_df <- df %>%
    group_by(area, species) %>%
    group_modify(~{
      sdf <- .x %>% select(all_of(c(xcol, ycol)))
      
      if (nrow(sdf) < 5) return(tibble())
      
      # Fit GAM with k = 3
      fit <- mgcv::gam(
        reformulate(termlabels = "s(bin_mean, k = 3)", response = ycol),
        data = sdf,
        method = "REML"
      )
      
      xr <- range(sdf[[xcol]], na.rm = TRUE)
      xseq <- seq(xr[1], xr[2], length.out = 200)
      nd <- tibble(bin_mean = xseq)
      
      pr <- predict(fit, newdata = nd, se.fit = TRUE)
      
      # For mgcv GAMs, using normal approx is standard
      zcrit <- qnorm(1 - alpha/2)
      
      nd %>%
        mutate(
          predicted = as.numeric(pr$fit),
          se = as.numeric(pr$se.fit),
          lower = predicted - zcrit * se,
          upper = predicted + zcrit * se
        )
    }) %>%
    ungroup()
  
  # -------------------------
  # Plot: TWO PANELS by area + confidence bands
  # -------------------------
  plot_clean <- ggplot() +
    
    geom_ribbon(
      data = pred_df,
      aes(x = bin_mean, ymin = lower, ymax = upper,
          fill = species, group = interaction(area, species)),
      alpha = 0.18,
      color = NA
    ) +
    
    geom_point(
      data = df,
      aes(x = .data[[xcol]],
          y = .data[[ycol]],
          color = species,
          shape = species,
          size = percentage),
      alpha = 0.7
    ) +
    
    geom_line(
      data = pred_df,
      aes(x = bin_mean,
          y = predicted,
          color = species),
      linewidth = 1
    ) +
    
    facet_wrap(
      ~ area,
      ncol = 2,
      labeller = as_labeller(c(
        Germany = "(a) Germany",
        species = "(b) species"
      ))
    ) +
    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17,
                                  "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit (mm)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
  
  print(plot_clean)
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = plot_clean,
           width = 10, height = 5, dpi = 300)
  }
  
  invisible(plot_clean)
}

# plot_TDiff_PSIbin_byArea(combined_all, "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_area_gam3_till2022.png")
plot_TDiff_PSIbin_byArea(DT_combined, "results_rootzone/test/TDiff_PSIbin_area_gam3_till2022.png")


# ---- TDiff_PSIbin Germany vs Species Functions  ----

TDiff_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  species_totals <- df %>%
    group_by(area, species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
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
      avg_value = mean(transpiration_deficit, na.rm = TRUE),
      count     = n(),
      .groups   = "drop"
    ) %>%
    left_join(species_totals, by = c("area", "species")) %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    select(area, species, PSI_bin, bin_mean, avg_value, count, total_pixels, percentage)
}

TDiff_diff_species_minus_Germany_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
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
      sp_mean  = mean(transpiration_deficit, na.rm = TRUE),
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
      de_mean = mean(transpiration_deficit, na.rm = TRUE),
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
      diff_TDiff = sp_mean - de_mean,
      pct_bin    = pmin(sp_pct, de_pct)
    ) %>%
    select(species, PSI_bin, bin_mean, sp_mean, de_mean, sp_n, de_n, sp_pct, de_pct, pct_bin, diff_TDiff)
}


# ---- TDiff_PSIbin Germany vs Species (Bar difference) ----

plot_TDiff_PSIbin_area_Bar_Diff <- function(
    df,
    out_file = "results_rootzone/main_PSI/F4_TDiff_PSIbin_area_Bar_Diff.png",
    bin_width = 50,
    min_count = 1000,
    k_left = 3,
    gam_method = "REML",
    # bs_left = "tp",
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
  
  d_all <- TDiff_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  d_sp <- d_all %>% filter(area == "species")
  d_de <- d_all %>% filter(area == "Germany")
  
  legend_cols   <- c(sp_cols, "Germany" = germany_grey)
  legend_fills  <- c(sp_cols, "Germany" = "grey70")
  legend_shapes <- c(sp_shapes, "Germany" = 32)
  
  p_left <- ggplot(d_sp, aes(bin_mean, avg_value)) +
    # --- Germany: Points + Smooth with CI ---
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
    
    # --- Species: Points + Smooth with CI ---
    geom_point(aes(colour = species, shape = species, size = percentage), alpha = 0.4) +
    geom_smooth(
      aes(colour = species, fill = species),
      method = "gam",
      formula = y ~ s(x, k = k_left),
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
    labs(x = "soil water potential (kPa)", y = "TDiff", title = "(a)") +
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
  df_diff <- TDiff_diff_species_minus_Germany_PSIbin(df, bin_width, min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  df_bar <- df_diff %>%
    group_by(species) %>%
    summarise(
      mean_unweighted = mean(diff_TDiff, na.rm = TRUE),
      mean_weighted   = weighted.mean(diff_TDiff, w = pct_bin, na.rm = TRUE),
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
    labs(x = "", y = expression(Delta*" TDiff"), title = "(b) weighted mean") +
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
    labs(x = "", y = expression(Delta*" TDiff"), title = "(c) unweighted mean") +
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

# ---- run ----
plot_TDiff_PSIbin_area_Bar_Diff(
  df = DT_combined,
  out_file = "results_rootzone/Figures_till2022/SI_PSI/SI_TDiff_PSIbin_area_Bar_Diff.png",
  bin_width = 50,   # 50 kPa bins
  min_count = 1000
)


# ---- TDiff_PSIbin Germany vs Species (Bar difference) 4 Panels ----
# ---- TDiff_PSIbin Germany vs Species (Bar difference) 4 Panels ----
plot_TDiff_PSIbin_area_4panel <- function(
    df,
    out_file = "results_rootzone/main_PSI/SI_TDiff_PSIbin_area_4panel.png",
    bin_width = 50,
    min_count = 1000,
    k_left = 3,
    gam_method = "REML",
    width = 12,
    height = 10,
    dpi = 300,
    sp_levels = c("Oak", "Beech", "Spruce", "Pine"),
    sp_cols   = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    sp_shapes = c(Oak=16, Beech=17, Spruce=15, Pine=18)
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(mgcv)
    library(scales)
    library(patchwork)
    library(latex2exp)
  })
  
  # ---- binned summaries (top panels) ----
  d_all <- TDiff_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  d_de <- d_all %>% filter(area == "Germany")
  d_sp <- d_all %>% filter(area == "species")
  
  # ---- shared y-limits (top) ----
  ylim_top <- range(c(d_de$avg_value, d_sp$avg_value), na.rm = TRUE)
  
  # IMPORTANT:
  # We remove Germany (grey) point+line entirely, and ensure only ONE size legend
  # by hiding size legend in panel (b).
  make_area_panel <- function(d_sub, panel_title, show_size_legend = TRUE) {
    
    ggplot(d_sub, aes(bin_mean, avg_value)) +
      
      # --- Species only: points + smooth with CI ---
      geom_point(
        aes(colour = species, shape = species, size = percentage),
        alpha = 0.4
      ) +
      geom_smooth(
        aes(colour = species, fill = species),
        method = "gam",
        formula = y ~ s(x, k = k_left),
        method.args = list(method = gam_method),
        se = TRUE,
        alpha = 0.25,
        linewidth = 1.2
      ) +
      
      scale_y_continuous(limits = ylim_top) +
      scale_colour_manual(values = sp_cols, name = "") +
      scale_fill_manual(values = sp_cols, name = "") +
      scale_shape_manual(values = sp_shapes, name = "") +
      scale_size_continuous(
        range = c(1, 8),
        labels = percent_format(accuracy = 1),
        name = ""
      ) +
      labs(x = "soil water potential (kPa)", y = NULL, title = panel_title) +
      base_theme +
      guides(
        # pixel-size legend: only keep it in ONE panel
        size = if (show_size_legend) guide_legend(
          order = 1,
          override.aes = list(
            shape = 16,
            colour = "grey80",
            alpha = 1,
            fill = NA,
            linewidth = 0
          )
        ) else "none",
        
        # species line legend (with implied ribbon) — ONLY four species
        colour = guide_legend(
          order = 2,
          override.aes = list(
            linewidth = 1,
            alpha = 0.4,
            fill = sp_cols
          )
        ),
        
        # remove extra legends
        shape = "none",
        fill  = "none"
      ) +
      theme(
        plot.title.position = "plot",
        plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        legend.position = "bottom",
        legend.box = "horizontal"
      )
  }
  
  # ---- top row: (a) and (b) ----
  # show size legend ONLY once (panel a)
  p_a <- make_area_panel(d_de, "(a) Germany", show_size_legend = TRUE) +
    labs(y = TeX("transpiration deficit (mm)"))
  
  p_b <- make_area_panel(d_sp, "(b) species", show_size_legend = FALSE)
  
  # ---- bottom bars ----
  df_diff <- TDiff_diff_species_minus_Germany_PSIbin(df, bin_width, min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  df_bar <- df_diff %>%
    group_by(species) %>%
    summarise(
      mean_unweighted = mean(diff_TDiff, na.rm = TRUE),
      mean_weighted   = weighted.mean(diff_TDiff, w = pct_bin, na.rm = TRUE),
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
  
  # ---- shared y-limits (bottom) ----
  ylim_bot <- range(df_bar$value, na.rm = TRUE)
  
  p_c <- ggplot(
    df_bar %>% filter(type == "Weighted (by min area support)"),
    aes(x = species, y = value, fill = species)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) +
    scale_y_continuous(limits = ylim_bot) +
    scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = TeX("$\\Delta\\, T_d\\ (mm)$"), title = "(c) weighted mean") +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  p_d <- ggplot(
    df_bar %>% filter(type == "Unweighted"),
    aes(x = species, y = value, fill = species)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) +
    scale_y_continuous(limits = ylim_bot) +
    scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = NULL, title = "(d) unweighted mean") +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 1),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  # ---- assemble 2x2 ----
  p <- ((p_a | p_b) / (p_c | p_d)) +
    plot_layout(heights = c(2.0, 1.0), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  return(p)
}


plot_TDiff_PSIbin_area_4panel(
  df = DT_combined,
  out_file = "results_rootzone/Figures_till2022/SI_PSI/SI_TDiff_PSIbin_area_4panel.png",
  bin_width = 50,
  min_count = 1000
)
