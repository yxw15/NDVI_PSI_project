setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

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

# ---- F1 ----
### F1 - time serires NDVI and PSI; linear regression between NDVI with PSI ###
df_average_year_species <- function(df_all) {
  stopifnot(all(c("year", "species", "Quantiles", "soil_water_potential") %in% names(df_all)))
  
  df_all %>%
    group_by(year, species) %>%
    summarise(
      avg_quantile = mean(Quantiles, na.rm = TRUE),
      avg_psi      = mean(soil_water_potential, na.rm = TRUE),
      .groups = "drop"
    )
}

calc_stats_by_species <- function(df_species) {
  library(dplyr)
  library(broom)
  
  # correlation r + p
  corr_tbl <- df_species %>%
    group_by(species) %>%
    reframe({
      x <- avg_psi
      y <- avg_quantile
      ok <- is.finite(x) & is.finite(y)
      x <- x[ok]; y <- y[ok]
      
      if (length(x) >= 3 && sd(x) > 0 && sd(y) > 0) {
        ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
        tibble(r = unname(ct$estimate), p = ct$p.value, n = length(x))
      } else {
        tibble(r = NA_real_, p = NA_real_, n = length(x))
      }
    }) %>%
    ungroup()
  
  # linear model summary + coefficients
  lm_overall <- df_species %>%
    group_by(species) %>%
    reframe({
      d <- cur_data() %>% filter(is.finite(avg_quantile), is.finite(avg_psi))
      if (nrow(d) < 3 || sd(d$avg_psi) == 0) return(tibble())
      
      m <- lm(avg_quantile ~ avg_psi, data = d)
      g <- glance(m)
      tibble(
        model_F_stat  = g$statistic,
        model_p_value = g$p.value,
        r_squared     = g$r.squared,
        adj_r_squared = g$adj.r.squared,
        df_residual   = g$df.residual
      )
    }) %>% ungroup()
  
  lm_coef <- df_species %>%
    group_by(species) %>%
    dplyr::reframe({
      d <- cur_data() %>% filter(is.finite(avg_quantile), is.finite(avg_psi))
      if (nrow(d) < 3 || sd(d$avg_psi) == 0) return(tibble())
      
      m <- lm(avg_quantile ~ avg_psi, data = d)
      tidy(m) %>%
        mutate(Coefficient = recode(term, `(Intercept)` = "a", `avg_psi` = "b")) %>%
        dplyr::select(Coefficient, Estimate = estimate, p_value = p.value)
    }) %>% ungroup()
  
  list(corr = corr_tbl, overall = lm_overall, coef = lm_coef)
}

plot_F1_time_series_NDVI_PSI <- function(
    df_all,
    output_path = NULL,
    years = 2003:2022,
    months = NULL,  # optional: e.g. c("July","August") if df_all has month column
    species_order = c("Oak", "Beech", "Spruce", "Pine"),
    cb_palette = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    base_theme = theme_bw(),
    legend_position = "bottom",
    width = 12,
    height = 14,
    dpi = 300
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
  })
  
  # --- aggregate to yearly means by species
  df_species <- df_all %>%
    filter(year %in% years) %>%
    { if (!is.null(months) && "month" %in% names(.)) dplyr::filter(., month %in% months) else . } %>%
    group_by(year, species) %>%
    summarise(
      avg_quantile = mean(Quantiles, na.rm = TRUE),
      avg_psi      = mean(soil_water_potential, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(species = factor(species, levels = species_order))
  
  x_scale <- scale_x_continuous(
    breaks = seq(min(years), max(years), by = 2),
    limits = c(min(years), max(years))
  )
  
  # (a) NDVI quantile time series
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_quantile, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(x = "year", y = expression(atop("NDVI quantiles", "(rank)")), color = "") +
    scale_color_manual(values = cb_palette) +
    x_scale + base_theme +
    ggtitle("(a)") +
    theme(
      plot.title.position = "panel",
      plot.title = element_text(size = 18, hjust = 0.0, vjust = -9, margin = margin(b = 5))
    )
  
  # (b) PSI time series
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(x = "year", y = expression(atop("soil water potential", "(kPa)"))) +
    scale_color_manual(values = cb_palette) +
    x_scale + base_theme +
    guides(color = "none") +
    ggtitle("(b)") +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  # (c) correlation: NDVI quantile ~ PSI
  p3 <- ggplot(df_species, aes(x = avg_psi, y = avg_quantile, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "soil water potential (kPa)", y = expression(atop("NDVI quantiles", "(rank)"))) +
    scale_color_manual(values = cb_palette) +
    base_theme +
    guides(color = "none") +
    ggtitle("(c)") +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  final_plot <- (p1 / p2 / p3) +
    plot_layout(guides = "collect") &
    theme(legend.position = legend_position)
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = final_plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  return(final_plot)
}

p <- plot_F1_time_series_NDVI_PSI(
  df_all = data,
  output_path = "results_rootzone/Figures_till2022/main_PSI/F1_time_series_NDVI_PSI_species_till2022.png",
  base_theme = base_theme,
  legend_position = "top"
)
print(p)

# ---- F2 ----
### F2 - NDVI-PSIbin with CI ###
NDVI_PSIbin <- function(df, bin_width = 50) {
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  # Bin the soil water potential values
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute statistics for each species and bin, including filtering out bins with count < 2000
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value  = mean(Quantiles, na.rm = TRUE),
      # bin_median = median(soil_water_potential, na.rm = TRUE),  # x-value for plotting
      bin_mean   = mean(soil_water_potential, na.rm = TRUE),  
      count      = n(),
      .groups    = "drop"
    ) %>%
    left_join(species_totals, by = "species") %>%
    filter(count >= 1000) %>%
    mutate(percentage = count / total_pixels) 
  return(meanNDVI_PSIbin_species)
}

NDVI_PSIbin_df <- NDVI_PSIbin(data)
NDVI_PSIbin_df %>%
  group_by(species) %>%
  summarise(
    n_bins     = n(),          # number of bins (nrow per species)
    min_pixels = min(count),   # minimum pixels across bins
    .groups = "drop"
  ) %>%
  print()

ribbon_line_key <- function(data, params, size) {
  grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(
        col  = NA,
        fill = scales::alpha(data$colour %||% data$fill, 0.25)
      )
    ),
    grid::linesGrob(
      x = c(0.08, 0.92), y = c(0.5, 0.5),
      gp = grid::gpar(
        col  = data$colour %||% data$color,
        lwd  = (data$linewidth %||% data$size %||% 0.5) * .pt,
        lty  = data$linetype %||% 1,
        alpha = 1
      )
    )
  )
}

### Delta method use bin_mean for CI
plot_NDVI_PSI_exp_linear_slope_coeff <- function(data, output_figure) {
  
  # 1. Data Preparation
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(nlme)
  library(car)
  library(MASS)
  library(scales)
  library(latex2exp)
  
  value_col <- "avg_value"
  linear_species <- c("Beech")
  exp_species <- c("Oak", "Spruce", "Pine")
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce"= "#009E73", "Pine"  = "#F0E442")
  
  x_col <- "bin_mean"
  data <- data %>% mutate(x = .data[[x_col]])
  
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  # 2. Model Fitting
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  models_linear <- list()
  for (sp in linear_species) {
    models_linear[[sp]] <- lm(avg_value ~ x, data = data_clean %>% filter(species == sp))
  }
  
  models_exp <- list()
  for (sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    models_exp[[sp]] <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), data = sp_data, start = start_list, control = ctrl)
    }, error = function(e) NULL)
  }
  
  # 3. Prediction and SYMMETRICAL Confidence Intervals
  pred_linear_ci <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if (is.null(mod)) return(NULL)
    
    # UPDATED: use observed bin_mean (or bin_median fallback via x_col), not a seq grid
    x_seq <- data_clean %>%
      filter(species == sp) %>%
      pull(.data[[x_col]]) %>%
      unique() %>%
      sort()
    
    pr   <- predict(mod, newdata = data.frame(x = x_seq), se.fit = TRUE)
    crit <- qt(0.975, df = df.residual(mod))
    
    data.frame(
      species = sp,
      x = x_seq,
      fit = pr$fit,
      lwr = pr$fit - crit * pr$se.fit,
      upr = pr$fit + crit * pr$se.fit
    )
  })
  
  pred_exp_ci <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    
    # UPDATED: use observed bin_mean (or bin_median fallback via x_col), not a seq grid
    x_seq <- data_clean %>%
      filter(species == sp) %>%
      pull(.data[[x_col]]) %>%
      unique() %>%
      sort()
    
    ci_points <- lapply(x_seq, function(val) {
      dm <- deltaMethod(mod, paste0("a + b * exp(c * ", val, ")"))
      data.frame(x = val, fit = dm$Estimate, se = dm$SE)
    }) %>% bind_rows()
    
    crit <- qt(0.975, df = summary(mod)$df[2])
    
    data.frame(
      species = sp,
      x = ci_points$x,
      fit = ci_points$fit,
      lwr = ci_points$fit - (crit * ci_points$se),
      upr = ci_points$fit + (crit * ci_points$se)
    )
  })
  
  pred_ci_all <- bind_rows(c(pred_linear_ci, pred_exp_ci))
  pred_ci_all$species <- factor(pred_ci_all$species, levels = levels(data$species))
  
  # 4. Plotting Panel A (REMOVE ribbon legend)
  p_line <- ggplot() +
    geom_point(
      data = data_clean,
      aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
      alpha = 0.4
    ) +
    geom_ribbon(
      data = pred_ci_all,
      aes(x = x, ymin = lwr, ymax = upr, fill = species),
      alpha = 0.25,
      show.legend = FALSE
    ) +
    geom_line(
      data = pred_ci_all,
      aes(x = x, y = fit, color = species),
      linewidth = 1.2,
      key_glyph = ribbon_line_key
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", alpha = 0.5) +
    annotate(
      "text",
      x = min(data_clean$x, na.rm = TRUE),
      y = threshold,
      label = "median",
      hjust = -0.1,
      vjust = -0.3,
      size = 6
    ) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette,  name = "") +
    scale_shape_manual(values = c("Oak"=16,"Beech"=17,"Spruce"=15,"Pine"=18), name = "") +
    scale_size_continuous(range = c(1, 8), labels = scales::percent_format(), name = "") +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(color = "grey80", shape = 16, alpha = 1)
      ),
      color = guide_legend(
        order = 2,
        override.aes = list(shape = NA, linewidth = 1.2)
      ),
      fill  = "none",
      shape = "none"
    ) +
    base_theme +
    theme(
      legend.position = "none",
      plot.title.position = "plot", 
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1),
    )
  
  
  # 5. Stats Analysis for B & C
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    summ <- summary(mod)$coefficients
    data.frame(
      species = sp,
      x50 = (threshold - summ[1, 1]) / summ[2, 1],
      abs_slope = abs(summ[2, 1]),
      se = summ[2, 2],
      p_val = summ[2, 4],
      r_squared = summary(mod)$r.squared
    )
  })
  
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    
    cc <- coef(mod)
    x50_v   <- log((threshold - cc["a"]) / cc["b"]) / cc["c"]
    slope50 <- cc["c"] * (threshold - cc["a"])
    dm      <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"))
    r2      <- 1 - sum(residuals(mod)^2) /
      sum((data_clean$avg_value[data_clean$species == sp] -
             mean(data_clean$avg_value[data_clean$species == sp]))^2)
    
    data.frame(
      species = sp,
      x50 = as.numeric(x50_v),
      abs_slope = abs(slope50),
      se = dm$SE,
      p_val = 2 * (1 - pt(abs(slope50 / dm$SE), summary(mod)$df[2])),
      r_squared = r2
    )
  })
  
  stats_all <- bind_rows(c(lin_stats, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = TeX("$\\Psi_{\\,\\, Q_{50}}$(kPa)")) +
    base_theme + theme(legend.position = "none") +
    ggtitle("(b)") +
    theme(
      legend.position = "none",
      plot.title.position = "plot", 
      plot.title = element_text(size = 18, hjust = 0.25, vjust = 1),
    )
  
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")),
                  y = abs_slope/2)) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = TeX("$|\\lambda\\, Q_{50}(\\Psi_{\\, \\, soil})|$")) +
    ggtitle("(c)") +
    base_theme +
    theme(
      legend.position = "none",
      plot.title.position = "plot", 
      plot.title = element_text(size = 18, hjust = 0.25, vjust = 1),
    )
  
  # Final Layout and Save
  # Final Layout and Save
  final_plot <- (p_line + (p_x50 / p_slope)) + 
    plot_layout(widths = c(2, 1), guides = "collect") & 
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(filename = output_figure, plot = final_plot, width = 12, height = 8, dpi = 300)
  
  print("Process Complete. Symmetrical ribbons applied via Delta Method. Ribbon legend removed.")
}

plot_NDVI_PSI_exp_linear_slope_coeff(data, 
                                     "results_rootzone/Figures_till2022/main_PSI/F2_NDVI_Q_PSIbin_exp_linear_slope.png")

# ---- F3 ----
## F3 - PSI vs PSI ###
##### data processing 
# prepare_psi_summary_data <- function(
#     years     = 2003:2022,
#     months    = c("July", "August"),
#     species   = c("Oak", "Beech", "Spruce", "Pine"),
#     base_dir  = "results_monthly_rootzone",
#     bin_width = 50,
#     min_count = 1000,
#     output_path = "results_rootzone/Data/psi_summary_pairwise.RData"
# ) {
#   suppressPackageStartupMessages({
#     library(terra)
#     library(tidyverse)
#   })
#   
#   # 1) Helper: Read Rasters 
#   read_psi <- function(year, month, species_name) {
#     file_path <- file.path(base_dir, month, species_name, paste0("psi_", year, ".tif"))
#     if(!file.exists(file_path)) return(NULL) # Safety check
#     
#     r <- rast(file_path)
#     df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
#     colnames(df)[3] <- "soil_water_potential"
#     
#     df %>%
#       mutate(
#         species = species_name,
#         month   = month,
#         year    = year
#       )
#   }
#   
#   message("Reading rasters...")
#   psi_all <- expand.grid(
#     year = years, month = months, species = species, stringsAsFactors = FALSE
#   ) %>%
#     pmap_dfr(~ read_psi(..1, ..2, ..3))
#   
#   # 2) Pivot Wide for Pairwise Processing 
#   message("Processing pairwise joins...")
#   psi_wide <- psi_all %>%
#     pivot_wider(names_from = species, values_from = soil_water_potential)
#   
#   # 3) Long + Pairwise Join 
#   psi_long <- psi_wide %>%
#     pivot_longer(cols = all_of(species), names_to = "species", values_to = "psi_value")
#   
#   psi_pairs <- psi_long %>%
#     inner_join(
#       psi_long,
#       by = c("x", "y", "month", "year"),
#       suffix = c("_x", "_y"),
#       relationship = "many-to-many"
#     ) %>%
#     filter(species_x != species_y) %>%
#     mutate(
#       species_x = factor(species_x, levels = species),
#       species_y = factor(species_y, levels = species)
#     )
#   
#   # 4) Binning and Summarising 
#   message("Binning data...")
#   psi_summary <- psi_pairs %>%
#     group_by(species_x) %>%
#     mutate(
#       PSI_bin = cut(
#         psi_value_x,
#         breaks = seq(
#           floor(min(psi_value_x, na.rm = TRUE)),
#           ceiling(max(psi_value_x, na.rm = TRUE)),
#           by = bin_width
#         ),
#         include.lowest = TRUE,
#         right = FALSE
#       )
#     ) %>%
#     group_by(species_x, species_y, PSI_bin) %>%
#     summarise(
#       bin_mean   = mean(psi_value_x, na.rm = TRUE),
#       mean_PSI_y = mean(psi_value_y, na.rm = TRUE),
#       count      = n(),
#       .groups    = "drop"
#     ) %>%
#     filter(count >= min_count) %>%
#     group_by(species_x, species_y) %>%
#     mutate(percentage = count / sum(count, na.rm = TRUE)) %>%
#     ungroup()
#   
#   # Save Results 
#   if (!dir.exists(dirname(output_path))) {
#     dir.create(dirname(output_path), recursive = TRUE)
#   }
#   
#   save(psi_summary, file = output_path)
#   message(paste("Process Complete. Data saved to:", output_path))
#   
#   return(psi_summary)
# }
#
# prepare_psi_summary_data(years = 2003:2022)
##### end 

load("results_rootzone/Data/psi_summary_pairwise.RData") 

plot_PSI_PSI_pairwise_GAM <- function(
    df_summary,
    k = 6,
    smooth_method = "REML",
    point_alpha = 0.4,
    smooth_alpha = 0.25,
    smooth_linewidth = 1,
    abline_lwd = 0.8,
    cb_palette = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    cb_shapes  = c(Oak=16, Beech=17, Spruce=15, Pine=18),
    save_path = NULL,
    save_width = 12,
    save_height = 8,
    save_dpi = 300
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(mgcv)
    library(scales)
    library(patchwork)
    library(dplyr)
  })
  
  # Ensure consistent order
  df_summary <- df_summary %>%
    mutate(
      species_x = factor(species_x, levels = names(cb_palette)),
      species_y = factor(species_y, levels = names(cb_palette))
    )
  
  p <- ggplot(
    df_summary,
    aes(x = bin_mean, y = mean_PSI_y)
  ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", alpha = 0.5, linewidth = abline_lwd
    ) +
    
    # Points: show ONLY size legend (no species legend from points)
    geom_point(
      aes(color = species_y, shape = species_y, size = percentage),
      alpha = point_alpha,
      show.legend = c(color = FALSE, shape = FALSE, size = TRUE)
    ) +
    
    ## --- (A) Ribbon layer: THIS drives species legend (fill) ---
    geom_smooth(
      aes(group = species_y, fill = species_y),
      method  = "gam",
      formula = y ~ s(x, k = k),
      method.args = list(method = smooth_method),
      se = TRUE,
      alpha = smooth_alpha,
      linewidth = 0,                # ribbon only
      color = NA,                   # no line
      show.legend = TRUE            # <-- IMPORTANT: allow legend
    ) +
    
    # --- (B) Line layer: drawn in plot, but not in legend ---
    geom_smooth(
      aes(group = species_y, color = species_y),
      method  = "gam",
      formula = y ~ s(x, k = k),
      method.args = list(method = smooth_method),
      se = FALSE,
      linewidth = smooth_linewidth,
      show.legend = FALSE
    ) +
    
    facet_wrap(~ species_x, ncol = 2, scales = "free") +
    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette,  name = "") +
    scale_shape_manual(values = cb_shapes,  name = "") +
    
    scale_size_continuous(
      name   = "",
      range  = c(1, 8),
      labels = percent_format(accuracy = 1)
    ) +
    
    labs(
      x = "soil water potential of panel species (kPa)",
      y = "soil water potential of line species (kPa)"
    ) +
    
    base_theme +
    
    # Legend control:
    #  - size legend first (grey points)
    #  - fill legend second (smooth key glyph shows ribbon+line)
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          colour = "grey70",
          shape = 16,
          fill = NA,
          linewidth = 0,
          alpha = 1
        )
      ),
      fill = guide_legend(
        order = 2,
        override.aes = list(
          # make sure legend key shows as "smooth" not points
          linewidth = smooth_linewidth,
          alpha = 0.4,
          colour = unname(cb_palette),  # line color per species
          shape = NA,
          size = 3
        )
      ),
      color = "none",
      shape = "none"
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.background = element_blank()
    )
  
  print(p)
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(save_path, p, width = save_width, height = save_height, dpi = save_dpi, bg = "white")
  }
  
  invisible(p)
}


out_plot <- plot_PSI_PSI_pairwise_GAM(
  df_summary = psi_summary, # Input the data you created/loaded
  save_path = "results_rootzone/Figures_till2022/main_PSI/F3_PSI_PSI_pairwise_4panel_GAM.png"
)

# ---- F4 Line ----
## F4 NDVI-PSIbin Germany vs Species-location Line Difference (option 1) ###
load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")
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

plot_F4_NDVI_PSIbin_area_Line_Diff <- function(
    df,
    out_file = "results_rootzone/main_PSI/F4_NDVI_PSIbin_area_Line_Diff.png",
    bin_width = 50,
    min_count = 1000,
    k_left  = 5,
    k_right = 3,
    bs = "cs",
    width = 16,
    height = 8,
    dpi = 300,
    sp_levels = c("Oak", "Beech", "Spruce", "Pine"),
    sp_cols   = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    sp_shapes = c(Oak=16, Beech=17, Spruce=15, Pine=18),
    show_germany = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(mgcv)
    library(scales)
    library(patchwork)
  })
  
  d_all <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  d_sp <- d_all %>% filter(area == "species")
  d_de <- d_all %>% filter(area == "Germany")
  
  # --- (a) left ---
  p_left <- ggplot(d_sp, aes(x = bin_mean, y = avg_value)) +
    
    # Germany reference: points + CI ribbon + line (all grey, NOT in legend)
    {if (show_germany) list(
      geom_point(
        data = d_de,
        aes(x = bin_mean, y = avg_value, size = percentage),
        colour = "grey50",
        alpha = 0.35,
        inherit.aes = FALSE,
        show.legend = FALSE
      ),
      geom_smooth(
        data = d_de,
        aes(x = bin_mean, y = avg_value, group = species),
        method  = "gam",
        formula = y ~ s(x, k = k_left, bs = bs),
        se = TRUE,
        alpha = 0.18,
        colour = NA,
        fill = "grey70",
        linewidth = 0,
        inherit.aes = FALSE,
        show.legend = FALSE
      ),
      geom_smooth(
        data = d_de,
        aes(x = bin_mean, y = avg_value, group = species),
        method  = "gam",
        formula = y ~ s(x, k = k_left, bs = bs),
        se = FALSE,
        colour = "grey40",
        linewidth = 1.0,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    )} +
    
    # species points (drives colour/shape + size legend)
    geom_point(
      aes(colour = species, shape = species, size = percentage),
      alpha = 0.40
    ) +
    
    # species CI ribbon (hidden from legend)
    geom_smooth(
      aes(fill = species, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_left, bs = bs),
      se = TRUE,
      alpha = 0.20,
      linewidth = 0,
      show.legend = FALSE
    ) +
    
    # species line (in legend)
    geom_smooth(
      aes(colour = species, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_left, bs = bs),
      se = FALSE,
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    
    facet_wrap(~ species, nrow = 2) +
    scale_colour_manual(values = sp_cols,  name = "") +
    scale_fill_manual(values = sp_cols,    name = "") +
    scale_shape_manual(values = sp_shapes, name = "") +
    
    # single size legend we keep (panel a)
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = percent_format(accuracy = 1)
    ) +
    
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)",
      title = "(a)"
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 0.25, linewidth = 1.2)),
      fill   = "none"
    ) +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  # --- (b) right ---
  df_diff <- NDVI_diff_species_minus_Germany_PSIbin(df, bin_width, min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  p_right <- ggplot(df_diff, aes(x = bin_mean, y = diff_NDVI)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    
    geom_point(
      aes(colour = species, shape = species, size = pct_bin),
      alpha = 0.40
    ) +
    
    geom_smooth(
      aes(fill = species, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_right, bs = bs),
      se = TRUE,
      alpha = 0.25,
      linewidth = 0,
      show.legend = FALSE
    ) +
    
    geom_smooth(
      aes(colour = species, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_right, bs = bs),
      se = FALSE,
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    
    scale_colour_manual(values = sp_cols,  name = "") +
    scale_fill_manual(values = sp_cols,    name = "") +
    scale_shape_manual(values = sp_shapes, name = "") +
    
    # keep same size scale, BUT hide it in panel (b)
    scale_size_continuous(
      name   = "pixels per bin (%)",   # keep consistent (doesn't matter if hidden)
      range  = c(1, 8),
      labels = percent_format(accuracy = 1)
    ) +
    
    labs(
      x = "soil water potential (kPa)",
      y = expression(Delta[NDVI]~"(species - Germany)"),
      title = "(b)"
    ) +
    guides(
      colour = guide_legend(override.aes = list(alpha = 0.25, linewidth = 1.2)),
      fill   = "none",
      size   = "none"   # <<< removes the 2nd percentage legend
    ) +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  p <- (p_left | p_right) +
    plot_layout(widths = c(1.25, 1), guides = "collect") &
    theme(
      legend.position = "bottom"
    ) +
    legend_theme_F2()
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
  invisible(p)
}

plot_F4_NDVI_PSIbin_area_Line_Diff(
  df = DT_combined,
  out_file = "results_rootzone/main_PSI/F4_NDVI_PSIbin_area_Line_Diff.png",
  bin_width = 50,
  min_count = 1000,
  k_left = 5,
  k_right = 3
)

# ---- F4 Bar ----
load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")

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
