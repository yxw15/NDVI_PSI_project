setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_wholeGermany.RData")
# load("results/Data/AllSpecies_AllMonths_rootzone.RData")
# 
# library(tidyverse)
# 
# data <- combined %>% 
#   filter(month %in% c("July", "August"),
#          Quantiles > 0,
#          year <= 2022) %>%
#   drop_na()
# 
# data$area <- "species"
# 
# DT_combined <- bind_rows(
#   DT_all,
#   data %>%
#     rename(NDVI = Quantiles) %>%     # match column name
#     select(
#       x, y, NDVI, soil_water_potential, transpiration_deficit,
#       species, month, year, area
#     )                               # drop rootzone + reorder
# )
# 
# save(DT_combined, file = "results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")

load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")

########### only species located data #######################
# -------------------------
# 1) BINNING: use bin_mean (mean soil_water_potential in each bin)
# -------------------------
NDVI_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # totals per species for percentage
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # global bin breaks so bins consistent across species
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
    group_by(species, PSI_bin) %>%
    summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),             # ✅ NEW
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = "drop"
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    select(species, PSI_bin, bin_mean, avg_value, count, total_pixels, percentage)
}

# -------------------------
# 2) Delta-method CI for exponential model: f(x)=a + b exp(c x)
# -------------------------
pred_exp_delta <- function(fit, xseq, conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  zcrit <- qnorm(1 - alpha/2)
  
  co <- coef(fit)
  V  <- vcov(fit)
  
  a <- unname(co["a"]); b <- unname(co["b"]); c <- unname(co["c"])
  ex <- exp(c * xseq)
  
  pred <- a + b * ex
  
  # gradient wrt (a,b,c): (1, exp(cx), b x exp(cx))
  se <- map_dbl(seq_along(xseq), \(i) {
    g <- c(1, ex[i], b * xseq[i] * ex[i])
    as.numeric(sqrt(t(g) %*% V %*% g))
  })
  
  tibble(
    x = xseq,
    pred = pred,
    lower = pred - zcrit * se,
    upper = pred + zcrit * se
  )
}

# -------------------------
# 3) Main plotting function
# -------------------------
plot_NDVI_PSI_exp_linear_slope_coeff <- function(data,
                                                 combined_coef_fig,
                                                 output_figure,
                                                 aic_barplot_fig,
                                                 bin_width = 50,
                                                 min_count = 1000,
                                                 conf_level = 0.95) {
  
  library(ggplot2)
  library(patchwork)
  library(car)
  
  # ---- process + bin
  data_b <- NDVI_PSIbin(data, bin_width = bin_width, min_count = min_count) %>%
    drop_na()
  
  value_col <- "avg_value"
  
  # model assignment
  linear_species <- c("Beech")
  exp_species <- c("Oak", "Spruce", "Pine")
  
  data_b <- data_b %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      x = bin_mean                                     # ✅ NEW: use bin_mean as x
    ) %>%
    filter(is.finite(x), !is.na(.data[[value_col]]))
  
  # palette
  cb_palette <- c("Oak" = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce" = "#009E73",
                  "Pine" = "#F0E442")
  
  # threshold
  threshold <- 11.5
  
  # -------------------------
  # AIC comparison (linear vs exponential) for each species
  # -------------------------
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  aic_df <- map_dfr(c(linear_species, exp_species), \(sp) {
    df_sp <- data_b %>% filter(species == sp)
    
    lm_mod <- tryCatch(lm(avg_value ~ x, data = df_sp), error = \(e) NULL)
    exp_mod <- tryCatch(
      nls(avg_value ~ a + b * exp(c * x),
          data = df_sp, start = start_list, control = ctrl),
      error = \(e) NULL
    )
    
    tibble(
      species = sp,
      AIC_linear = if (!is.null(lm_mod)) AIC(lm_mod) else NA_real_,
      AIC_exponential = if (!is.null(exp_mod)) AIC(exp_mod) else NA_real_
    )
  }) %>%
    mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))
  
  aic_long <- aic_df %>%
    pivot_longer(c(AIC_linear, AIC_exponential), names_to = "model", values_to = "AIC") %>%
    mutate(model = dplyr::recode(model,
                                 AIC_linear = "linear",
                                 AIC_exponential = "exponential"))
  
  p_aic <- ggplot(aic_long, aes(x = species, y = AIC, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("linear" = "dodgerblue", "exponential" = "orange"), name = "") +
    labs(title = "", x = "", y = "AIC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  ggsave(filename = aic_barplot_fig, plot = p_aic, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------
  # Fit models
  # -------------------------
  models_linear <- map(linear_species, \(sp) {
    df_sp <- data_b %>% filter(species == sp)
    lm(avg_value ~ x, data = df_sp)
  }) %>% set_names(linear_species)
  
  models_exp <- map(exp_species, \(sp) {
    df_sp <- data_b %>% filter(species == sp)
    if (nrow(df_sp) < 5) return(NULL)
    tryCatch(
      nls(avg_value ~ a + b * exp(c * x),
          data = df_sp, start = start_list, control = ctrl),
      error = \(e) NULL
    )
  }) %>% set_names(exp_species)
  
  # -------------------------
  # Coefficients plots (kept close to your original)
  # -------------------------
  lin_coef_df <- map_dfr(linear_species, \(sp) {
    mod <- models_linear[[sp]]
    summ <- summary(mod)$coefficients
    tibble(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"]
    ) %>%
      mutate(Coefficient = dplyr::recode(Coefficient, "(Intercept)" = "a", "x" = "b"))
  }) %>%
    mutate(
      species = factor(species, levels = linear_species),
      label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue))
    )
  
  p_coeff_linear <- ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = label, y = Value/2), vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models",
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "", y = "coefficient value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  exp_coef_df <- map_dfr(exp_species, \(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(tibble())
    summ <- summary(mod)$coefficients
    tibble(
      species = sp,
      Coefficient = tolower(rownames(summ)),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"]
    )
  }) %>%
    mutate(
      species = factor(species, levels = exp_species),
      label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue))
    )
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = label, y = Value/2), vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models",
         subtitle = expression(NDVI == a + b * e^{c * italic(x)} + epsilon),
         x = "", y = "coefficient value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0.5))
  
  combined_coeff <- p_coeff_linear + p_coeff_exp +
    plot_layout(widths = c(0.6, 1.4), guides = "collect")
  
  ggsave(filename = combined_coef_fig, plot = combined_coeff, device = "png",
         width = 14, height = 8, dpi = 300)
  
  # -------------------------
  # Panel (a): predictions + DELTA METHOD confidence bands
  # -------------------------
  # linear prediction CI
  pred_linear_ci <- map_dfr(linear_species, \(sp) {
    sp_data <- data_b %>% filter(species == sp)
    fit <- models_linear[[sp]]
    
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 200)
    nd <- tibble(x = x_seq)
    
    pr <- predict(fit, newdata = nd, interval = "confidence", level = conf_level)
    
    tibble(
      species = sp,
      x = x_seq,
      pred = pr[, "fit"],
      lower = pr[, "lwr"],
      upper = pr[, "upr"]
    )
  })
  
  # exponential prediction CI via delta method
  pred_exp_ci <- map_dfr(exp_species, \(sp) {
    sp_data <- data_b %>% filter(species == sp)
    fit <- models_exp[[sp]]
    if (is.null(fit)) return(tibble())
    
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 200)
    ci_tbl <- pred_exp_delta(fit, x_seq, conf_level = conf_level)
    
    ci_tbl %>% mutate(species = sp)
  })
  
  pred_all_ci <- bind_rows(pred_linear_ci, pred_exp_ci) %>%
    mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))
  
  p_combined <- ggplot() +
    # ✅ confidence band (behind)
    geom_ribbon(
      data = pred_all_ci,
      aes(x = x, ymin = lower, ymax = upper, fill = species),
      alpha = 0.18,
      color = NA
    ) +
    geom_point(
      data = data_b,
      aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
      alpha = 0.7
    ) +
    geom_line(
      data = pred_all_ci,
      aes(x = x, y = pred, color = species),
      linewidth = 1
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text",
             x = min(data_b$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -0.1, vjust = -0.3,
             fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "Pixels per bin (%)", range = c(1, 8),
                          labels = scales::percent_format(accuracy = 1)) +
    guides(color = guide_legend(order = 1),
           size  = guide_legend(order = 2)) +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 18, hjust = 0.1, vjust = 1),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title  = element_text(face = "bold", size = 16),
      axis.text   = element_text(color = "black", size = 14),
      legend.position = "bottom",
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  # -------------------------
  # Panels (b) + (c): x50 & slope (kept; uses your deltaMethod)
  # -------------------------
  lin_stats <- map_dfr(linear_species, \(sp) {
    mod <- models_linear[[sp]]
    summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]
    slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]
    p_val <- summ["x", "Pr(>|t|)"]
    x50 <- (threshold - intercept) / slope
    r2 <- summary(mod)$r.squared
    
    tibble(species = sp, model = "linear", x50 = x50, abs_slope = abs(slope),
           r_squared = r2, se = se_slope, p_val = p_val)
  })
  
  exp_stats <- map_dfr(exp_species, \(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(tibble())
    
    coefs <- coef(mod)
    a_val <- coefs["a"]; b_val <- coefs["b"]; c_val <- coefs["c"]
    
    x50 <- ifelse((threshold - a_val) > 0 & b_val > 0,
                  log((threshold - a_val)/b_val) / c_val,
                  NA_real_)
    
    slope50 <- c_val * (threshold - a_val)
    abs_slope <- abs(slope50)
    
    sp_data <- data_b %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = sp_data)
    r2 <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) /
      sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    
    dm <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                      parameterNames = c("a", "b", "c"))
    
    se <- dm$SE
    df_exp <- summary(mod)$df[2]
    tval <- slope50 / se
    p_val <- 2 * (1 - pt(abs(tval), df_exp))
    
    tibble(species = sp, model = "exponential", x50 = x50, abs_slope = abs_slope,
           r_squared = r2, se = se, p_val = p_val)
  })
  
  stats_all <- bind_rows(lin_stats, exp_stats) %>%
    mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))
  
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 18, hjust = 0.25, vjust = 1),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title  = element_text(face = "bold", size = 16),
      axis.text   = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2),
              color = "black", size = 5) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 18, hjust = 0.25, vjust = 1),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title  = element_text(face = "bold", size = 16),
      axis.text   = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  final_plot <- (p_combined + (p_x50 / p_slope)) +
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
  
  # plateau point (same)
  plateau_df <- map_dfr(exp_species, \(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(tibble())
    c_val <- coef(mod)["c"]
    tibble(species = sp, x_plateau = -log(0.05) / c_val)
  })
  
  print("Estimated soil water potential (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

# run
plot_NDVI_PSI_exp_linear_slope_coeff(
  data,
  "results_rootzone/test/NDVI_Q_PSIbin_exp_linear_coeff.png",
  "results_rootzone/test/NDVI_Q_PSIbin_exp_linear_slope.png",
  "results_rootzone/test/NDVI_Q_PSIbin_exp_linear_aic.png",
  bin_width = 50,
  min_count = 1000,
  conf_level = 0.95
)

##################### Germany and species-located data #################
library(tidyverse)
# -------------------------
# 1) BINNING: mean PSI per bin + mean NDVI per bin
#    (bins are consistent globally using global psi_min/max)
# -------------------------

NDVI_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  species_totals <- df %>%
    dplyr::group_by(area, species) %>%
    dplyr::summarise(total_pixels = dplyr::n(), .groups = "drop")
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df %>%
    dplyr::mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    dplyr::group_by(area, species, PSI_bin) %>%
    dplyr::summarise(
      bin_mean  = mean(soil_water_potential, na.rm = TRUE),
      avg_value = mean(NDVI, na.rm = TRUE),
      count     = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::left_join(species_totals, by = c("area", "species")) %>%
    dplyr::mutate(percentage = count / total_pixels) %>%
    dplyr::filter(count >= min_count) %>%
    dplyr::select(area, species, PSI_bin, bin_mean, avg_value, count, total_pixels, percentage)
}

### loess ###
plot_NDVI_Q_PSIbin_area <- function(
    df,
    out_file,
    bin_width = 50,
    min_count = 1000,
    width = 12,
    height = 5,
    dpi = 300
) {
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  facet_labs <- c(
    "Germany" = "(a) Germany",
    "species" = "(b) species"
  )
  
  p <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(area = factor(area, levels = c("Germany", "species"))) %>%
    ggplot(aes(bin_mean, avg_value)) +
    geom_point(aes(color = species, shape = species, size = percentage), alpha = 0.8) +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    geom_smooth(
      aes(color = species, fill = species),
      method = "loess",
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    facet_wrap(~ area, nrow = 1, labeller = as_labeller(facet_labs)) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
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
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
}

plot_NDVI_Q_PSIbin_area(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_Q_PSIbin_area_loess.png"
)

### gam k=6 bs=cs (area as panel) ###
plot_NDVI_Q_PSIbin_area <- function(
    df,
    out_file,
    bin_width = 50,
    min_count = 1000,
    k = 6,                 # GAM flexibility (increase for more wiggle)
    bs = "cs",             # shrinkage spline ("cs" is a good default)
    width = 12,
    height = 5,
    dpi = 300
) {
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  facet_labs <- c(
    "Germany" = "(a) Germany",
    "species" = "(b) species"
  )
  
  p <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(area = factor(area, levels = c("Germany", "species"))) %>%
    ggplot(aes(bin_mean, avg_value)) +
    geom_point(aes(color = species, shape = species, size = percentage), alpha = 0.8) +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    geom_smooth(
      aes(color = species, fill = species),
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    facet_wrap(~ area, nrow = 1, labeller = as_labeller(facet_labs)) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
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
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
  invisible(p)
}

plot_NDVI_Q_PSIbin_area(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_Q_PSIbin_area_gam.png",
  k = 5
)

### gam k=6 bs=cs (species as panel) ###
plot_NDVI_Q_PSIbin_species4_area2 <- function(
    df,
    out_file,
    bin_width = 50,
    min_count = 1000,
    k = 6,
    bs = "cs",
    width = 12,
    height = 8,
    dpi = 300
) {
  
  # species order + colors for the two "area" curves
  sp_levels <- c("Oak", "Beech", "Spruce", "Pine")
  area_cols <- c("Germany" = "blue", "species" = "red")
  
  # optional nicer facet text; keep simple labels
  facet_labs <- c("Germany" = "Germany", "species" = "species")
  
  p <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(
      area    = factor(area, levels = c("Germany", "species")),
      species = factor(species, levels = sp_levels)
    ) %>%
    ggplot(aes(bin_mean, avg_value)) +
    
    # points (optional): colored by area, shaped by area
    geom_point(
      aes(color = area, shape = area, size = percentage),
      alpha = 0.8
    ) +
    
    # GAM smooths (two lines per facet), ribbon follows area color
    geom_smooth(
      aes(color = area, fill = area),
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.20,
      linewidth = 1
    ) +
    
    # 4 panels, ordered Oak -> Beech -> Spruce -> Pine
    facet_wrap(~ species, nrow = 2) +
    
    # blue/red for areas
    scale_color_manual(values = area_cols, name = "") +
    scale_fill_manual(values = area_cols, guide = "none") +
    
    # shapes for the two lines/points
    scale_shape_manual(values = c("Germany" = 16, "species" = 17), name = "") +
    
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
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
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
  invisible(p)
}

plot_NDVI_Q_PSIbin_species4_area2(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_Q_PSIbin_4species_2areas_gam.png",
  k = 5
)

# ---------- 1) NDVI difference along PSI bins ----------
NDVI_diff_species_minus_Germany_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df_binned <- df %>%
    dplyr::mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  sp_tbl <- df_binned %>%
    dplyr::filter(area == "species") %>%
    dplyr::group_by(species, PSI_bin) %>%
    dplyr::summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      sp_mean  = mean(NDVI, na.rm = TRUE),
      sp_n     = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::filter(sp_n >= min_count)
  
  de_tbl <- df_binned %>%
    dplyr::filter(area == "Germany") %>%
    dplyr::group_by(species, PSI_bin) %>%
    dplyr::summarise(
      de_mean = mean(NDVI, na.rm = TRUE),
      de_n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(de_n >= min_count)
  
  sp_tbl %>%
    dplyr::inner_join(de_tbl, by = c("species", "PSI_bin")) %>%
    dplyr::mutate(diff_NDVI = sp_mean - de_mean) %>%
    dplyr::select(species, PSI_bin, bin_mean, sp_mean, de_mean, sp_n, de_n, diff_NDVI)
}


plot_NDVI_diff_PSI <- function(data, figure_output = NULL,
                               bin_width = 50, min_count = 1000,
                               k = 3, bs = "cs",
                               width = 8, height = 6, dpi = 300) {
  
  df_diff <- NDVI_diff_species_minus_Germany_PSIbin(data, bin_width, min_count) %>%
    dplyr::mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  p <- ggplot2::ggplot(df_diff,
                       ggplot2::aes(x = bin_mean, y = diff_NDVI, color = species)) +
    
    ggplot2::geom_hline(yintercept = 0,
                        linewidth = 0.6,
                        linetype = "dashed",
                        color = "grey40") +
    
    ggplot2::geom_point(alpha = 0.5, size = 2) +
    
    ggplot2::geom_smooth(
      ggplot2::aes(fill = species),
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    
    ggplot2::scale_color_manual(values = cb_palette, name = "") +
    ggplot2::scale_fill_manual(values = cb_palette, guide = "none") +
    
    ggplot2::labs(
      x = "soil water potential (kPa)",
      y = expression(Delta[NDVI]~"(species - Germany)")
    ) +
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title  = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  print(p)
  
  if (!is.null(figure_output)) {
    dir.create(dirname(figure_output), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(figure_output, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  invisible(p)
}

plot_NDVI_diff_PSI(
  data = DT_combined,
  figure_output = "results_rootzone/test/NDVI_PSIbin_area_diff_till2022.png",
  bin_width = 50,
  min_count = 1000,
  k = 3
)

### difference bar plot 
plot_NDVI_diff_PSI_bar <- function(data, figure_output = NULL,
                                   bin_width = 50, min_count = 1000,
                                   weighted = TRUE,
                                   width = 6, height = 5, dpi = 300) {
  
  library(dplyr)
  library(ggplot2)
  
  df_bins <- NDVI_diff_species_minus_Germany_PSIbin(data, bin_width, min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      w = pmin(sp_n, de_n)
    )
  
  df_bar <- df_bins %>%
    group_by(species) %>%
    summarise(
      mean_diff = if (weighted) weighted.mean(diff_NDVI, w = w, na.rm = TRUE)
      else mean(diff_NDVI, na.rm = TRUE),
      
      # SE across bins (optionally weighted)
      se_diff   = if (weighted) {
        # effective sample size
        w2 <- w / sum(w)
        mu <- sum(w2 * diff_NDVI, na.rm = TRUE)
        var_w <- sum(w2 * (diff_NDVI - mu)^2, na.rm = TRUE) / (1 - sum(w2^2))
        sqrt(var_w) / sqrt(n())
      } else {
        sd(diff_NDVI, na.rm = TRUE) / sqrt(n())
      },
      .groups = "drop"
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  p <- ggplot(df_bar, aes(x = species, y = mean_diff, fill = species)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "grey40") +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
                  width = 0.2, linewidth = 0.6) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = expression(mean~Delta[NDVI]~"(species - Germany)")) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  if (!is.null(figure_output)) {
    dir.create(dirname(figure_output), recursive = TRUE, showWarnings = FALSE)
    ggsave(figure_output, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  invisible(p)
}

plot_NDVI_diff_PSI_bar(
  DT_combined,
  figure_output = "results_rootzone/test/NDVI_PSIbin_area_diff_bar_noW_till2022.png",
  bin_width = 50,
  min_count = 1000,
  weighted = FALSE
)

plot_NDVI_diff_PSI_bar_both <- function(
    data,
    figure_output = NULL,
    bin_width = 50,
    min_count = 1000,
    width = 9,
    height = 5,
    dpi = 300
) {
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # ---- bin-level differences ----
  df_bins <- NDVI_diff_species_minus_Germany_PSIbin(data, bin_width, min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      w = pmin(sp_n, de_n)
    )
  
  # ---- summaries ----
  df_weighted <- df_bins %>%
    group_by(species) %>%
    summarise(
      mean_diff = weighted.mean(diff_NDVI, w = w, na.rm = TRUE),
      se_diff   = sd(diff_NDVI, na.rm = TRUE) / sqrt(n()),
      type = "(a) weighted",
      .groups = "drop"
    )
  
  df_unweighted <- df_bins %>%
    group_by(species) %>%
    summarise(
      mean_diff = mean(diff_NDVI, na.rm = TRUE),
      se_diff   = sd(diff_NDVI, na.rm = TRUE) / sqrt(n()),
      type = "(b) non-weighted",
      .groups = "drop"
    )
  
  df_bar <- bind_rows(df_weighted, df_unweighted)
  
  # ---- colors ----
  cb_palette <- c(
    "Oak"   = "#E69F00",
    "Beech" = "#0072B2",
    "Spruce"= "#009E73",
    "Pine"  = "#F0E442"
  )
  
  # ---- plot ----
  p <- ggplot(df_bar, aes(x = species, y = mean_diff, fill = species)) +
    
    geom_hline(
      yintercept = 0,
      linewidth = 0.6,
      linetype = "dashed",
      color = "grey40"
    ) +
    
    geom_col(width = 0.7, alpha = 0.9) +
    
    geom_errorbar(
      aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
      width = 0.2, linewidth = 0.6
    ) +
    
    facet_wrap(~ type, nrow = 1) +
    
    scale_fill_manual(values = cb_palette, guide = "none") +
    
    labs(
      x = "",
      y = expression(mean~Delta[NDVI]~"(species - Germany)")
    ) +
    
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  if (!is.null(figure_output)) {
    dir.create(dirname(figure_output), recursive = TRUE, showWarnings = FALSE)
    ggsave(
      figure_output,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  invisible(p)
}

plot_NDVI_diff_PSI_bar_both(
  DT_combined,
  figure_output = "results_rootzone/test/NDVI_PSIbin_area_diff_bar_weighted_vs_unweighted.png",
  bin_width = 50,
  min_count = 1000
)


### combined Germany vs Species-locations and Difference ###
library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)

# ----- shared theme -----
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

# palettes
sp_levels  <- c("Oak", "Beech", "Spruce", "Pine")
sp_cols    <- c("Oak"="#E69F00", "Beech"="#0072B2", "Spruce"="#009E73", "Pine"="#F0E442")
de_col     <- "grey40"  # Germany curve color (change if you want)

NDVI_diff_species_minus_Germany_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df_binned <- df %>%
    dplyr::mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # totals per species in each area (for percentages)
  totals <- df_binned %>%
    dplyr::group_by(area, species) %>%
    dplyr::summarise(total_n = dplyr::n(), .groups = "drop")
  
  sp_tbl <- df_binned %>%
    dplyr::filter(area == "species") %>%
    dplyr::group_by(species, PSI_bin) %>%
    dplyr::summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      sp_mean  = mean(NDVI, na.rm = TRUE),
      sp_n     = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::filter(sp_n >= min_count) %>%
    dplyr::left_join(
      totals %>% dplyr::filter(area == "species") %>% dplyr::select(species, sp_total = total_n),
      by = "species"
    ) %>%
    dplyr::mutate(sp_pct = sp_n / sp_total)
  
  de_tbl <- df_binned %>%
    dplyr::filter(area == "Germany") %>%
    dplyr::group_by(species, PSI_bin) %>%
    dplyr::summarise(
      de_mean = mean(NDVI, na.rm = TRUE),
      de_n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(de_n >= min_count) %>%
    dplyr::left_join(
      totals %>% dplyr::filter(area == "Germany") %>% dplyr::select(species, de_total = total_n),
      by = "species"
    ) %>%
    dplyr::mutate(de_pct = de_n / de_total)
  
  sp_tbl %>%
    dplyr::inner_join(de_tbl, by = c("species", "PSI_bin")) %>%
    dplyr::mutate(
      diff_NDVI = sp_mean - de_mean,
      pct_bin  = pmin(sp_pct, de_pct)   # or use (sp_pct + de_pct)/2
    ) %>%
    dplyr::select(species, PSI_bin, bin_mean, sp_mean, de_mean, sp_n, de_n, sp_pct, de_pct, pct_bin, diff_NDVI)
}

# ---------- (a) left plot ----------
plot_NDVI_Q_PSIbin_species4_area2_g <- function(
    df,
    bin_width = 50,
    min_count = 1000,
    k = 6,
    bs = "cs",
    show_germany = TRUE   # set TRUE if you want a grey reference line
) {
  
  d_all <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    dplyr::mutate(species = factor(species, levels = sp_levels))
  
  d_sp <- d_all %>%
    dplyr::filter(area == "species")
  
  d_de <- d_all %>%
    dplyr::filter(area == "Germany")
  
  ggplot(d_sp, aes(x = bin_mean, y = avg_value, colour = species, fill = species)) +
    
    # optional Germany reference (grey; NOT in legend)
    {if (show_germany)
      geom_smooth(
        data = d_de,
        aes(x = bin_mean, y = avg_value, group = species),
        method = "gam",
        formula = y ~ s(x, k = k, bs = bs),
        se = FALSE,
        colour = "grey50",
        linewidth = 0.8,
        inherit.aes = FALSE
      )
    } +
    
    geom_point(aes(size = percentage), alpha = 0.8) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.20,
      linewidth = 1
    ) +
    facet_wrap(~ species, nrow = 2) +
    scale_colour_manual(values = sp_cols, name = "") +
    scale_fill_manual(values = sp_cols, name = "") +
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
    ) +
    # line + ribbon in legend key
    guides(
      colour = guide_legend(
        override.aes = list(alpha = 0.20, linewidth = 1.2)
      ),
      fill = "none"
    ) +
    base_theme
}

# ---------- (b) right plot ----------
plot_NDVI_diff_PSI_g <- function(
    data,
    bin_width = 50,
    min_count = 1000,
    k = 3,
    bs = "cs"
) {
  
  df_diff <- NDVI_diff_species_minus_Germany_PSIbin(data, bin_width, min_count) %>%
    dplyr::mutate(species = factor(species, levels = sp_levels))
  
  ggplot(df_diff, aes(x = bin_mean, y = diff_NDVI, colour = species, fill = species)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_point(aes(size = pct_bin), alpha = 0.6) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    scale_colour_manual(values = sp_cols, name = "") +
    scale_fill_manual(values = sp_cols, name = "") +
    scale_shape_manual(values = c(Oak = 16, Beech = 17, Spruce = 15, Pine = 18), name ="") +
    scale_size_continuous(
      name   = "bin % (min area support)",
      range  = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = expression(Delta[NDVI]~"(species - Germany)")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(alpha = 0.25, linewidth = 1.2)
      ),
      fill = "none"
    ) +
    base_theme +
    theme(legend.position = "none")
}

# ---------- combine (a) + (b) ----------
plot_combined_big <- function(
    df,
    out_file,
    bin_width = 50,
    min_count = 1000,
    k_left = 5,
    k_right = 3,
    bs = "cs",
    width = 16,
    height = 8,
    dpi = 300
) {
  p_left  <- plot_NDVI_Q_PSIbin_species4_area2_g(
    df, bin_width, min_count, k = k_left, bs = bs, show_germany = FALSE
  )
  p_right <- plot_NDVI_diff_PSI_g(df, bin_width, min_count, k = k_right, bs = bs)
  
  p <- (p_left | p_right) +
    patchwork::plot_layout(widths = c(1.25, 1), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom")
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
  invisible(p)
}

plot_combined_big(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_combined_a_b.png",
  bin_width = 50,
  min_count = 1000,
  k_left = 5,
  k_right = 3
)



### combined Germany vs Species-locations and Bar difference ###
# change ONLY the LEFT plot mapping:
# - species points/lines colored by species (sp_cols)
# - Germany shown in grey (single colour, not in legend)
# - keep size legend for percentage; keep species colour legend like your exp_linear plot

plot_NDVI_left4_right2_combined <- function(
    df,
    out_file = "results_rootzone/test/NDVI_left4_right2_combined.png",
    bin_width = 50,
    min_count = 1000,
    k_left = 6,
    bs_left = "cs",
    width = 16,
    height = 9,
    dpi = 300,
    sp_levels = c("Oak", "Beech", "Spruce", "Pine"),
    sp_cols   = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    sp_shapes = c(Oak=16, Beech=17, Spruce=15, Pine=18),
    germany_grey = "grey50"
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(mgcv)
    library(scales)
    library(patchwork)
  })
  
  # ---- helper: NDVI diff species - Germany by PSI bin ----
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
  
  # -------------------------
  # (LEFT) 4 panels:
  #   - species coloured by species (sp_cols) + shapes by species
  #   - Germany as grey reference line (not in legend)
  # -------------------------
  d_all <- NDVI_PSIbin(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  d_sp <- d_all %>% filter(area == "species")
  d_de <- d_all %>% filter(area == "Germany")
  
  p_left <- ggplot(d_sp, aes(bin_mean, avg_value)) +
    
    # Germany reference (grey) per facet/species, NOT in legend
    geom_point(
      data = d_de,
      aes(x = bin_mean, y = avg_value, size = percentage),
      colour = germany_grey,
      alpha = 0.35,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_smooth(
      data = d_de,
      aes(x = bin_mean, y = avg_value, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_left, bs = bs_left),
      se = TRUE,
      alpha = 0.18,
      colour = NA,          # no border line from this layer
      fill   = "grey70",
      linewidth = 0,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_smooth(
      data = d_de,
      aes(x = bin_mean, y = avg_value, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_left, bs = bs_left),
      se = FALSE,
      colour = germany_grey,
      linewidth = 1,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    # species smooth line (this is what appears in legend)
    geom_smooth(
      aes(colour = species, group = species),
      method  = "gam",
      formula = y ~ s(x, k = k_left, bs = bs_left),
      se = FALSE,
      linewidth = 1.2
    ) +
    
    facet_wrap(~ species, nrow = 2) +
    scale_colour_manual(values = sp_cols, name = "") +
    scale_fill_manual(values = sp_cols,  name = "") +
    scale_shape_manual(values = sp_shapes, name = "") +
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
      fill = "none"
    ) +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1),
      legend.position = "bottom"
    )
  
  # -------------------------
  # (RIGHT) bars: weighted vs unweighted (unchanged)
  # -------------------------
  df_diff <- NDVI_diff_species_minus_Germany_PSIbin(df, bin_width, min_count) %>%
    mutate(species = factor(species, levels = sp_levels))
  
  df_bar <- df_diff %>%
    group_by(species) %>%
    summarise(
      mean_unweighted = mean(diff_NDVI, na.rm = TRUE),
      mean_weighted   = weighted.mean(diff_NDVI, w = pct_bin, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean_weighted, mean_unweighted),
                 names_to = "type", values_to = "value") %>%
    mutate(
      type = recode(type,
                    mean_weighted   = "Weighted (by min area support)",
                    mean_unweighted = "Unweighted")
    )
  
  p_right_top <- ggplot(df_bar %>% filter(type == "Weighted (by min area support)"),
                        aes(x = species, y = value, fill = species)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) +
    scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = expression(Delta[NDVI]~"(species - Germany)"), title = "(b) Weighted") +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  p_right_bottom <- ggplot(df_bar %>% filter(type == "Unweighted"),
                           aes(x = species, y = value, fill = species)) +
    geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.7) +
    scale_fill_manual(values = sp_cols, guide = "none") +
    labs(x = "", y = expression(Delta[NDVI]~"(species - Germany)"), title = "(c) Unweighted") +
    base_theme +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  p_right <- p_right_top / p_right_bottom
  
  # -------------------------
  # combine
  # -------------------------
  p <- (p_left | p_right) +
    plot_layout(widths = c(1.6, 1), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  print(p)
  invisible(p)
}

plot_NDVI_left4_right2_combined(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_left4_right2_combined.png",
  bin_width = 50,
  min_count = 1000
)
