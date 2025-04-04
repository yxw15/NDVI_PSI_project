# =============================================================================
# Organized NDVI_PSI_project Script
# =============================================================================

# Set working directory and load data
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# =============================================================================
# Load Libraries
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(nlme)
library(car)
library(purrr)
library(tibble)

# =============================================================================
# Global Custom Theme for Consistent Figure Style
# =============================================================================

my_theme <- theme_minimal() +
  theme(
    legend.position = "top",    # Place legend at the top
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
    plot.caption = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14, color = "black"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


# =============================================================================
# Data Summarizing Functions
# =============================================================================

df_average_year <- function(df_all) {
  # Determine correct column for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) "Quantiles" else if ("Proportions" %in% names(df_all)) "Proportions" else stop("Neither Quantiles nor Proportions column found.")
  
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  df_avg <- df_all %>%
    group_by(year) %>%
    summarise(
      !!avg_column_name := mean(.data[[value_column]], na.rm = TRUE),
      avg_psi = mean(soil_water_potential, na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_avg)
}

df_average_year_species <- function(df_all) {
  value_column <- if ("Quantiles" %in% names(df_all)) "Quantiles" else if ("Proportions" %in% names(df_all)) "Proportions" else stop("Neither Quantiles nor Proportions column found.")
  
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  df_avg <- df_all %>%
    group_by(year, species) %>%
    summarise(
      !!avg_column_name := mean(.data[[value_column]], na.rm = TRUE),
      avg_psi = mean(soil_water_potential, na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_avg)
}

# =============================================================================
# Binning Functions
# =============================================================================

NDVI_PSIbin <- function(df, bin_width = 50) {
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  species_totals <- df %>% group_by(species) %>% summarise(total_pixels = n(), .groups = "drop")
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.001) %>%
    select(species, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_PSIbin_species)
}

NDVI_TDiffbin <- function(df, bin_width = 3) {
  value_col <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  species_totals <- df %>% group_by(species) %>% summarise(total_pixels = n(), .groups = "drop")
  
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df <- df %>%
    mutate(TDiff_bin = cut(transpiration_deficit, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  get_bin_median <- function(bin_label) {
    nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
    mean(nums)
  }
  
  meanNDVI_TDiffbin_species <- df %>%
    group_by(species, TDiff_bin) %>%
    summarise(
      avg_value = mean(.data[[value_col]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(bin_median = sapply(as.character(TDiff_bin), get_bin_median)) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.0001) %>%
    select(species, TDiff_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_TDiffbin_species)
}

TDiff_PSIbin <- function(df, bin_width = 50) {
  # Here the value is transpiration_deficit
  value_column <- "transpiration_deficit"
  
  species_totals <- df %>% group_by(species) %>% summarise(total_pixels = n(), .groups = "drop")
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  meanTDiff_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_transpiration_deficit = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage > 0.001) %>%
    select(species, PSI_bin, bin_median, avg_transpiration_deficit, count, total_pixels, percentage)
  
  return(meanTDiff_PSIbin_species)
}

# =============================================================================
# Plotting Functions: Time Series and Correlation
# =============================================================================

plot_time_series_and_correlation_combined <- function(df_all, output_path) {
  library(patchwork)
  library(cowplot)
  
  df_species <- df_average_year_species(df_all)
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) "avg_quantile" else if ("avg_proportion" %in% names(df_species)) "avg_proportion" else stop("No NDVI column found.")
  ndvi_label <- if (ndvi_column == "avg_quantile") "NDVI quantiles" else "NDVI proportions"
  x_scale <- scale_x_continuous(breaks = seq(2003, 2024, by = 2), limits = c(2003, 2024))
  
  # Time Series Plots
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(y = "NDVI quantiles", color = "Species", caption = "(a)") +
    scale_color_manual(values = cb_palette) +
    my_theme +
    theme(legend.position = "bottom") +
    x_scale
  
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(y = "Soil water potential", caption = "(b)") +
    scale_color_manual(values = cb_palette) +
    my_theme +
    theme(legend.position = "none") +
    x_scale
  
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(y = "Transpiration deficit", caption = "(c)") +
    scale_color_manual(values = cb_palette) +
    my_theme +
    theme(legend.position = "none") +
    x_scale
  
  ts_plot <- p1 / p2 / p3
  
  # Correlation Plots
  annotations_d <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_psi, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  annotations_e <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_tdiff, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  text_d <- paste0(species_order, ": r = ", annotations_d$corr[match(species_order, annotations_d$species)], collapse = "\n")
  text_e <- paste0(species_order, ": r = ", annotations_e$corr[match(species_order, annotations_e$species)], collapse = "\n")
  
  p_d <- ggplot(df_species, aes(x = avg_psi, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Soil water potential", y = ndvi_label, caption = "(d)") +
    scale_color_manual(values = cb_palette) +
    my_theme +
    theme(legend.position = "none") +
    annotate("text", x = -Inf, y = Inf, label = text_d, hjust = 0, vjust = 1, size = 5, color = "black")
  
  p_e <- ggplot(df_species, aes(x = avg_tdiff, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Transpiration deficit", y = "", caption = "(e)") +
    scale_color_manual(values = cb_palette) +
    my_theme +
    theme(legend.position = "none") +
    annotate("text", x = Inf, y = Inf, label = text_e, hjust = 1, vjust = 1, size = 5, color = "black")
  
  corr_plot <- p_d | p_e
  combined_main <- ts_plot / corr_plot & theme(legend.position = "top")
  
  legend <- suppressWarnings(cowplot::get_legend(p1))
  final_plot <- cowplot::plot_grid(combined_main, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  print(final_plot)
  ggsave(filename = output_path, plot = final_plot, width = 12, height = 14, dpi = 300)
}

# =============================================================================
# Plotting Functions: Yearly Mean Linear Coefficients
# =============================================================================

plot_yearly_mean_linear_coeffs <- function(df_all, output_path) {
  df_species <- df_average_year_species(df_all)
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) "avg_quantile" else if ("avg_proportion" %in% names(df_species)) "avg_proportion" else stop("No NDVI column found.")
  ndvi_label <- if (ndvi_column == "avg_quantile") "NDVI quantiles" else "NDVI proportions"
  
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Regression Coefficients for NDVI ~ Soil Water Potential
  coeffs_psi <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_psi")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, term = names(coef_vals), estimate = coef_vals, p.value = p_vals, stringsAsFactors = FALSE)
  }))
  coeffs_psi$term <- dplyr::recode(coeffs_psi$term, "(Intercept)" = "a", "avg_psi" = "b")
  coeffs_psi$species <- factor(coeffs_psi$species, levels = species_order)
  coeffs_psi$label <- ifelse(coeffs_psi$p.value < 0.05, "*", sprintf("%.2f", coeffs_psi$p.value))
  
  p_coeff_psi <- ggplot(coeffs_psi, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2), position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI ~ Soil Water Potential", y = "Coefficient Value", caption = "(a)") +
    my_theme
  
  # Regression Coefficients for NDVI ~ Transpiration Deficit
  coeffs_tdiff <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_tdiff")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, term = names(coef_vals), estimate = coef_vals, p.value = p_vals, stringsAsFactors = FALSE)
  }))
  coeffs_tdiff$term <- dplyr::recode(coeffs_tdiff$term, "(Intercept)" = "a", "avg_tdiff" = "b")
  coeffs_tdiff$species <- factor(coeffs_tdiff$species, levels = species_order)
  coeffs_tdiff$label <- ifelse(coeffs_tdiff$p.value < 0.05, "*", sprintf("%.2f", coeffs_tdiff$p.value))
  
  p_coeff_tdiff <- ggplot(coeffs_tdiff, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2), position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI ~ Transpiration Deficit", caption = "(b)") +
    my_theme
  
  combined_plot <- p_coeff_psi + p_coeff_tdiff + plot_layout(guides = "collect") + theme(legend.position = "top")
  
  print(combined_plot)
  ggsave(filename = output_path, plot = combined_plot, width = 14, height = 8, dpi = 300)
}

# =============================================================================
# Plotting Functions: Combined AIC and R2 Comparisons
# =============================================================================

plot_combined_AIC_R2 <- function(data, save_combined_fig) {
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Panel A: NDVI ~ PSIbin
  data_a <- NDVI_PSIbin(data)
  data_a <- na.omit(data_a)
  data_a$species <- factor(data_a$species, levels = species_order)
  data_a <- data_a %>% mutate(x = -bin_median)
  
  start_list_a <- list(a = 5, b = 3, c = 0.001)
  control_params_a <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  aic_results_a <- list()
  for (sp in levels(data_a$species)) {
    sp_data <- data_a %>% filter(species == sp)
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    aic_exp <- NA; r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data, start = start_list_a, control = control_params_a)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data$avg_value - mean(sp_data$avg_value))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(species = sp, Model = c("Linear", "Exponential"),
                             AIC = c(aic_linear, aic_exp), R2 = c(r2_linear, r2_exp))
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_a[[sp]] <- sp_results
  }
  aic_df_a <- do.call(rbind, aic_results_a)
  aic_df_a$species <- factor(aic_df_a$species, levels = species_order)
  
  model_palette_shared <- c("Linear" = "orange", "Exponential" = "dodgerblue")
  
  p_a <- ggplot(aic_df_a, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI ~ PSIbin", caption = "(a)") +
    scale_fill_manual(values = model_palette_shared) +
    my_theme
  
  # Panel B: NDVI ~ TDiffbin
  data_b <- NDVI_TDiffbin(data)
  data_b <- na.omit(data_b)
  data_b$species <- factor(data_b$species, levels = species_order)
  data_b <- data_b %>% mutate(x = bin_median)
  
  start_list_b <- list(a = 5, b = 7, c = 0.04)
  control_params_b <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  aic_results_b <- list()
  for (sp in levels(data_b$species)) {
    sp_data <- data_b %>% filter(species == sp)
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    aic_exp <- NA; r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data, start = start_list_b, control = control_params_b)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data$avg_value - mean(sp_data$avg_value))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(species = sp, Model = c("Linear", "Exponential"),
                             AIC = c(aic_linear, aic_exp), R2 = c(r2_linear, r2_exp))
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_b[[sp]] <- sp_results
  }
  aic_df_b <- do.call(rbind, aic_results_b)
  aic_df_b$species <- factor(aic_df_b$species, levels = species_order)
  
  p_b <- ggplot(aic_df_b, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI ~ TDiffbin", caption = "(b)") +
    scale_fill_manual(values = model_palette_shared) +
    my_theme
  
  # Panel C: TDiff ~ PSIbin (with polynomial fits)
  data_c <- TDiff_PSIbin(data)
  data_c <- na.omit(data_c)
  data_c$species <- factor(data_c$species, levels = species_order)
  data_c <- data_c %>% mutate(x = bin_median)
  
  aic_results_c <- list()
  for (sp in levels(data_c$species)) {
    sp_data <- data_c %>% filter(species == sp)
    lm_linear <- lm(avg_transpiration_deficit ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    lm_poly2 <- lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    r2_poly2 <- summary(lm_poly2)$r.squared
    
    lm_poly3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    r2_poly3 <- summary(lm_poly3)$r.squared
    
    sp_results <- data.frame(species = sp, Model = c("Linear", "Poly2", "Poly3"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3),
                             R2 = c(r2_linear, r2_poly2, r2_poly3))
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_c[[sp]] <- sp_results
  }
  aic_df_c <- do.call(rbind, aic_results_c)
  aic_df_c$species <- factor(aic_df_c$species, levels = species_order)
  
  model_palette_c <- c("Linear" = "orange", "Poly2" = "dodgerblue", "Poly3" = "green4")
  
  p_c <- ggplot(aic_df_c, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "TDiff ~ PSIbin", caption = "(c)") +
    scale_fill_manual(values = model_palette_c) +
    my_theme
  
  combined_top <- (p_a | p_b) + plot_layout(guides = "collect") + theme(legend.position = "top")
  combined_plot <- combined_top / p_c
  
  dir.create(dirname(save_combined_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_combined_fig, plot = combined_plot, width = 14, height = 12, dpi = 300)
  print(combined_plot)
}

# =============================================================================
# Plotting Function: NDVI Q ~ PSIbin Exponential Model and Slope at Negative PSI
# =============================================================================

plot_NDVI_Q_PSIbin_exp_slope_negPSI <- function(data, save_coeff_fig, save_slope_fig) {
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  threshold <- 11.5
  
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        log((threshold - a)/b) / c,
                        NA),
           slope50 = c * (threshold - a))
  
  # Panel A: Observed Data with Fitted Curves
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  x_scale <- scale_x_continuous()
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x), y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles", caption = "(a)") +
    my_theme +
    coord_cartesian(clip = "off")
  
  # Panel B: Bar Plot of x50 Values
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(y = "Soil Water Potential", caption = "(b)") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    my_theme +
    theme(legend.position = "none")
  
  # Panel C: Absolute Slope at x50 with Error Bars and Annotations
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    t_val <- slope50_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) / sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    
    tibble(species = sp,
           slope50 = slope50_est,
           slope_abs = abs(slope50_est),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05, sprintf("%.2f*", r_squared), sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 5) +
    labs(y = "Absolute Slope", caption = "(c)") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    my_theme +
    theme(legend.position = "none")
  
  final_slope_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") +
    theme(legend.position = "top", legend.title = element_blank())
  
  print(final_slope_plot)
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  # Coefficient Plot
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  coeff_stats <- coeff_stats %>%
    mutate(label = ifelse(`Pr(>|t|)` < 0.05, "*", sprintf("%.2f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"), names_to = "Coefficient", values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2), color = "black", size = 4, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(x = "", y = "Coefficient Value") +
    my_theme
  
  print(p_coeffs)
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, width = 10, height = 8, dpi = 300)
}

# =============================================================================
# Execute Plot Functions with Output Paths
# =============================================================================

plot_time_series_and_correlation_combined(all_results_df, "results/key_displays/time_series_Quantiles_PSI_TDiff_species.png")
plot_yearly_mean_linear_coeffs(all_results_df, "results/key_displays/yearly_mean_linear_coeffs.png")
plot_NDVI_Q_PSIbin_exp_slope_negPSI(all_results_df, 
                                    "results/key_displays/NDVI_Q_PSIbin_exp_coeff_negPSI.png",
                                    "results/key_displays/NDVI_Q_PSIbin_exp_slope_negPSI.png")
plot_combined_AIC_R2(all_results_df, "results/key_displays/AIC.png")
