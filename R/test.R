# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

source("R/functions_NDVI.R")
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

plot_box_TDiff_PSI(NDVI_PSI_TDiff_species, "results/Figures/mean_box_TDiff_PSI.png")
plot_box_merge_TDiff_PSI(NDVI_PSI_TDiff_species, "results/Figures/mean_box_merge_TDiff_PSI.png")
plot_mean_box_TDiff_PSI(NDVI_PSI_TDiff_species, "results/Figures/mean_TDiff_PSI.png")

plot_mean_box_NDVI_PSI_with_slope(NDVI_PSI_TDiff_species, "results/Figures/mean_box_NDVI_PSI_with_slope.png")
plot_mean_box_NDVI_PSI_with_slope_bin(NDVI_PSI_TDiff_species, "results/Figures/mean_box_NDVI_PSI_with_slope_bin.png")

plot_mean_box_TDiff_PSI_with_slope(NDVI_PSI_TDiff_species, "results/Figures/mean_box_TDiff_PSI_with_slope.png")
plot_mean_box_TDiff_PSI_with_slope_bin(NDVI_PSI_TDiff_species, "results/Figures/mean_box_TDiff_PSI_with_slope_bin.png")

plot_mean_box_TDiff_PSI_with_slope_filter(NDVI_PSI_TDiff_species, "results/Figures/mean_box_TDiff_PSI_with_slope_filter.png")
plot_mean_box_TDiff_PSI_with_slope_bin_filter(NDVI_PSI_TDiff_species, "results/Figures/mean_box_TDiff_PSI_with_slope_bin_filter.png")

plot_NDVI_PDM_PSIbin_slope(all_results_df, "results/Figures/Proportions_PSIbin_slope.png")

plot_box_NDVI_TDiff(NDVI_PSI_TDiff_species, "results/Figures/box_NDVI_TDiff.png")
plot_box_merge_NDVI_TDiff(NDVI_PSI_TDiff_species, "results/Figures/box_merge_NDVI_TDiff.png")
plot_mean_box_NDVI_TDiff(NDVI_PSI_TDiff_species, "results/Figures/box_mean_NDVI_TDiff.png")
plot_mean_box_NDVI_TDiff_with_slope(NDVI_PSI_TDiff_species, "results/Figures/box_mean_NDVI_TDiff_slope.png")
plot_mean_box_NDVI_TDiff_with_slope_bin(NDVI_PSI_TDiff_species, "results/Figures/box_mean_NDVI_TDiff_slope_bin.png")

plot_NDVI_PSIbin_slope(all_results_df, 
                       save_coeff_fig = "results/Figures/Quantiles_PSIbin_Coeffs.png", 
                       save_slope_fig = "results/Figures/Quantiles_PSIbin.png")

plot_NDVI_PSIbin_slope_each(all_results_df, 
                            save_coeff_fig = "results/Figures/Quantiles_PSIbin_Coeffs_each.png", 
                            save_slope_fig = "results/Figures/Quantiles_PSIbin_each.png")

plot_correlation_NDVI_TDiff_avg(all_results_df,
                                plot_corr_path = "results/Figures/Quantiles_TDiff_corr.png")


plot_series_NDVI_TDiff_same_month(all_results_df, 
                                  plot_series_path = "results/Figures/Quantiles_TDiff_series.png")


plot_series_NDVI_TDiff_same_month_species(all_results_df, 
                                          plot_species_path = "results/Figures/Proportions_TDiff_series_species.png")

plot_correlation_NDVI_TDiff_species_avg(all_results_df, 
                                        plot_corr_path = "results/Figures/Proportions_TDiff_corr_species.png")


plot_correlation_NDVI_PSI_avg(all_results_df, "results/Figures/meanQuantiles_PSIbin_correlation.png")
plot_correlation_NDVI_PSI_species_avg(all_results_df, "results/Figures/meanQuantiles_PSIbin_correlation_species.png")
plot_series_NDVI_PSI_same_month(all_results_df, "results/Figures/meanQuantiles_PSIbin_time_series.png")
plot_series_NDVI_PSI_same_month_species(all_results_df, "results/Figures/meanQuantiles_PSIbin_time_series_species.png")
plot_NDVI_PSIbin_slope(all_results_df, "results/Figures/meanQuantiles_PSIbin_non_linear.png")

# load Proportions
plot_correlation_NDVI_PSI_avg(all_results_df, "results/Figures/meanProportions_PSIbin_correlation.png")
plot_correlation_NDVI_PSI_species_avg(all_results_df, "results/Figures/meanProportions_PSIbin_correlation_species.png")
plot_series_NDVI_PSI_same_month(all_results_df, "results/Figures/meanProportions_PSIbin_time_series.png")
plot_series_NDVI_PSI_same_month_species(all_results_df, "results/Figures/meanProportions_PSIbin_time_series_species.png")
plot_NDVI_PDM_PSIbin_slope(all_results_df, "results/Figures/meanProportions_PSIbin_non_linear.png")


# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")


plot_NDVI_PSIbin_slope_each <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept in case of future needs)
  library(tidyr)    # for pivot_longer
  library(ggpattern)  # for patterned bars in Panel C
  
  # Identify the value column and order species; define the color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create a positive soil water potential (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum soil water potential for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients and join medians, x_max, and x_min
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  coef_df <- left_join(coef_df, medians_df, by = "species")
  coef_df <- left_join(coef_df, max_x_df, by = "species")
  coef_df <- left_join(coef_df, min_x_df, by = "species")
  
  # Calculate x50: the x value corresponding to the species-specific median NDVI
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((median_value - a) > 0 & b > 0,
                        -log((median_value - a)/b) / c,
                        NA))
  
  #### PANEL A: Observed Data, Fitted Nonlinear Curves, and Dashed Linear Regression Lines ####
  # Get nonlinear model predictions
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
  
  # Compute endpoints for median-to-maximum regression line
  linear_endpoints_max <- coef_df %>%
    mutate(x_upper = x_max,
           y_upper = a + b * exp(-c * x_max)) %>%
    select(species, x50, median_value, x_upper, y_upper)
  
  linear_lines_max <- linear_endpoints_max %>%
    pivot_longer(cols = c(x50, x_upper), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_upper))
  
  # Compute endpoints for median-to-minimum regression line
  linear_endpoints_min <- coef_df %>%
    mutate(x_lower = x_min,
           y_lower = a + b * exp(-c * x_min)) %>%
    select(species, x50, median_value, x_lower, y_lower)
  
  linear_lines_min <- linear_endpoints_min %>%
    pivot_longer(cols = c(x50, x_lower), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_lower))
  
  # Reverse the x-axis so that soil water potential is shown as positive
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed regression lines from median to maximum and median to minimum
    geom_line(data = linear_lines_max, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = linear_lines_min, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    # Add a larger black symbol at the median point (x50, median NDVI)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes (Both Regressions) ####
  # Compute slope statistics for median-to-maximum regression
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_max
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_max <- bind_rows(stats_list_max)
  stats_df_max$species <- factor(stats_df_max$species, levels = species_order)
  stats_df_max <- stats_df_max %>%
    mutate(RegType = "Max",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Compute slope statistics for median-to-minimum regression
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_min
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_min <- bind_rows(stats_list_min)
  stats_df_min$species <- factor(stats_df_min$species, levels = species_order)
  stats_df_min <- stats_df_min %>%
    mutate(RegType = "Min",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Combine the two sets of slope statistics
  stats_df_combined <- bind_rows(stats_df_max, stats_df_min)
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope_lm, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope_lm - se, ymax = slope_lm + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope_lm/2), 
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and the combined Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_slope_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
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
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.3f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2), 
              color = "black", size = 3.5, 
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles",
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_slope_each <- function(data, coef_fig, output_figure) {
  
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept for future needs, e.g., deltaMethod)
  library(broom)
  library(tidyr)
  library(ggpattern)  # for patterned bars in Panel C
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum transpiration deficit for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients and join medians, x_max, and x_min
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  coef_df <- left_join(coef_df, medians_df, by = "species")
  coef_df <- left_join(coef_df, max_x_df, by = "species")
  coef_df <- left_join(coef_df, min_x_df, by = "species")
  
  # Calculate x50: the x value corresponding to the species-specific median NDVI
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((median_value - a) > 0 & b > 0,
                        -log((median_value - a)/b) / c,
                        NA))
  
  #### PANEL A: Observed Data, Fitted Nonlinear Curves, and Dashed Linear Regression Lines ####
  # Get nonlinear model predictions
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
  
  # Compute endpoints for median-to-maximum regression line
  linear_endpoints_max <- coef_df %>%
    mutate(x_upper = x_max,
           y_upper = a + b * exp(-c * x_max)) %>%
    select(species, x50, median_value, x_upper, y_upper)
  
  linear_lines_max <- linear_endpoints_max %>%
    pivot_longer(cols = c(x50, x_upper), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_upper))
  
  # Compute endpoints for median-to-minimum regression line
  linear_endpoints_min <- coef_df %>%
    mutate(x_lower = x_min,
           y_lower = a + b * exp(-c * x_min)) %>%
    select(species, x50, median_value, x_lower, y_lower)
  
  linear_lines_min <- linear_endpoints_min %>%
    pivot_longer(cols = c(x50, x_lower), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_lower))
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed regression lines from median to maximum and median to minimum
    geom_line(data = linear_lines_max, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = linear_lines_min, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    # Add a larger black symbol at the median point (x50, median NDVI)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes (Both Regressions) ####
  # Compute slope statistics for median-to-maximum regression
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_max
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_max <- bind_rows(stats_list_max)
  stats_df_max$species <- factor(stats_df_max$species, levels = species_order)
  stats_df_max <- stats_df_max %>%
    mutate(RegType = "Max",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Compute slope statistics for median-to-minimum regression
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_min
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_min <- bind_rows(stats_list_min)
  stats_df_min$species <- factor(stats_df_min$species, levels = species_order)
  stats_df_min <- stats_df_min %>%
    mutate(RegType = "Min",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Combine the two sets of slope statistics
  stats_df_combined <- bind_rows(stats_df_max, stats_df_min)
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope_lm, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope_lm - se, ymax = slope_lm + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope_lm/2), 
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
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
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.3f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2), 
              color = "black", size = 3.5, 
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles",
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(coef_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = coef_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_PSIbin_slope_each(all_results_df, 
                            save_coeff_fig = "results/Figures/Quantiles_PSIbin_Coeffs.png", 
                            save_slope_fig = "results/Figures/Quantiles_PSIbin.png")

plot_Quantiles_TDiff_slope_each(all_results_df, 
                                coef_fig = "results/Figures/Quantiles_TDiffbin_Coeffs.png", 
                                output_figure = "results/Figures/Quantiles_TDiffbin.png")

load("results/Data/All_Species_Proportions_PSI_TDiff.RData")

plot_NDVI_PDM_PSIbin_slope_each <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept in case of future needs)
  library(tidyr)    # for pivot_longer
  library(ggpattern)  # for patterned bars in Panel C
  
  # Identify the value column and order species; define the color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create a positive soil water potential (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum soil water potential for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 1, b = 0.01, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients and join medians, x_max, and x_min
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  coef_df <- left_join(coef_df, medians_df, by = "species")
  coef_df <- left_join(coef_df, max_x_df, by = "species")
  coef_df <- left_join(coef_df, min_x_df, by = "species")
  
  # Calculate x50: the x value corresponding to the species-specific median NDVI
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((median_value - a) > 0 & b > 0,
                        -log((median_value - a)/b) / c,
                        NA))
  
  #### PANEL A: Observed Data, Fitted Nonlinear Curves, and Dashed Linear Regression Lines ####
  # Get nonlinear model predictions
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
  
  # Compute endpoints for median-to-maximum regression line
  linear_endpoints_max <- coef_df %>%
    mutate(x_upper = x_max,
           y_upper = a + b * exp(-c * x_max)) %>%
    select(species, x50, median_value, x_upper, y_upper)
  
  linear_lines_max <- linear_endpoints_max %>%
    pivot_longer(cols = c(x50, x_upper), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_upper))
  
  # Compute endpoints for median-to-minimum regression line
  linear_endpoints_min <- coef_df %>%
    mutate(x_lower = x_min,
           y_lower = a + b * exp(-c * x_min)) %>%
    select(species, x50, median_value, x_lower, y_lower)
  
  linear_lines_min <- linear_endpoints_min %>%
    pivot_longer(cols = c(x50, x_lower), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_lower))
  
  # Reverse the x-axis so that soil water potential is shown as positive
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed regression lines from median to maximum and median to minimum
    geom_line(data = linear_lines_max, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = linear_lines_min, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    # Add a larger black symbol at the median point (x50, median NDVI)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Proportions") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes (Both Regressions) ####
  # Compute slope statistics for median-to-maximum regression
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_max
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_max <- bind_rows(stats_list_max)
  stats_df_max$species <- factor(stats_df_max$species, levels = species_order)
  stats_df_max <- stats_df_max %>%
    mutate(RegType = "Max",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Compute slope statistics for median-to-minimum regression
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_min
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_min <- bind_rows(stats_list_min)
  stats_df_min$species <- factor(stats_df_min$species, levels = species_order)
  stats_df_min <- stats_df_min %>%
    mutate(RegType = "Min",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Combine the two sets of slope statistics
  stats_df_combined <- bind_rows(stats_df_max, stats_df_min)
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope_lm, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope_lm - se, ymax = slope_lm + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope_lm/2), 
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and the combined Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_slope_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
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
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.3f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2), 
              color = "black", size = 3.5, 
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Proportions",
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Proportions_TDiff_slope_each <- function(data, coef_fig, output_figure) {
  
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept for future needs, e.g., deltaMethod)
  library(broom)
  library(tidyr)
  library(ggpattern)  # for patterned bars in Panel C
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum transpiration deficit for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES WITH SPECIES-SPECIFIC START VALUES ####
  # Define different starting values for each species:
  start_list <- list(
    Oak   = list(a = 1, b = 0.05, c = 0.001),
    Beech = list(a = 1, b = 0.06, c = 0.001),
    Spruce= list(a = 0, b = 0.05, c = 0.001),
    Pine  = list(a = 0, b = 0.05, c = 0.002)
  )
  
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients and join medians, x_max, and x_min
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  coef_df <- left_join(coef_df, medians_df, by = "species")
  coef_df <- left_join(coef_df, max_x_df, by = "species")
  coef_df <- left_join(coef_df, min_x_df, by = "species")
  
  # Calculate x50: the x value corresponding to the species-specific median NDVI
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((median_value - a) > 0 & b > 0,
                        -log((median_value - a)/b) / c,
                        NA))
  
  #### PANEL A: Observed Data, Fitted Nonlinear Curves, and Dashed Linear Regression Lines ####
  # Get nonlinear model predictions
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
  
  # Compute endpoints for median-to-maximum regression line
  linear_endpoints_max <- coef_df %>%
    mutate(x_upper = x_max,
           y_upper = a + b * exp(-c * x_max)) %>%
    select(species, x50, median_value, x_upper, y_upper)
  
  linear_lines_max <- linear_endpoints_max %>%
    pivot_longer(cols = c(x50, x_upper), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_upper))
  
  # Compute endpoints for median-to-minimum regression line
  linear_endpoints_min <- coef_df %>%
    mutate(x_lower = x_min,
           y_lower = a + b * exp(-c * x_min)) %>%
    select(species, x50, median_value, x_lower, y_lower)
  
  linear_lines_min <- linear_endpoints_min %>%
    pivot_longer(cols = c(x50, x_lower), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_lower))
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed regression lines from median to maximum and median to minimum
    geom_line(data = linear_lines_max, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = linear_lines_min, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    # Add a larger black symbol at the median point (x50, median NDVI)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Proportions") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes (Both Regressions) ####
  # Compute slope statistics for median-to-maximum regression
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_max
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_max <- bind_rows(stats_list_max)
  stats_df_max$species <- factor(stats_df_max$species, levels = species_order)
  stats_df_max <- stats_df_max %>%
    mutate(RegType = "Max",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Compute slope statistics for median-to-minimum regression
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_min
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_min <- bind_rows(stats_list_min)
  stats_df_min$species <- factor(stats_df_min$species, levels = species_order)
  stats_df_min <- stats_df_min %>%
    mutate(RegType = "Min",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Combine the two sets of slope statistics
  stats_df_combined <- bind_rows(stats_df_max, stats_df_min)
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope_lm, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope_lm - se, ymax = slope_lm + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope_lm/2), 
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
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
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.3f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2), 
              color = "black", size = 3.5, 
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Proportions",
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(coef_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = coef_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_PDM_PSIbin_slope_each(all_results_df, 
                                save_coeff_fig = "results/Figures/Proportions_PSIbin_Coeffs.png", 
                                save_slope_fig = "results/Figures/Proportions_PSIbin.png")

plot_Proportions_TDiff_slope_each(all_results_df, 
                                  coef_fig = "results/Figures/Proportions_TDiffbin_Coeffs.png", 
                                  output_figure = "results/Figures/Proportions_TDiffbin.png")


plot_Proportions_TDiff_slope_linear_each <- function(data, coef_fig, output_figure) {
  
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)    # for lmList
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(ggpattern)  # for patterned bars in Panel C
  library(broom)      # to tidy lm objects
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x)) %>% droplevels()
  
  #### Compute species-specific summary statistics ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(avg_value),
              x_median = median(x))
  
  max_x_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(x_max = max(x))
  
  min_x_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(x_min = min(x))
  
  #### LINEAR MODELING PER SPECIES (Overall) ####
  # Fit an overall linear model for each species
  lm_models <- lmList(avg_value ~ x | species, data = data_clean)
  print(summary(lm_models))
  
  # Extract overall coefficients: intercept (a) and slope (b)
  coef_df <- as.data.frame(coef(lm_models), optional = TRUE) %>%
    rownames_to_column(var = "species") %>%
    filter(!is.na(`(Intercept)`)) %>%
    rename(a = `(Intercept)`, b = x)
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Join with summary statistics
  coef_df <- coef_df %>% 
    left_join(medians_df, by = "species") %>%
    left_join(max_x_df, by = "species") %>%
    left_join(min_x_df, by = "species")
  
  # Calculate x50: x value where the overall predicted value equals the species-specific median
  coef_df <- coef_df %>%
    mutate(x50 = ifelse(b != 0, (median_value - a) / b, NA))
  
  #### PANEL A: Observed Data, Fitted Overall Linear Curves, and Dashed Segment Regressions ####
  # Get predictions for the overall linear model for each species
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_coef <- coef_df %>% filter(species == sp)
      pred <- sp_coef$a + sp_coef$b * x_seq
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  #### Segment Regressions: Fit separate linear models on each side of the median ####
  # Max side: for x >= x50
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    df_sp <- data_clean %>% filter(species == sp)
    x50_val <- coef_df %>% filter(species == sp) %>% pull(x50)
    df_subset <- df_sp %>% filter(x >= x50_val)
    if(nrow(df_subset) < 2) {
      row <- coef_df %>% filter(species == sp)
      x_pts <- c(row$x50, row$x_max)
      y_pts <- c(row$median_value, row$a + row$b * row$x_max)
      lm_fit <- lm(y_pts ~ x_pts)
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
    }
    tibble(species = sp,
           intercept_max = coef(lm_fit)[["(Intercept)"]],
           slope_max = coef(lm_fit)[["x"]],
           se_max = summary(lm_fit)$coefficients["x", "Std. Error"],
           p_val_max = summary(lm_fit)$coefficients["x", "Pr(>|t|)"])
  })
  stats_df_max <- bind_rows(stats_list_max)
  
  # Min side: for x <= x50
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    df_sp <- data_clean %>% filter(species == sp)
    x50_val <- coef_df %>% filter(species == sp) %>% pull(x50)
    df_subset <- df_sp %>% filter(x <= x50_val)
    if(nrow(df_subset) < 2) {
      row <- coef_df %>% filter(species == sp)
      x_pts <- c(row$x50, row$x_min)
      y_pts <- c(row$median_value, row$a + row$b * row$x_min)
      lm_fit <- lm(y_pts ~ x_pts)
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
    }
    tibble(species = sp,
           intercept_min = coef(lm_fit)[["(Intercept)"]],
           slope_min = coef(lm_fit)[["x"]],
           se_min = summary(lm_fit)$coefficients["x", "Std. Error"],
           p_val_min = summary(lm_fit)$coefficients["x", "Pr(>|t|)"])
  })
  stats_df_min <- bind_rows(stats_list_min)
  
  # Create data frames for the dashed line segments using the segment regression results
  seg_max_list <- lapply(levels(data_clean$species), function(sp) {
    row_coef <- coef_df %>% filter(species == sp)
    row_seg <- stats_df_max %>% filter(species == sp)
    x_seq <- seq(row_coef$x50, row_coef$x_max, length.out = 50)
    tibble(species = sp, x = x_seq, y = row_seg$intercept_max + row_seg$slope_max * x_seq)
  })
  seg_max_df <- bind_rows(seg_max_list)
  
  seg_min_list <- lapply(levels(data_clean$species), function(sp) {
    row_coef <- coef_df %>% filter(species == sp)
    row_seg <- stats_df_min %>% filter(species == sp)
    x_seq <- seq(row_coef$x_min, row_coef$x50, length.out = 50)
    tibble(species = sp, x = x_seq, y = row_seg$intercept_min + row_seg$slope_min * x_seq)
  })
  seg_min_df <- bind_rows(seg_min_list)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed segments from the segment regressions
    geom_line(data = seg_max_df, aes(x = x, y = y, group = species),
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = seg_min_df, aes(x = x, y = y, group = species),
              color = "black", size = 0.5, linetype = "dashed") +
    # Emphasize the median point with a triangle marker (shape = 17)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Proportions") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit (x50)") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes of Segment Regressions ####
  # Combine slope stats from the two segment regressions
  stats_df_combined <- bind_rows(
    stats_df_max %>% rename(slope = slope_max, se = se_max, p_val = p_val_max) %>% mutate(RegType = "Max"),
    stats_df_min %>% rename(slope = slope_min, se = se_min, p_val = p_val_min) %>% mutate(RegType = "Min")
  ) %>%
    mutate(label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope - se, ymax = slope + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope/2),
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  right_side <- p_x50 / p_slope_combined
  final_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_plot)
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (Intercept and Slope) AS A SEPARATE FIGURE ####
  # Use broom to tidy the individual lm models so we can extract p-values.
  tidy_lm <- map_df(lm_models, broom::tidy, .id = "species") %>%
    mutate(species = factor(species, levels = species_order),
           Coefficient = if_else(term == "(Intercept)", "a", "b"),
           label = if_else(p.value < 0.05, "*", sprintf("%.3f", p.value)))
  
  # Create a coefficient dataframe for plotting
  coeffs_long <- coef_df %>%
    select(species, a, b) %>%
    pivot_longer(cols = c("a", "b"), names_to = "Coefficient", values_to = "Value")
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  # Join with tidy lm to get p-value labels
  coeffs_long <- left_join(coeffs_long, tidy_lm %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 3.5, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Proportions (Linear Model)",
         x = "Species", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text = element_text(color = "black", size = 10),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
          axis.title = element_text(face = "bold", size = 12),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12))
  print(p_coeffs)
  
  dir.create(dirname(coef_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = coef_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Proportions_TDiff_slope_linear_each(all_results_df, 
                                         coef_fig = "results/Figures/Proportions_TDiffbin_Coeffs.png", 
                                         output_figure = "results/Figures/Proportions_TDiffbin.png")
