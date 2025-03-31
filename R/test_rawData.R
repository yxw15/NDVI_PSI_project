setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

plot_orig_quantiles_PSI_exp_slope <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Use the original data (no binning) and remove missing values
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(tidyr)    # for pivot_longer
  
  # Set the dependent variable column to the original "Quantiles"
  value_col <- "Quantiles"
  
  # Define species order (modify if needed)
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Define a color palette for the species
  cb_palette <- c("Oak"    = "#E69F00",   # Orange
                  "Beech"  = "#0072B2",   # Deep blue
                  "Spruce" = "#009E73",   # Bluish-green
                  "Pine"   = "#F0E442")   # Yellow
  
  # Create a positive soil water potential variable by multiplying by -1
  data <- data %>% mutate(x = -soil_water_potential)
  
  # Clean data: remove rows with missing or non-finite values in Quantiles or x
  data <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set a threshold manually (adjust if needed)
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  # Starting parameters and control settings for the nls models
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Fit the nonlinear model for each species using the raw data
  nls_models <- nlsList(Quantiles ~ a + b * exp(-c * x) | species,
                        data = data,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients for each species and enforce species order
  # Only keep species with valid (non-NA) coefficients
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Identify valid species that converged
  valid_species <- coef_df$species
  
  # Calculate x50 (the x value where Quantiles equals threshold) and the slope at x50
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a)/b) / c,
                        NA),
           slope50 = -c * (threshold - a))
  
  #### PANEL A: Observed Data and Fitted Curves ####
  # Subset data to only include valid species
  data_valid <- data %>% filter(species %in% valid_species)
  
  # Generate predicted values over a sequence of x values for each valid species
  pred_list <- data_valid %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  # Reverse the x-axis so the displayed values reflect the original (negative) water potential
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    geom_point(data = data_valid, aes(x = x, y = Quantiles, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 2000, y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Annotations ####
  # Use only valid species for statistics
  stats_list <- lapply(valid_species, function(sp) {
    mod <- nls_models[[as.character(sp)]]
    coefs <- coef(mod)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    
    # Compute standard error using deltaMethod for c*(a - threshold)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    t_val <- slope50_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
    # Compute R² for the species model
    df_sp <- data_valid %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    
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
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
  
  #### COMBINE PANELS A, B, and C INTO THE SLOPE FIGURE ####
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### COEFFICIENT PLOT: Bar Plot of Model Coefficients ####
  coeff_stats_list <- lapply(valid_species, function(sp) {
    mod <- nls_models[[as.character(sp)]]
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
                          sprintf("%.2f", `Pr(>|t|)`)))
  
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
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles with PSI",
         subtitle = expression(italic(NDVI) == a + b * e^(-c * italic(x))),
         x = "",
         y = "Coefficient Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
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
  
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_orig_quantiles_PSI_exp_slope(all_results_df,
                                  "results/key_displays/NDVI_Q_PSI_exp_coeff.png",
                                  "results/key_displays/NDVI_Q_PSI_exp_slope.png")

plot_orig_quantiles_PSI_linear_slope <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Use the original data (no binning) and remove missing values
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)      # for lmList
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)       # for deltaMethod (not needed here but kept for consistency)
  library(tidyr)     # for pivot_longer
  
  # Set the dependent variable column to the original "Quantiles"
  value_col <- "Quantiles"
  
  # Define species order (modify if needed)
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Define a color palette for the species
  cb_palette <- c("Oak"    = "#E69F00",   # Orange
                  "Beech"  = "#0072B2",   # Deep blue
                  "Spruce" = "#009E73",   # Bluish-green
                  "Pine"   = "#F0E442")   # Yellow
  
  # Create a positive soil water potential variable by multiplying by -1
  data <- data %>% mutate(x = -soil_water_potential)
  
  # Clean data: remove rows with missing or non-finite values in Quantiles or x
  data <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set a threshold manually (adjust if needed)
  threshold <- 11.5
  
  #### LINEAR MODELING PER SPECIES ####
  # Fit a linear model for each species using lmList
  lm_models <- lmList(Quantiles ~ x | species, data = data)
  print(summary(lm_models))
  
  # Extract coefficients for each species and enforce species order
  # The coefficients are named "(Intercept)" and "x"; we rename them to a and b.
  coef_df <- as.data.frame(coef(lm_models)) %>% 
    rownames_to_column(var = "species")
  names(coef_df)[names(coef_df) == "(Intercept)"] <- "a"
  names(coef_df)[names(coef_df) == "x"] <- "b"
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Identify valid species (all species with available coefficients)
  valid_species <- coef_df$species
  
  # Calculate x50 (the x value where Quantiles equals threshold) and the slope at x50
  # For the linear model: a + b*x50 = threshold, so x50 = (threshold - a) / b.
  # The original soil water potential is -x50.
  coef_df <- coef_df %>%
    mutate(x50 = (threshold - a) / b,
           slope50 = b)
  
  #### PANEL A: Observed Data and Fitted Lines ####
  # Subset data to only include valid species
  data_valid <- data %>% filter(species %in% valid_species)
  
  # Generate predicted values over a sequence of x values for each valid species
  pred_list <- data_valid %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- lm_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  # Reverse the x-axis so the displayed values reflect the original (negative) water potential
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    # geom_point(data = data_valid, aes(x = x, y = Quantiles, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 2000, y = threshold, label = "median", 
             hjust = -1, vjust = -0.3, fontface = "italic", size = 4) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  # Plot original soil water potential at threshold, which is -x50.
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Bar Plot of Slope Estimates with Error Bars and Annotations ####
  # For the linear model, the slope is simply b.
  stats_list <- lapply(valid_species, function(sp) {
    mod <- lm_models[[as.character(sp)]]
    s <- summary(mod)$coefficients
    slope_est <- s["x", "Estimate"]
    se <- s["x", "Std. Error"]
    df_mod <- mod$df.residual
    t_val <- slope_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    r_squared <- summary(mod)$r.squared
    tibble(species = sp,
           slope_est = slope_est,
           slope_abs = abs(slope_est),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
  
  #### COMBINE PANELS A, B, and C INTO THE SLOPE FIGURE ####
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### COEFFICIENT PLOT: Bar Plot of Model Coefficients ####
  # Extract coefficient statistics from each valid model
  coeff_stats_list <- lapply(valid_species, function(sp) {
    mod <- lm_models[[as.character(sp)]]
    s <- summary(mod)$coefficients
    s <- as.data.frame(s)
    s$Coefficient <- rownames(s)
    s$species <- sp
    # Rename coefficients for consistency: (Intercept) -> a, x -> b
    s$Coefficient <- ifelse(s$Coefficient == "(Intercept)", "a", "b")
    s
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  
  coeff_stats <- coeff_stats %>%
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.2f", `Pr(>|t|)`)))
  
  # Pivot the coefficient data into long format and merge labels
  coeffs_long <- coef_df %>%
    select(species, a, b) %>%
    pivot_longer(cols = c("a", "b"),
                 names_to = "Coefficient",
                 values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles with PSI (Linear Model)",
         subtitle = expression(italic(NDVI) == a + b %.% italic(x)),
         x = "",
         y = "Coefficient Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
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
  
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_orig_quantiles_PSI_linear_slope(all_results_df,
                                     "results/key_displays/NDVI_Q_PSI_linear_coeff.png",
                                     "results/key_displays/NDVI_Q_PSI_linear_slope.png")

plot_orig_quantiles_TDiff_exp_slope_coeff <- function(data, coef_fig, output_figure) {
  # Use raw data and remove missing values
  data <- na.omit(data)
  
  # Load required libraries (assumed to be installed)
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(broom)
  library(tidyr)
  
  # Set the dependent variable column to the raw "Quantiles"
  value_col <- "Quantiles"
  
  # Define species order and convert species to factor
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use the raw value: set x = transpiration_deficit
  data <- data %>% mutate(x = transpiration_deficit)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set a threshold (using a fixed value here; adjust as needed)
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  # Fit the exponential model for each species
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  nls_models <- nlsList(Quantiles ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficient summaries for a, b, and c from the model
  model_summary <- summary(nls_models)
  coeff_a <- model_summary$coefficients[, , "a"]
  coeff_b <- model_summary$coefficients[, , "b"]
  coeff_c <- model_summary$coefficients[, , "c"]
  
  # Create a data frame with species names, estimates, and p-values
  results_df <- data.frame(
    species    = rownames(coeff_a),
    a_estimate = coeff_a[, "Estimate"],
    a_pvalue   = coeff_a[, "Pr(>|t|)"],
    b_estimate = coeff_b[, "Estimate"],
    b_pvalue   = coeff_b[, "Pr(>|t|)"],
    c_estimate = coeff_c[, "Estimate"],
    c_pvalue   = coeff_c[, "Pr(>|t|)"]
  )
  results_df$species <- factor(results_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Reshape the coefficient results into long format and merge p-values for labeling
  coeffs_long <- results_df %>%
    pivot_longer(
      cols = c(a_estimate, b_estimate, c_estimate),
      names_to = "Coefficient",
      values_to = "Value",
      names_pattern = "([abc])_.*"
    ) %>%
    left_join(
      results_df %>%
        pivot_longer(
          cols = c(a_pvalue, b_pvalue, c_pvalue),
          names_to = "Coefficient",
          values_to = "pvalue",
          names_pattern = "([abc])_.*"
        ),
      by = c("species", "Coefficient")
    ) %>%
    mutate(
      label = if_else(pvalue < 0.05, "*", sprintf("%.3f", pvalue))
    )
  
  # Create the faceted coefficient bar plot with p-value labels
  p_coeff <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles with TDiff (Raw Data)",
         subtitle = expression(italic(NDVI) == a + b * e^(-c * italic(x))),
         x = "",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
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
  print(p_coeff)
  
  # Save the coefficient plot to file
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # Additional analysis and combined panels
  # -------------------------------
  
  # Extract coefficients for a, b, and c per species from the model; keep only species with valid models
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Calculate x50 and slope50 from exponential parameters
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a)/b) / c,
                        NA),
           slope50 = -c * (threshold - a))
  
  # Generate predictions for each species (Panel A) – skip species if model is NULL.
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      sp_model <- nls_models[[as.character(sp)]]
      if (is.null(sp_model)) {
        return(data.frame())
      } else {
        x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
        pred <- predict(sp_model, newdata = data.frame(x = x_seq))
        data.frame(x = x_seq, pred = pred)
      }
    })
  pred_all <- bind_rows(pred_list)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = Quantiles, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -2.1, vjust = -0.3, fontface = "italic", size = 4) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # Panel B: Bar plot of x50 values
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Panel C: Bar plot of absolute slope at x50 with error bars and R² annotations.
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    if (is.null(mod)) return(NULL)
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
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    
    tibble(species = sp,
           slope50 = slope50_est,
           slope_abs = abs(slope50_est),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_list <- stats_list[!sapply(stats_list, is.null)]
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  stats_df <- stats_df %>%
    mutate(label_text = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
  
  # Combine Panels: Panel (a) on the left; Panels (b) and (c) in the right column.
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot to file
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_orig_quantiles_TDiff_exp_slope_coeff(all_results_df,
                                          "results/key_displays/NDVI_Q_TDiff_exp_coeff.png",
                                          "results/key_displays/NDVI_Q_TDiff_exp_slope.png")

plot_orig_quantiles_TDiff_linear_slope_coeff <- function(data, coef_fig, output_figure) {
  # Use raw data and remove missing values
  data <- na.omit(data)
  
  # Load required libraries (assumed to be installed)
  library(ggplot2)
  library(nlme)      # for lmList
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)
  library(broom)
  library(tidyr)
  
  # Set the dependent variable column to the raw "Quantiles"
  value_col <- "Quantiles"
  
  # Define species order and convert species to factor
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use the raw value: set x = transpiration_deficit
  data <- data %>% mutate(x = transpiration_deficit)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set a threshold (using a fixed value here; adjust as needed)
  threshold <- 11.5
  
  #### LINEAR MODELING PER SPECIES ####
  # Fit a linear model: Quantiles ~ x, by species
  lm_models <- lmList(Quantiles ~ x | species, data = data_clean)
  print(summary(lm_models))
  
  # Extract coefficients for each species and compute p-values for intercept and slope.
  # We'll build a results_df by looping over valid species.
  valid_species <- levels(data_clean$species)
  
  results_df <- do.call(rbind, lapply(valid_species, function(sp) {
    mod <- lm_models[[sp]]
    summ <- summary(mod)$coefficients
    data.frame(
      species    = sp,
      a_estimate = summ["(Intercept)", "Estimate"],
      a_pvalue   = summ["(Intercept)", "Pr(>|t|)"],
      b_estimate = summ["x", "Estimate"],
      b_pvalue   = summ["x", "Pr(>|t|)"],
      b_se       = summ["x", "Std. Error"],
      r_squared  = summary(mod)$r.squared,
      stringsAsFactors = FALSE
    )
  }))
  results_df$species <- factor(results_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Reshape the coefficient results into long format for plotting
  coeffs_long <- results_df %>%
    pivot_longer(
      cols = c(a_estimate, b_estimate),
      names_to = "Coefficient",
      values_to = "Value",
      names_pattern = "([ab])_.*"
    ) %>%
    left_join(
      results_df %>%
        pivot_longer(
          cols = c(a_pvalue, b_pvalue),
          names_to = "Coefficient",
          values_to = "pvalue",
          names_pattern = "([ab])_.*"
        ),
      by = c("species", "Coefficient")
    ) %>%
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.3f", pvalue)))
  
  # Create the faceted coefficient bar plot with p-value labels
  p_coeff <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles with TDiff (Raw Data)",
         subtitle = expression(italic(NDVI) == a + b %.% italic(x)),
         x = "",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
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
  print(p_coeff)
  
  # Save the coefficient plot to file
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # Additional analysis and combined panels
  # -------------------------------
  
  # For the linear model, calculate x50 and use the slope directly.
  # a + b*x50 = threshold  ⟹  x50 = (threshold - a)/b
  results_df <- results_df %>% mutate(x50 = (threshold - a_estimate) / b_estimate,
                                      slope = b_estimate)
  
  # Generate predictions for each species (Panel A)
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- lm_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  # Reverse x-axis if desired (here we leave x as transpiration deficit; adjust if needed)
  p_combined <- ggplot() +
    # geom_point(data = data_clean, aes(x = x, y = Quantiles, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -3, vjust = -0.3, fontface = "italic", size = 4) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # Panel B: Bar plot of x50 values (here, x50 represents the transpiration deficit at threshold)
  p_x50 <- ggplot(results_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit at Threshold") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Panel C: Bar plot of absolute slope values with error bars and R² annotations
  stats_df <- results_df %>% mutate(slope_abs = abs(slope))
  
  # Add p-values and standard errors from the slope coefficient for annotation
  stats_df <- stats_df %>%
    rowwise() %>%
    mutate(
      slope_se = summary(lm_models[[as.character(species)]])$coefficients["x", "Std. Error"],
      df_mod = lm_models[[as.character(species)]]$df.residual,
      t_val = slope / slope_se,
      p_val = 2 * (1 - pt(abs(t_val), df_mod))
    ) %>%
    ungroup() %>%
    mutate(label_text = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - slope_se), ymax = slope_abs + slope_se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot to file
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_orig_quantiles_TDiff_linear_slope_coeff(all_results_df,
                                             "results/key_displays/NDVI_Q_TDiff_linear_coeff.png",
                                             "results/key_displays/NDVI_Q_TDiff_linear_slope.png")

plot_orig_TDiff_PSI_poly_3_slope <- function(data, coef_output, figure_output) {
  # Use raw data and remove missing values
  df <- na.omit(data)
  
  # Rename raw columns to match what the original function expects.
  # Here we assume that the original "soil_water_potential" becomes "bin_median"
  # and "transpiration_deficit" becomes "avg_transpiration_deficit".
  df <- df %>% 
    dplyr::rename(bin_median = soil_water_potential,
                  avg_transpiration_deficit = transpiration_deficit)
  
  # Load required libraries (if not already loaded)
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(gridExtra) # For arranging multiple plots
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  
  # Define a custom color palette and set species order.
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_levels)
  
  ## Panel A: Mixed-Effects Model Plot
  
  # Fit the mixed-effects model with a third-order polynomial for bin_median.
  model <- lmer(avg_transpiration_deficit ~ poly(bin_median, 3) + 
                  (poly(bin_median, 3) | species),
                data = df)
  
  # Create prediction data for each species.
  pred_data <- df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median)) %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Predict the fitted values using the model (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Calculate the overall mean and median of transpiration deficit from the raw data.
  line_val <- mean(data$transpiration_deficit)
  line_median <- median(data$transpiration_deficit)
  
  # Create the mixed-effects model plot with horizontal lines for mean and median.
  plot_mixed <- ggplot(df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), linewidth = 1) +
    geom_hline(yintercept = line_val, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(df$bin_median, na.rm = TRUE), 
             y = line_val, 
             label = paste0("mean: ", round(line_val, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(df$bin_median, na.rm = TRUE), 
             y = line_median, 
             label = paste0("median: ", round(line_median, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential", y = "Transpiration Deficit") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 16)
    )
  
  ## Panel B: Bar Plot of x-values at Mean Transpiration Deficit
  
  x_at_mean <- pred_data %>%
    group_by(species) %>%
    summarise(x_at_mean = bin_median[which.min(abs(predicted - line_val))])
  x_at_mean$species <- factor(x_at_mean$species, levels = species_levels)
  
  p_bar_x <- ggplot(x_at_mean, aes(x = species, y = x_at_mean, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "PSI at Mean TDiff") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## Panel C: Bar Plot of Local Slopes at 30% of Maximum Transpiration Deficit
  
  local_slope_data <- pred_data %>%
    group_by(species) %>%
    do({
      df_sub <- .
      max_pred <- max(df_sub$predicted, na.rm = TRUE)
      thresh_val <- 0.3 * max_pred
      idx <- which.min(abs(df_sub$predicted - thresh_val))
      win_start <- max(1, idx - 5)
      win_end <- min(nrow(df_sub), idx + 5)
      df_window <- df_sub[win_start:win_end, ]
      lm_local <- lm(predicted ~ bin_median, data = df_window)
      summ <- summary(lm_local)
      slope_val <- coef(lm_local)["bin_median"]
      slope_se <- summ$coefficients["bin_median", "Std. Error"]
      p_val <- summ$coefficients["bin_median", "Pr(>|t|)"]
      r2_val <- summ$r.squared
      data.frame(threshold_value = thresh_val,
                 slope = slope_val,
                 slope_se = slope_se,
                 p_value = p_val,
                 r2 = r2_val)
    }) %>%
    ungroup()
  
  local_slope_data <- local_slope_data %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               sprintf("%.2f*", r2),
                               sprintf("p = %.2f\nR² = %.2f", p_value, r2)))
  
  p_bar <- ggplot(local_slope_data, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se),
                  width = 0.2, color = "black") +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 5, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "Absolute Slope at Mean PSI") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## Combine Panels: Panel A on the left, Panels B (top) and C (bottom) on the right.
  right_panel <- p_bar_x / p_bar
  final_plot <- plot_mixed + right_panel + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined plot with width = 10 and height = 8.
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  ## New Figure: Grouped Bar Plot for All Model Coefficients per Species
  # Extract the species-specific (conditional) coefficients.
  species_coef <- coef(model)$species
  
  # Convert rownames to a proper column and pivot the data to long format.
  coeff_data <- species_coef %>%
    tibble::rownames_to_column("species") %>%
    pivot_longer(cols = -species, names_to = "term", values_to = "value") %>%
    mutate(species = factor(species, levels = species_levels),
           term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 3)1" = "b",
                                "poly(bin_median, 3)2" = "c",
                                "poly(bin_median, 3)3" = "d"))
  
  # Set the order of the term factor.
  coeff_data$term <- factor(coeff_data$term, levels = c("a", "b", "c", "d"))
  
  ## Compute p-values for coefficients by fitting separate linear models for each species.
  species_list <- species_levels
  coeff_stats_list <- lapply(species_list, function(sp) {
    subdata <- subset(df, species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 3), data = subdata)
    summ <- summary(mod_sp)$coefficients
    df_mod <- as.data.frame(summ)
    df_mod$term <- rownames(df_mod)
    df_mod$species <- sp
    df_mod
  })
  coeff_stats <- do.call(rbind, coeff_stats_list)
  coeff_stats <- coeff_stats %>%
    mutate(term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 3)1" = "b",
                                "poly(bin_median, 3)2" = "c",
                                "poly(bin_median, 3)3" = "d")) %>%
    dplyr::select(species, term, p_value = `Pr(>|t|)`)
  
  # Join the p-values with the coefficient data.
  coeff_data <- left_join(coeff_data, coeff_stats, by = c("species", "term"))
  
  # Create a label: if p < 0.05 then "*" else the p_value with 2 decimals.
  coeff_data <- coeff_data %>%
    mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)))
  
  coeff_data$species <- factor(coeff_data$species, levels = species_levels)
  
  # Create the grouped bar plot with p-value labels centered in each bar.
  plot_coeff <- ggplot(coeff_data, aes(x = species, y = value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = label_text, y = value/2), 
              color = "black", size = 5, 
              position = position_dodge(width = 0.8)) +
    facet_wrap(~ term, scales = "free_y") +
    scale_fill_manual(values = cb_palette) +
    labs(title = "Coefficients of Transpiration Deficit (y) with Soil Water Potential (x)", 
         subtitle = expression(hat(Y) == a + b*x + c*x^2 + d*x^3),
         x = "Coefficient Term", 
         y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.title = element_text(face = "bold", size = 20),
          axis.text = element_text(color = "black", size = 16),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 16),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold", size = 16))
  
  print(plot_coeff)
  
  # Save the coefficient plot with width = 10 and height = 8.
  ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
  
  print(coeff_data)
  # Return both plots as a list (if desired)
  return(list(combined_plot = final_plot, coeff_plot = plot_coeff))
}

plot_orig_TDiff_PSI_poly_3_slope(all_results_df,
                                 "results/key_displays/PSI_TDiff_exp_coeff.png",
                                 "results/key_displays/PSI_TDiff_exp_slope.png")

plot_orig_TDiff_PSI_linear_slope <- function(data, coef_output, figure_output) {
  # Use raw data and remove missing values
  df <- na.omit(data)
  
  # Rename raw columns to match what the original function expects.
  # Here we assume that "soil_water_potential" becomes "bin_median"
  # and "transpiration_deficit" becomes "avg_transpiration_deficit".
  df <- df %>% 
    dplyr::rename(bin_median = soil_water_potential,
                  avg_transpiration_deficit = transpiration_deficit)
  
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(gridExtra) # For arranging multiple plots
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  
  # Define a custom color palette and set species order.
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_levels)
  
  ## Panel A: Linear Mixed-Effects Model Plot
  
  # Fit a linear mixed-effects model with a single predictor (bin_median)
  model <- lmer(avg_transpiration_deficit ~ bin_median + (bin_median | species),
                data = df)
  
  # Create prediction data for each species.
  pred_data <- df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median)) %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Predict fitted values using the linear mixed model (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Calculate overall mean and median of transpiration deficit from raw data.
  line_val <- mean(data$transpiration_deficit)
  line_median <- median(data$transpiration_deficit)
  
  # Create the mixed-effects model plot.
  plot_mixed <- ggplot(df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
    #geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), linewidth = 1) +
    geom_hline(yintercept = line_val, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(df$bin_median, na.rm = TRUE), 
             y = line_val, 
             label = paste0("mean: ", round(line_val, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(df$bin_median, na.rm = TRUE), 
             y = line_median, 
             label = paste0("median: ", round(line_median, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential", y = "Transpiration Deficit") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text = element_text(face = "bold", size = 16)
    )
  
  ## Panel B: Bar Plot of PSI Values at Mean Transpiration Deficit
  
  # For each species, find the bin_median value where the predicted value is closest to the overall mean.
  x_at_mean <- pred_data %>%
    group_by(species) %>%
    summarise(x_at_mean = bin_median[which.min(abs(predicted - line_val))])
  x_at_mean$species <- factor(x_at_mean$species, levels = species_levels)
  
  p_bar_x <- ggplot(x_at_mean, aes(x = species, y = x_at_mean, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "PSI at Mean TDiff") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## Panel C: Bar Plot of Slopes from Species–Specific Linear Models
  # For each species, fit an ordinary linear model and extract the slope, its standard error, p–value, and R².
  # Additionally, compute the x value at which the predicted value equals line_val based on the model,
  # verifying that the slope (slope = (line_val - intercept) / x_at_mean) is indeed the model's slope.
  slope_stats <- df %>%
    group_by(species) %>%
    do({
      mod <- lm(avg_transpiration_deficit ~ bin_median, data = .)
      summ <- summary(mod)$coefficients
      intercept <- summ["(Intercept)", "Estimate"]
      slope <- summ["bin_median", "Estimate"]
      # Compute x_at_mean for this species based on the regression line:
      # line_val = intercept + slope * x_at_mean  =>  x_at_mean = (line_val - intercept) / slope
      x_at_mean_model <- (line_val - intercept) / slope
      data.frame(slope = slope,
                 slope_se = summ["bin_median", "Std. Error"],
                 p_val = summ["bin_median", "Pr(>|t|)"],
                 r2 = summary(mod)$r.squared,
                 x_at_mean_model = x_at_mean_model)
    }) %>% ungroup()
  
  # For plotting, we use the absolute value of the slope. Note that because the model is linear,
  # the slope at y = line_val is the same as the estimated regression slope.
  slope_stats <- slope_stats %>% 
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r2),
                               sprintf("p = %.2f\nR² = %.2f", p_val, r2)))
  
  p_slope <- ggplot(slope_stats, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se),
                  width = 0.2, color = "black") +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 5, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Species", y = "Absolute Slope") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## Combine Panels: Panel A on the left; Panels B (top) and C (bottom) on the right.
  right_panel <- p_bar_x / p_slope
  final_plot <- plot_mixed + right_panel + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined plot.
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  ## New Figure: Grouped Bar Plot for All Model Coefficients per Species
  # Extract species-specific coefficients from the mixed model.
  species_coef <- coef(model)$species
  
  # Convert rownames to a column and pivot the data to long format.
  coeff_data <- species_coef %>%
    tibble::rownames_to_column("species") %>%
    pivot_longer(cols = -species, names_to = "term", values_to = "value") %>%
    mutate(species = factor(species, levels = species_levels),
           term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "bin_median" = "b"))
  
  # Set the order of the term factor.
  coeff_data$term <- factor(coeff_data$term, levels = c("a", "b"))
  
  ## Compute p-values for coefficients by fitting separate linear models for each species.
  species_list <- species_levels
  coeff_stats_list <- lapply(species_list, function(sp) {
    subdata <- subset(df, species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ bin_median, data = subdata)
    summ <- summary(mod_sp)$coefficients
    data.frame(species = sp,
               term = c("a", "b"),
               p_value = c(summ["(Intercept)", "Pr(>|t|)"],
                           summ["bin_median", "Pr(>|t|)"]))
  })
  coeff_stats <- do.call(rbind, coeff_stats_list)
  coeff_stats$species <- factor(coeff_stats$species, levels = species_levels)
  
  # Join the p-values with the coefficient data.
  coeff_data <- left_join(coeff_data, coeff_stats, by = c("species", "term"))
  
  # Create a label: if p < 0.05 then "*" else show the p-value with 2 decimals.
  coeff_data <- coeff_data %>%
    mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)))
  
  # Create the grouped bar plot with p-value labels.
  plot_coeff <- ggplot(coeff_data, aes(x = species, y = value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = label_text, y = value/2), 
              color = "black", size = 5, 
              position = position_dodge(width = 0.8)) +
    facet_wrap(~ term, scales = "free_y") +
    scale_fill_manual(values = cb_palette) +
    labs(title = "Coefficients of Transpiration Deficit (y) with Soil Water Potential (x)", 
         subtitle = expression(hat(Y) == a + b*x),
         x = "Coefficient Term", 
         y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 18, face = "italic"),
          axis.title = element_text(face = "bold", size = 20),
          axis.text = element_text(color = "black", size = 16),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 16),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold", size = 16))
  
  print(plot_coeff)
  
  # Save the coefficient plot.
  ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
  
  print(coeff_data)
  # Return both plots as a list (if desired)
  return(list(combined_plot = final_plot, coeff_plot = plot_coeff))
}

plot_orig_TDiff_PSI_linear_slope(all_results_df,
                                 "results/key_displays/PSI_TDiff_linear_coeff.png",
                                 "results/key_displays/PSI_TDiff_linear_slope.png")
