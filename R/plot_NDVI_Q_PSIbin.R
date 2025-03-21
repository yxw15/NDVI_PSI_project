# Script: NDVI Quantiles vs Soil Water Potential Analysis
#
# This script performs the following tasks:
#
# 1. Data Binning:
#    - The function NDVI_PSIbin() bins the soil water potential values into intervals
#      (with a default bin width of 100) and computes the average NDVI quantiles (or proportions)
#      for each bin and species. It also calculates the median soil water potential for each bin.
#
# 2. Model Fitting:
#    - The script fits an exponential decay model of the form:
#         avg_value = a + b * exp(-c * x)
#      where x is defined as the positive soil water potential (x = -bin_median).
#
# 3. Calculation of x50 and slope50:
#    - x50: This value represents the soil water potential (x) at which the fitted NDVI quantile
#           reaches a specified threshold (set manually to 11.5). It is derived by solving:
#
#              threshold = a + b * exp(-c * x)
#
#           for x. Rearranging and solving gives:
#
#              x = -log((threshold - a) / b) / c
#
#           This calculation is performed only when (threshold - a) > 0 and b > 0;
#           otherwise, x50 is set to NA.
#
#    - slope50: This represents the slope (rate of change) of the model at x50.
#           The derivative of the model with respect to x is:
#
#              d(avg_value)/dx = -c * b * exp(-c * x)
#
#           At x = x50, note that b * exp(-c * x50) equals (threshold - a) (from the x50 calculation).
#           Therefore, the slope at x50 is calculated as:
#
#              slope50 = -c * (threshold - a)
#
# 4. Visualization:
#    - The script creates a composite figure consisting of three panels:
#        Panel A: Scatter plot of observed data points with fitted model curves.
#        Panel B: Bar plot of the computed x50 values (soil water potential at the threshold NDVI quantile).
#        Panel C: Bar plot showing the absolute slope at x50 along with error bars and annotations.
#
#    - Additionally, a separate coefficient plot is generated to display the model parameters (a, b, and c)
#      for each species, including significance annotations.
#
# 5. Output:
#    - The final figures (composite slope figure and coefficient figure) are saved to specified file paths.
#
# ---------------------------------------------------------------------------------------

# Load the required data and define functions:
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

NDVI_PSIbin <- function(df, bin_width = 100) {
  
  library(dplyr)
  
  # Identify the correct column dynamically
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Define bin breaks dynamically based on soil water potential range
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean for either Quantiles or Proportions per species
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),  # Dynamic selection
      count = n(),
      .groups = 'drop'
    ) %>%
    filter(count >= 2000) %>%  # Filter bins with at least 2000 observations
    mutate(bin_median = sapply(as.character(PSI_bin), function(bin_label) {
      nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
      mean(nums)
    })) %>%
    select(species, PSI_bin, bin_median, avg_value)  # Keep necessary columns
  
  return(meanNDVI_PSIbin_species)
}

plot_NDVI_Q_PSIbin_exp_slope <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
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
  
  # Set threshold manually to 11.5 instead of computing the median from the data
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  # Starting parameters and control parameters for the nls models
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients for each species into a data frame and ensure species order
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Calculate x50 and the slope at x50 (for Panel B)
  # x50 is computed as the soil water potential value where the model reaches the threshold:
  #   x50 = -log((threshold - a)/b) / c
  # This is valid only if (threshold - a) > 0 and b > 0.
  # The slope (slope50) at x50 is derived from the derivative of the model:
  #   slope50 = -c * (threshold - a)
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a)/b) / c,
                        NA),
           slope50 = -c * (threshold - a))
  
  #### PANEL A: Observed Data and Fitted Curves ####
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
  
  # Reverse the x-axis so soil water potential is positive
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
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
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
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
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  # Create annotation labels for Panel C:
  # If p < 0.05, display R² with an appended star; otherwise, display R² and p-value on a new line.
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
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
  # Get coefficient summary statistics (including p-values) from each nls model
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients  # matrix with Estimate, Std. Error, t value, Pr(>|t|)
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  
  # Create a label: if p < 0.05, label is "*"; else "p = <pvalue>"
  coeff_stats <- coeff_stats %>%
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.2f", `Pr(>|t|)`)))
  
  # Pivot the coefficient data (from coef_df) into long format and merge the labels
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  
  coeffs_long$species <- factor(coeffs_long$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles",
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
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
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

plot_NDVI_Q_PSIbin_exp_slope(all_results_df,
                             "results/key_displays/NDVI_Q_PSIbin_exp_coff.png",
                             "results/key_displays/NDVI_Q_PSIbin_exp_slope.png")
