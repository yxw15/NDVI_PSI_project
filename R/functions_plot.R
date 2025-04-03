setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

df_average_year <- function(df_all) {
  # Determine the correct column name for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Dynamically name the output column
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average quantile/proportion and soil water potential (PSI) per year
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
  # Determine the correct column name for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Dynamically name the output column
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average quantile/proportion and soil water potential (PSI) per year and species
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

plot_time_series_and_correlation_combined <- function(df_all, output_path) {
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  
  # Get species-level averaged data
  df_species <- df_average_year_species(df_all)
  
  # Set species order and assign factor levels
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  # Define the color palette for the species
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Determine NDVI column and label
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in species data.")
  }
  
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "NDVI quantiles"
  } else {
    "NDVI proportions"
  }
  
  # Define a custom theme with consistent text sizes
  custom_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      legend.text = element_text(color = "black", size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption.position = "plot",
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ### TIME SERIES PLOTS ###
  # Panel (a): NDVI Time Series – keep the legend here
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = "NDVI quantiles", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    labs(caption = "(a)")
  
  # Panel (b): Soil Water Potential Time Series – remove legend
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = "Soil water potential", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    labs(caption = "(b)")
  
  # Panel (c): Transpiration Deficit Time Series – remove legend
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = "Transpiration deficit", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    labs(caption = "(c)")
  
  # Combine time series panels vertically
  ts_plot <- p1 / p2 / p3
  
  ### CORRELATION PLOTS ###
  # Compute correlations for panel (d): NDVI vs Soil Water Potential
  annotations_d <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_psi, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  # Compute correlations for panel (e): NDVI vs Transpiration Deficit
  annotations_e <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_tdiff, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  # Create annotation text (ordered by species_order)
  text_d <- paste0(species_order, ": r = ", 
                   annotations_d$corr[match(species_order, annotations_d$species)],
                   collapse = "\n")
  
  text_e <- paste0(species_order, ": r = ", 
                   annotations_e$corr[match(species_order, annotations_e$species)],
                   collapse = "\n")
  
  # Panel (d): NDVI vs Soil Water Potential – remove legend
  p_d <- ggplot(df_species, aes(x = avg_psi, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Soil water potential", y = ndvi_label, caption = "(d)", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = -Inf, y = Inf, label = text_d, hjust = 0, vjust = 1, size = 5, color = "black")
  
  # Panel (e): NDVI vs Transpiration Deficit – remove legend
  p_e <- ggplot(df_species, aes(x = avg_tdiff, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Transpiration deficit", y = "", caption = "(e)", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = Inf, y = Inf, label = text_e, hjust = 1, vjust = 1, size = 5, color = "black")
  
  # Combine correlation panels side-by-side
  corr_plot <- p_d | p_e
  
  # Combine time series (ts_plot) on top and correlation (corr_plot) on bottom;
  # use patchwork's guide collection to pull a single shared legend from p1.
  combined_plot <- ts_plot / corr_plot + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  print(combined_plot)
  
  # Save the combined figure
  ggsave(filename = output_path, plot = combined_plot, width = 10, height = 14, dpi = 300)
}

plot_time_series_and_correlation_combined(all_results_df,
                                          "results/key_displays/time_series_Quantiles_PSI_TDiff_species.png")


NDVI_PSIbin <- function(df, bin_width = 50) {
  
  library(dplyr)
  
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

plot_NDVI_Q_PSIbin_exp_slope_negPSI <- function(data, save_coeff_fig, save_slope_fig) {
  
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
  
  # Keep original (negative) soil water potential
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold manually to 11.5 instead of computing the median from the data
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Compute x50 and slope at x50
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        log((threshold - a)/b) / c,
                        NA),
           slope50 = c * (threshold - a))
  
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
  
  # No transformation — use natural negative axis
  x_scale <- scale_x_continuous()
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x), y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot",
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14) 
    ) +
    labs(caption = "(a)") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
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
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 5) +
    labs(x = "", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### COMBINE PANELS ####
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot ####
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
    scale_fill_manual(values = cb_palette, , name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "",
         # subtitle = expression(italic(NDVI) == a + b * e^{c * italic(x)}),
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
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_exp_slope_negPSI(all_results_df, 
                                    "results/key_displays/NDVI_Q_PSIbin_exp_coeff_negPSI.png",
                                    "results/key_displays/NDVI_Q_PSIbin_exp_slope_negPSI.png")

plot_NDVI_Q_PSIbin_exp_slope_negPSI <- function(data, save_coeff_fig, save_slope_fig) {
  
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
  
  # Keep original (negative) soil water potential
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold manually to 11.5 instead of computing the median from the data
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Compute x50 and slope at x50
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        log((threshold - a)/b) / c,
                        NA),
           slope50 = c * (threshold - a))
  
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
  
  # No transformation — use natural negative axis
  x_scale <- scale_x_continuous()
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x), y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles", caption = "(a)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption.position = "plot",
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
      # Do not set legend.position here; it will be overridden in the combined plot
    ) +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential", caption = "(b)") +
    scale_fill_manual(values = cb_palette, name = "") +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption.position = "plot"
    )
  
  #### PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Annotations ####
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
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 5) +
    labs(x = "", y = "Absolute Slope", caption = "(c)") +
    scale_fill_manual(values = cb_palette, name = "") +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption.position = "plot"
    )
  
  #### COMBINE PANELS ####
  # Here, patchwork collects the legend (from p_combined) and we override the overall theme
  # to place the legend at the bottom center of the final big plot.
  final_slope_plot <- p_combined + (p_x50 / p_slope) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.justification = "center")
  
  print(final_slope_plot)
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot ####
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
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(x = "", y = "Coefficient Value") +
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
      strip.text = element_text(face = "bold", size = 12),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption.position = "plot"
    )
  
  print(p_coeffs)
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

