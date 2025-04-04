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
  library(cowplot)  # For get_legend and plot_grid
  
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
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
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
  
  # Create an x-axis scale that shows breaks every 2 years from 2003 to 2024
  x_scale <- scale_x_continuous(breaks = seq(2003, 2024, by = 2), limits = c(2003, 2024))
  
  ### TIME SERIES PLOTS ###
  # Panel (a): NDVI Time Series – keep the legend (with horizontal layout) for extraction later
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "", guide = guide_legend(nrow = 1)) +
    custom_theme +
    labs(caption = "(a)") +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    x_scale
  
  # Panel (b): Soil Water Potential Time Series – remove legend
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("Soil water potential"), bold("(kPa)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    labs(caption = "(b)") +
    x_scale
  
  # Panel (c): Transpiration Deficit Time Series – remove legend
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("Transpiration deficit"), bold("(mm)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    labs(caption = "(c)") +
    x_scale
  
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
    labs(x = "Soil water potential (kPa)", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), caption = "(d)", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = -Inf, y = Inf, label = text_d, hjust = 0, vjust = 1, size = 5, color = "black")
  
  # Panel (e): NDVI vs Transpiration Deficit – remove legend
  p_e <- ggplot(df_species, aes(x = avg_tdiff, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Transpiration deficit (mm)", y = "", caption = "(e)", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = Inf, y = Inf, label = text_e, hjust = 1, vjust = 1, size = 5, color = "black")
  
  # Combine correlation panels side-by-side
  corr_plot <- p_d | p_e
  
  # Combine time series on top and correlations on bottom, and remove legends from this combination
  combined_main <- ts_plot / corr_plot
  combined_main <- combined_main & theme(legend.position = "top")
  
  # Extract the legend from p1 (now a single horizontal guide)
  legend <- suppressWarnings(cowplot::get_legend(p1))
  
  # Combine the main plot and the legend with the legend underneath
  final_plot <- cowplot::plot_grid(combined_main, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  print(final_plot)
  ggsave(filename = output_path, plot = final_plot, width = 12, height = 14, dpi = 300)
}

plot_yearly_mean_linear_coeffs <- function(df_all, output_path) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  # Get species-level averaged data (assuming df_average_year_species is available)
  df_species <- df_average_year_species(df_all)
  
  # Set species order and assign factor levels
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  # Determine NDVI column and label
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in species data.")
  }
  
  ndvi_label <- if (ndvi_column == "avg_quantile") "NDVI quantiles" else "NDVI proportions"
  
  # Define the color palette for the species
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Define a custom theme similar to your provided theme.
  custom_theme <- theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
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
  
  ### NDVI ~ Soil Water Potential Regression Coefficients (Panel a) ###
  coeffs_psi <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_psi")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, 
               term = names(coef_vals), 
               estimate = coef_vals, 
               p.value = p_vals,
               stringsAsFactors = FALSE)
  }))
  # Rename coefficient terms: intercept -> a, slope -> b
  coeffs_psi$term <- dplyr::recode(coeffs_psi$term, "(Intercept)" = "a", "avg_psi" = "b")
  # Ensure species factor order is maintained
  coeffs_psi$species <- factor(coeffs_psi$species, levels = species_order)
  coeffs_psi$label <- ifelse(coeffs_psi$p.value < 0.05, "*", sprintf("%.2f", coeffs_psi$p.value))
  
  p_coeff_psi <- ggplot(coeffs_psi, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI quantiles ~ soil water potential",
         subtitle = expression(italic(NDVI) == a + b * x + epsilon[i]),
         x = "", y = "Coefficient Value") +
    custom_theme +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ### NDVI ~ Transpiration Deficit Regression Coefficients (Panel b) ###
  coeffs_tdiff <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_tdiff")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, 
               term = names(coef_vals), 
               estimate = coef_vals, 
               p.value = p_vals,
               stringsAsFactors = FALSE)
  }))
  # Rename coefficient terms: intercept -> a, slope -> b
  coeffs_tdiff$term <- dplyr::recode(coeffs_tdiff$term, "(Intercept)" = "a", "avg_tdiff" = "b")
  # Ensure species factor order is maintained
  coeffs_tdiff$species <- factor(coeffs_tdiff$species, levels = species_order)
  coeffs_tdiff$label <- ifelse(coeffs_tdiff$p.value < 0.05, "*", sprintf("%.2f", coeffs_tdiff$p.value))
  
  p_coeff_tdiff <- ggplot(coeffs_tdiff, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI quantile ~ transpiration deficit",
         subtitle = expression(italic(NDVI) == a + b * x + epsilon[i]),
         x = "", y = "") +
    custom_theme +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ### Combine the Two Panels Side-by-Side with a Shared Legend and Overall Caption ###
  combined_plot <- p_coeff_psi + p_coeff_tdiff + 
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(combined_plot)
  ggsave(filename = output_path, plot = combined_plot, width = 14, height = 8, dpi = 300)
}

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

NDVI_TDiffbin <- function(df, bin_width = 3) {
  
  library(dplyr)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Compute total pixels per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks for transpiration_deficit
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  # Create TDiff_bin column
  df <- df %>%
    mutate(TDiff_bin = cut(transpiration_deficit,
                           breaks = bin_breaks,
                           include.lowest = TRUE,
                           right = FALSE))
  
  # Helper function to compute bin median
  get_bin_median <- function(bin_label) {
    nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
    mean(nums)
  }
  
  # Group, summarize, calculate percentage, and filter
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
  
  library(dplyr)
  
  # Identify the correct column
  value_column <- if ("transpiration_deficit" %in% names(df)) "transpiration_deficit"
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks dynamically
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean transpiration_deficit and count per bin per species
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

plot_combined_AIC_R2 <- function(data, save_combined_fig) {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(patchwork)
  
  # Define common species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  #### Panel A: NDVI ~ PSIbin ####
  data_a <- NDVI_PSIbin(data)
  data_a <- na.omit(data_a)
  value_col_a <- "avg_value"
  data_a$species <- factor(data_a$species, levels = species_order)
  data_a <- data_a %>% mutate(x = -bin_median)
  
  start_list_a <- list(a = 5, b = 3, c = 0.001)
  control_params_a <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  aic_results_a <- list()
  for (sp in levels(data_a$species)) {
    sp_data <- data_a %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    # Exponential model
    aic_exp <- NA
    r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data,
          start = start_list_a,
          control = control_params_a)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col_a]] - mean(sp_data[[value_col_a]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Exponential"),
      AIC = c(aic_linear, aic_exp),
      R2 = c(r2_linear, r2_exp)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_a[[sp]] <- sp_results
  }
  aic_df_a <- do.call(rbind, aic_results_a)
  aic_df_a$species <- factor(aic_df_a$species, levels = species_order)
  
  # Shared color palette for Panels A and B
  model_palette_shared <- c("Linear" = "orange",
                            "Exponential" = "dodgerblue")
  
  p_a <- ggplot(aic_df_a, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI quantiles ~ soil water potential") +
    scale_fill_manual(values = model_palette_shared, name = "") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background= element_rect(fill = "white", color = "white"),
      panel.background=element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title     = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      axis.title     = element_text(face = "bold", size = 16),
      axis.text      = element_text(color = "black", size = 14),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text    = element_text(size = 14)
    )+
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  #### Panel B: NDVI ~ TDiffbin ####
  data_b <- NDVI_TDiffbin(data)
  data_b <- na.omit(data_b)
  value_col_b <- "avg_value"
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
    
    aic_exp <- NA
    r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data,
          start = start_list_b,
          control = control_params_b)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col_b]] - mean(sp_data[[value_col_b]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Exponential"),
      AIC = c(aic_linear, aic_exp),
      R2 = c(r2_linear, r2_exp)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_b[[sp]] <- sp_results
  }
  aic_df_b <- do.call(rbind, aic_results_b)
  aic_df_b$species <- factor(aic_df_b$species, levels = species_order)
  
  p_b <- ggplot(aic_df_b, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI quantiles ~ transpiration deficit") +
    scale_fill_manual(values = model_palette_shared, name = "") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background= element_rect(fill = "white", color = "white"),
      panel.background=element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title     = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      axis.title     = element_text(face = "bold", size = 16),
      axis.text      = element_text(color = "black", size = 14),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text    = element_text(size = 14)
    )+
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  #### Panel C: TDiff ~ PSIbin ####
  data_c <- TDiff_PSIbin(data)
  data_c <- na.omit(data_c)
  value_col_c <- "avg_transpiration_deficit"
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
    
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Poly2", "Poly3"),
      AIC = c(aic_linear, aic_poly2, aic_poly3),
      R2 = c(r2_linear, r2_poly2, r2_poly3)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_c[[sp]] <- sp_results
  }
  aic_df_c <- do.call(rbind, aic_results_c)
  aic_df_c$species <- factor(aic_df_c$species, levels = species_order)
  
  model_palette_c <- c("Linear" = "orange",
                       "Poly2"  = "dodgerblue",
                       "Poly3"  = "green4")
  
  p_c <- ggplot(aic_df_c, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "transpiration deficit ~ soil water potential") +
    scale_fill_manual(values = model_palette_c, name = "") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background= element_rect(fill = "white", color = "white"),
      panel.background=element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title     = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      axis.title     = element_text(face = "bold", size = 16),
      axis.text      = element_text(color = "black", size = 14),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text    = element_text(size = 14)
    )+
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine Panels A and B with a shared legend placed below (centered)
  combined_top <- (p_a | p_b) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "top", legend.box = "horizontal")
  
  # Combine with Panel C below
  combined_plot <- combined_top / p_c
  
  # Save the combined figure to file
  dir.create(dirname(save_combined_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_combined_fig, plot = combined_plot, width = 14, height = 12, dpi = 300)
  
  print(combined_plot)
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
    labs(x = "Soil Water Potential (kPa)", y = "NDVI Quantiles (rank)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
      legend.position = "none",
      legend.text = element_text(size = 14) 
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))+
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential (kPa)") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
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
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.caption.position = "plot")
  
  #### COMBINE PANELS ####
  # Add legend to panel c only, remove from others
  p_combined <- p_combined + theme(legend.position = "none")
  p_x50 <- p_x50 + theme(legend.position = "bottom")
  p_slope <- p_slope + theme(legend.position = "none")  # Keep legend here
  
  # Combine plots and collect legend
  final_slope_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  
  # Save and print
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
    labs(title = "NDVI quantiles ~ soil water potential",
         subtitle = expression(italic(NDVI) == a + b * e^{c * italic(x)} + epsilon[i]),
         x = "",
         y = "Coefficient Value") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size =14),
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

plot_Quantiles_TDiff_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(nlme)    # for consistency (even though lm is used for linear)
  library(car)     # for deltaMethod
  library(broom)
  
  # Define the value column and species groups.
  # Use linear regression for Oak and Beech,
  # and the exponential model for Spruce and Pine.
  value_col <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species <- c("Spruce", "Pine")
  data$species <- factor(data$species, levels = c(linear_species, exp_species))
  
  # Define the color palette (same for all species)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, set x = bin_median
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  
  # Linear models for Oak and Beech
  models_linear <- list()
  for(sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    lm_model <- lm(avg_value ~ x, data = sp_data)
    models_linear[[sp]] <- lm_model
  }
  
  # Exponential models for Spruce and Pine
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    exp_model <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      NULL
    })
    models_exp[[sp]] <- exp_model
  }
  
  ##############################################
  # Coefficient Extraction and Plotting - Linear
  ##############################################
  
  lin_coef_list <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Rename coefficients for clarity using dplyr::recode
    df$Coefficient <- dplyr::recode(df$Coefficient, "(Intercept)" = "a", "x" = "b")
    return(df)
  })
  lin_coef_df <- do.call(rbind, lin_coef_list)
  lin_coef_df$species <- factor(lin_coef_df$species, levels = linear_species)
  lin_coef_df <- lin_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_linear <- ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Linear Model (Oak & Beech)",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * x + epsilon[i]),
         subtitle = expression(italic(NDVI) == a + b * x + epsilon[i]),
         x = "",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ##############################################
  # Coefficient Extraction and Plotting - Exponential
  ##############################################
  
  exp_coef_list <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Convert coefficient names to lowercase (a, b, c)
    df$Coefficient <- tolower(df$Coefficient)
    return(df)
  })
  exp_coef_df <- do.call(rbind, exp_coef_list)
  exp_coef_df$species <- factor(exp_coef_df$species, levels = exp_species)
  exp_coef_df <- exp_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Exponential Model (Spruce & Pine)",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * e^{-c * italic(x)} + epsilon[i]),
         subtitle = expression(italic(NDVI) == a + b * e^{-c * italic(x)} + epsilon[i]),
         x = "",
         y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine the coefficient plots into one figure with panel labels (a) and (b)
  combined_coeff <- p_coeff_linear + p_coeff_exp
  print(combined_coeff)
  
  # Save the combined coefficient plot
  ggsave(filename = combined_coef_fig, plot = combined_coeff, device = "png", width = 10, height = 8, dpi = 300)
  
  #################################################
  # Combined Final Plot with Model Predictions (Panel A)
  #################################################
  
  # Generate predictions for linear species (Oak & Beech)
  pred_linear <- lapply(linear_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    lm_mod <- models_linear[[sp]]
    pred <- predict(lm_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_linear_df <- do.call(rbind, pred_linear)
  
  # Generate predictions for exponential species (Spruce & Pine)
  pred_exp <- lapply(exp_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    nls_mod <- models_exp[[sp]]
    if(is.null(nls_mod)) return(NULL)
    pred <- predict(nls_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_exp_df <- do.call(rbind, pred_exp)
  
  # Combine all predictions for Panel A
  pred_all <- bind_rows(pred_linear_df, pred_exp_df)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -2.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = "Transpiration Deficit (mm)", y = "NDVI Quantiles (rank)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
      legend.position = "none",
      legend.text = element_text(size = 14) 
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5)) +
    coord_cartesian(clip = "off")
  
  #################################################
  # Additional Analysis for Panel B & C: x50 and Slope for All Species
  #################################################
  
  # For linear models: calculate x50 and absolute slope (and p-value from slope)
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]
    slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]
    p_val_lin <- summ["x", "Pr(>|t|)"]
    x50_lin <- (threshold - intercept) / slope
    abs_slope_lin <- abs(slope)
    r2_lin <- summary(mod)$r.squared
    data.frame(species = sp, model = "Linear", x50 = x50_lin, abs_slope = abs_slope_lin,
               r_squared = r2_lin, se = se_slope, p_val = p_val_lin)
  })
  
  # For exponential models: calculate x50 and absolute slope with deltaMethod and compute p-value
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    coefs <- coef(mod)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    c_val <- coefs["c"]
    x50_exp <- ifelse((threshold - a_val) > 0 & b_val > 0,
                      -log((threshold - a_val)/b_val) / c_val,
                      NA)
    slope50_exp <- -c_val * (threshold - a_val)
    abs_slope_exp <- abs(slope50_exp)
    sp_data <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = sp_data)
    r2_exp <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) / 
      sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se_exp <- dm_result$SE
    df_exp <- summary(mod)$df[2]
    t_val_exp <- slope50_exp / se_exp
    p_val_exp <- 2 * (1 - pt(abs(t_val_exp), df_exp))
    data.frame(species = sp, model = "Exponential", x50 = x50_exp, abs_slope = abs_slope_exp,
               r_squared = r2_exp, se = se_exp, p_val = p_val_exp)
  })
  
  stats_all <- rbind(do.call(rbind, lin_stats), do.call(rbind, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: Bar plot of x50 (Transpiration Deficit) for all species
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    # geom_text(aes(label = sprintf("%.2f", x50), y = x50/2),
    #          color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "Transpiration Deficit (mm)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
      legend.position = "none",
      legend.text = element_text(size = 14) 
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Panel C: Bar plot of absolute slope for all species with error bars and R² annotations.
  # The annotation includes an asterisk if p < 0.05.
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2),
              color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "Absolute Slope") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
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
      legend.position = "none",
      legend.text = element_text(size = 14) 
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),         # <- No background
      legend.box.background = element_blank()      # <- No box border
    )
  print(final_plot)
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_combined_AIC_R2(all_results_df,  "results/key_displays/AIC.png")

plot_time_series_and_correlation_combined(all_results_df,
                                          "results/key_displays/time_series_Quantiles_PSI_TDiff_species.png")

plot_yearly_mean_linear_coeffs(all_results_df, "results/key_displays/yearly_mean_linear_coeffs.png")

plot_NDVI_Q_PSIbin_exp_slope_negPSI(all_results_df, 
                                    "results/key_displays/NDVI_Q_PSIbin_exp_coeff_negPSI.png",
                                    "results/key_displays/NDVI_Q_PSIbin_exp_slope_negPSI.png")

plot_Quantiles_TDiff_exp_linear_slope_coeff(all_results_df,
                                            "results/key_displays/NDVI_Q_TDiffbin_exp_linear_coeff.png",
                                            "results/key_displays/NDVI_Q_TDiffbin_exp_linear_slope.png")


plot_TDiff_PSIbin_poly_2_slope <- function(data, coef_output, figure_output) {
  # Load required libraries (if not already loaded)
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(gridExtra) # For arranging multiple plots
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  
  # Process the data and remove any rows with missing values.
  TDiff_PSIbin_df <- TDiff_PSIbin(data)
  TDiff_PSIbin_df <- na.omit(TDiff_PSIbin_df)
  
  # Define a custom color palette and set species order.
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  TDiff_PSIbin_df$species <- factor(TDiff_PSIbin_df$species, levels = species_levels)
  
  ## Panel A: Mixed-Effects Model Plot
  
  # Fit the mixed-effects model with a second-order polynomial for bin_median.
  model <- lmer(avg_transpiration_deficit ~ poly(bin_median, 2) + 
                  (poly(bin_median, 2) | species),
                data = TDiff_PSIbin_df)
  
  # Create prediction data for each species.
  pred_data <- TDiff_PSIbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median)) %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Predict the fitted values using the model (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Keep only the segment of the fitted line starting from the maximum fitted value.
  pred_data <- pred_data %>%
    group_by(species) %>%
    filter(bin_median >= bin_median[which.max(predicted)]) %>%
    ungroup()
  
  # Calculate the overall mean and median of transpiration deficit.
  line_val <- mean(data$transpiration_deficit)
  line_median <- median(data$transpiration_deficit)
  
  # Create the mixed-effects model plot with horizontal lines for mean and median.
  plot_mixed <- ggplot(TDiff_PSIbin_df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), linewidth = 1) +
    geom_hline(yintercept = line_val, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_val, 
             label = paste0("mean: ", round(line_val, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_median, 
             label = paste0("median: ", round(line_median, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (bin_median)",
         y = "Average Transpiration Deficit") +
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
      df <- .
      max_pred <- max(df$predicted, na.rm = TRUE)
      thresh_val <- 0.3 * max_pred
      idx <- which.min(abs(df$predicted - thresh_val))
      win_start <- max(1, idx - 5)
      win_end <- min(nrow(df), idx + 5)
      df_window <- df[win_start:win_end, ]
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
  
  # Format the label for Panel C: if p < 0.05, append an asterisk to the R².
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
                                "poly(bin_median, 2)1" = "b",
                                "poly(bin_median, 2)2" = "c"))
  
  # Set the order of the term factor.
  coeff_data$term <- factor(coeff_data$term, levels = c("a", "b", "c"))
  
  ## Compute p-values for coefficients by fitting separate linear models for each species.
  species_list <- species_levels
  coeff_stats_list <- lapply(species_list, function(sp) {
    subdata <- subset(TDiff_PSIbin_df, species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 2), data = subdata)
    summ <- summary(mod_sp)$coefficients
    df <- as.data.frame(summ)
    df$term <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- do.call(rbind, coeff_stats_list)
  coeff_stats <- coeff_stats %>%
    mutate(term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 2)1" = "b",
                                "poly(bin_median, 2)2" = "c")) %>%
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
         subtitle = expression(hat(Y) == a + b*x + c*x^2),
         x = "Coefficient Term", 
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
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
      strip.text = element_text(face = "bold", size = 16)
    )
  
  print(plot_coeff)
  
  # Save the coefficient plot with width = 10 and height = 8.
  ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
  
  print(coeff_data)
  # Return both plots as a list (if desired)
  return(list(combined_plot = final_plot, coeff_plot = plot_coeff))
}


plot_TDiff_PSIbin_poly_2_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_2_coeff_each.png",
                               "results/key_displays/TDiff_PSIbin_poly_2_slope_each.png")