setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)

data <- combined %>% filter(month %in% c("July", "August"))
data <- data %>% filter(Quantiles > 0)
data <- na.omit(data)

### Figure 1 - time series with correlation of NDVI, PSI, TDiff ###
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
  
  # 1) Prepare species‐averaged data
  df_species <- df_average_year_species(df_all)
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  # 2) Color palette
  cb_palette <- c(
    Oak    = "#E69F00",
    Beech  = "#0072B2",
    Spruce = "#009E73",
    Pine   = "#F0E442"
  )
  
  # 3) Pick the NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion found in df_species.")
  }
  
  # 3b) Print correlation coefficients by species
  cat("Pearson correlations by species:\n")
  
  df_species %>%
    group_by(species) %>%
    summarise(
      cor_psi   = cor(.data[[ndvi_column]], avg_psi,   use = "complete.obs"),
      cor_tdiff = cor(.data[[ndvi_column]], avg_tdiff, use = "complete.obs")
    ) %>%
    ungroup() %>%
    print()
  
  # 4) Common theme: transparent backgrounds, no grid
  base_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      legend.text = element_text(color = "black", size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # 5) x‐axis scale: breaks every 2 years
  x_scale <- scale_x_continuous(breaks = seq(2003, 2024, by = 2),
                                limits = c(2003, 2024))
  
  ### TIME‐SERIES PANELS (a), (b), (c)
  
  p1 <- ggplot(df_species,
               aes(x = as.numeric(year),
                   y = .data[[ndvi_column]],
                   color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(
      x = "year",
      y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))),
      color = ""
    ) +
    scale_color_manual(values = cb_palette) +
    x_scale + base_theme +
    ggtitle("(a)") +
    theme(
      plot.title.position = "panel",  # place title below collected legend
      plot.title          = element_text(
        face   = "bold",
        size   = 18,
        hjust  = 0.0,
        vjust  = -9,
        margin = margin(b = 5)
      )
    )
  
  p2 <- ggplot(df_species,
               aes(x = as.numeric(year),
                   y = avg_psi,
                   color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(
      x = "year",
      y = expression(atop(bold("soil water potential"), bold("(kPa)")))
    ) +
    scale_color_manual(values = cb_palette) +
    x_scale + base_theme +
    guides(color = "none") +
    ggtitle("(b)") +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.1,
        vjust = 1
      )
    )
  
  p3 <- ggplot(df_species,
               aes(x = as.numeric(year),
                   y = avg_tdiff,
                   color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(
      x = "year",
      y = expression(atop(bold("transpiration deficit"), bold("(mm)")))
    ) +
    scale_color_manual(values = cb_palette) +
    x_scale + base_theme +
    guides(color = "none") +
    ggtitle("(c)") +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.1,
        vjust = 1
      )
    )
  
  ts_plot <- p1 / p2 / p3
  
  ### CORRELATION PANELS (d), (e)
  
  p4 <- ggplot(df_species,
               aes(x = avg_psi,
                   y = .data[[ndvi_column]],
                   color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "soil water potential (kPa)",
      y = expression(atop(bold("NDVI quantiles"), bold("(rank)")))
    ) +
    scale_color_manual(values = cb_palette) +
    base_theme +
    guides(color = "none") +
    ggtitle("(d)") +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.2,
        vjust = 1
      )
    )
  
  p5 <- ggplot(df_species,
               aes(x = avg_tdiff,
                   y = .data[[ndvi_column]],
                   color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      x = "transpiration deficit (mm)",
      y = NULL,
      color = "Species"
    ) +
    scale_color_manual(values = cb_palette) +
    base_theme +
    guides(color = "none") +
    ggtitle("(e)") +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.06,
        vjust = 1
      )
    )
  
  corr_plot <- p4 | p5
  
  # 6) Combine everything, collect one legend at top
  final_plot <- (ts_plot / corr_plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  # 7) Draw and save
  print(final_plot)
  ggsave(output_path,
         plot   = final_plot,
         width  = 12,
         height = 14,
         dpi    = 300)
}

plot_time_series_and_correlation_combined(combined,
                                          "results_rootzone/Figures/time_series_Quantiles_PSI_TDiff_species_rootzone.png")

### Figure 2 - NDVI - PSI ###
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

plot_NDVI_PSI_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure, aic_barplot_fig) {
  # Process data using the custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_PSIbin(data)
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
  linear_species <- c("Beech")
  exp_species <- c("Oak", "Spruce", "Pine")
  species_all <- cbind(linear_species, exp_species)
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Define the color palette (same for all species)
  cb_palette <- c("Oak" = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For transpiration deficit, set x = bin_median
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
  # Fit linear models and exponential models for each species
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  aic_list <- lapply(species_all, function(sp) {
    df <- data_clean %>% filter(species == sp)
    # Linear fit
    lm_mod <- tryCatch(lm(avg_value ~ x, data = df), error = function(e) NULL)
    aic_l  <- if (!is.null(lm_mod)) AIC(lm_mod) else NA
    # Exponential fit
    exp_mod <- tryCatch(
      nls(avg_value ~ a + b * exp(c * x), data = df,
          start = start_list, control = ctrl),
      error = function(e) NULL
    )
    aic_e <- if (!is.null(exp_mod)) AIC(exp_mod) else NA
    tibble(species = sp,
           AIC_linear      = aic_l,
           AIC_exponential = aic_e)
  })
  aic_df <- bind_rows(aic_list)
  message("AIC values for each species:")
  print(aic_df)
  
  # Compute AIC for each model
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  aic_list <- lapply(c(linear_species, exp_species), function(sp) {
    df_sp <- data_clean %>% filter(species == sp)
    lm_mod  <- tryCatch(lm(avg_value ~ x, data = df_sp), error = function(e) NULL)
    exp_mod <- tryCatch(nls(avg_value ~ a + b * exp(c * x), data = df_sp, start = start_list, control = ctrl), error = function(e) NULL)
    tibble(
      species = sp,
      AIC_linear      = if (!is.null(lm_mod)) AIC(lm_mod) else NA,
      AIC_exponential = if (!is.null(exp_mod)) AIC(exp_mod) else NA
    )
  })
  aic_df <- bind_rows(aic_list)
  
  # Bar plot of AIC values for each species
  aic_long <- aic_df %>%
    pivot_longer(cols = c(AIC_linear, AIC_exponential), names_to = "model", values_to = "AIC") %>%
    mutate(
      # use dplyr::recode to avoid masking by car::recode
      model = dplyr::recode(model,
                            AIC_linear = "linear",
                            AIC_exponential = "exponential")
    )
  
  aic_long$species   <- factor(aic_long$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  p_aic <- ggplot(aic_long, aes(x = species, y = AIC, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("linear" = "dodgerblue", "exponential" = "orange"), name = "") +
    labs(title = "", x = "", y = "AIC") +
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
    )
  
  print(p_aic)
  ggsave(filename = aic_barplot_fig, plot = p_aic, device = "png", width = 8, height = 6, dpi = 300)
  
  
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
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Check that sp_data has enough data
    if(nrow(sp_data) < 5) {
      message("Not enough data for species ", sp)
      models_exp[[sp]] <- NULL
      next
    }
    
    exp_model <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x),
          data = sp_data,
          start = start_list,
          control = control_params)
    }, error = function(e) {
      message("Error for species ", sp, ": ", e$message)
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
    geom_bar(
      stat     = "identity",
      position = position_dodge(width = 0.8),  # how far apart the centers sit
      width    = 0.5                            # bar thickness
    ) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models",
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "",
         y = "coefficient value") +
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
      legend.position = "none",
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
    if (is.null(mod)) return(NULL)
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
  # Remove any NULL entries:
  exp_coef_list <- exp_coef_list[!sapply(exp_coef_list, is.null)]
  exp_coef_df <- do.call(rbind, exp_coef_list)
  exp_coef_df$species <- factor(exp_coef_df$species, levels = exp_species)
  exp_coef_df <- exp_coef_df %>%
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(
      stat     = "identity",
      position = position_dodge(width = 0.8),  # how far apart the centers sit
      width    = 0.8                            # bar thickness
    ) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.8),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models",
         subtitle = expression(NDVI == a + b * e^{c * italic(x)} + epsilon),
         x = "",
         y = "coefficient value") +
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
      legend.position = "none",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0.5))
  
  # Combine the coefficient plots into one figure with panel labels (a) and (b)
  # Adjusted spacing between facets for better readability
  combined_coeff <- p_coeff_linear + p_coeff_exp +
    plot_layout(widths = c(0.6, 1.4), guides = "collect") &
    theme(
      plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"), # Add some margin around plots
      strip.placement = "outside" # Place strip labels outside the plot area
    )
  print(combined_coeff)
  
  # Save the combined coefficient plot
  ggsave(filename = combined_coef_fig, plot = combined_coeff, device = "png", width = 14, height = 8, dpi = 300) # Increased width
  
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
    # Incorporate point aesthetics from plot_Quantiles_PSI_exp_linear_slope_coeff_2
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species, shape = species, size = percentage), alpha = 0.7) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") + # Use "Species" as legend title
    # Add shape manual scale
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18), name = "Species", guide = FALSE) +
    # Add size continuous scale with percentage labels
    scale_size_continuous(name = "Pixels per bin (%)", range = c(1, 8), labels = scales::percent_format(accuracy = 1)) +
    # Arrange guides to have color/shape first, then size
    guides(
      color = guide_legend(order = 1),
      # shape = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      # position the title at the very top-left of the entire plot
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.1,    # left-align against the y-axis
        vjust = 1     # push into top margin
      ),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "bottom",
      legend.text         = element_text(size = 14),
      legend.title        = element_text(size = 14, face = "bold"),
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
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
    data.frame(species = sp, model = "linear", x50 = x50_lin, abs_slope = abs_slope_lin,
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
                      log((threshold - a_val)/b_val) / c_val,
                      NA)
    slope50_exp <- c_val * (threshold - a_val)
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
    data.frame(species = sp, model = "exponential", x50 = x50_exp, abs_slope = abs_slope_exp,
               r_squared = r2_exp, se = se_exp, p_val = p_val_exp)
  })
  
  stats_all <- rbind(do.call(rbind, lin_stats), do.call(rbind, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Panel B: Bar plot of x50 (transpiration deficit) for all species
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      # position the title at the very top-left of the entire plot
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.25,    # left-align against the y-axis
        vjust = 1     # push into top margin
      ),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "bottom",
      legend.text         = element_text(size = 14),
      legend.title        = element_text(size = 14, face = "bold"),
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # Panel C: Bar plot of absolute slope for all species with error bars and R² annotations.
  # The annotation includes an asterisk if p < 0.05.
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
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
      # position the title at the very top-left of the entire plot
      plot.title.position = "plot",
      plot.title          = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.25,    # left-align against the y-axis
        vjust = 1     # push into top margin
      ),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "bottom",
      legend.text         = element_text(size = 14),
      legend.title        = element_text(size = 14, face = "bold"),
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "bottom", # Moved to bottom as in _2 function
      legend.title = element_blank(), # Kept blank as in _2 function, but can be customized
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  print(final_plot)
  
  #### Plateau point calculation for exponential models (Spruce & Pine only)
  plateau_df <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    coefs <- coef(mod)
    c_val <- coefs["c"]
    x_plateau <- -log(0.05) / c_val
    data.frame(species = sp, x_plateau = x_plateau)
  }) %>% bind_rows()
  
  print("Estimated soil water potential (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_PSI_exp_linear_slope_coeff(data, 
                                     "results_rootzone/Figures/NDVI_Q_PSIbin_exp_linear_coeff.png",
                                     "results_rootzone/Figures/NDVI_Q_PSIbin_exp_linear_slope.png",
                                     "results_rootzone/Figures/NDVI_Q_PSIbin_exp_linear_aic.png")
### Figure 2 - NDVI - TDiff ###
NDVI_TDiffbin <- function(df, bin_width = 3) {
  
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
  
  library(dplyr); library(tibble)
  
  df <- as_tibble(df)
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

plot_NDVI_TDiff_exp_slope_coeff <- function(data, combined_coef_fig, output_figure, aic_barplot_fig) {
  # Process data and remove missing rows
  data <- NDVI_TDiffbin(data) %>% na.omit()
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(nlme)
  library(car)
  library(broom)
  
  # Define species and palette
  value_col   <- "avg_value"
  exp_species <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = exp_species)
  cb_palette  <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # Prepare data
  data_clean <- data %>%
    mutate(x = bin_median) %>%
    filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  # Define species and palette
  value_col   <- "avg_value"
  exp_species <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = exp_species)
  cb_palette  <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # Prepare data
  data_clean <- data %>%
    mutate(x = bin_median) %>%
    filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  ##########################
  # AIC Comparison         #
  ##########################
  start_nls <- list(a = 5, b = 7, c = 0.04)
  ctrl      <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  aic_list <- lapply(exp_species, function(sp) {
    df_sp <- filter(data_clean, species == sp)
    lm_mod <- tryCatch(lm(avg_value ~ x, data = df_sp), error = function(e) NULL)
    aic_l  <- if (!is.null(lm_mod)) AIC(lm_mod) else NA
    exp_mod <- tryCatch(
      nls(avg_value ~ a + b * exp(-c * x), data = df_sp,
          start = start_nls, control = ctrl),
      error = function(e) NULL
    )
    aic_e <- if (!is.null(exp_mod)) AIC(exp_mod) else NA
    tibble(species = sp,
           AIC_linear = aic_l,
           AIC_exponential = aic_e)
  }) %>% bind_rows()
  
  aic_long <- aic_list %>%
    pivot_longer(cols = c(AIC_linear, AIC_exponential),
                 names_to = "model", values_to = "AIC") %>%
    mutate(model = dplyr::recode(model,
                                 AIC_linear = "linear",
                                 AIC_exponential = "exponential"))
  aic_long$species <- factor(aic_long$species, levels = exp_species)
  
  p_aic <- ggplot(aic_long, aes(x = species, y = AIC, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c(linear = "dodgerblue", exponential = "orange"), name = "") +
    labs(x = "", y = "AIC") +
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
    )
  print(p_aic)
  ggsave(filename = aic_barplot_fig, plot = p_aic, width = 8, height = 6, dpi = 300)
  
  # Model fitting parameters
  start_list <- list(a = 5, b = 7, c = 0.04)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  # Fit exponential models for all species
  models_exp <- set_names(
    map(exp_species, function(sp) {
      sp_data <- filter(data_clean, species == sp)
      tryCatch(
        nls(avg_value ~ a + b * exp(-c * x), data = sp_data,
            start = start_list, control = ctrl),
        error = function(e) NULL
      )
    }),
    exp_species
  )
  
  # Extract exponential coefficients
  exp_coef_df <- map_dfr(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    tibble(
      species     = sp,
      Coefficient = tolower(rownames(summ)),
      Value       = summ[, "Estimate"],
      pvalue      = summ[, "Pr(>|t|)"]
    )
  }) %>%
    mutate(
      species = factor(species, levels = exp_species),
      label   = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue))
    )
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_col(position = position_dodge(0.9)) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(0.9), vjust = 0.5) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models",
         subtitle = expression(NDVI == a + b * e^{-~c * x} + epsilon),
         x = NULL, y = "coefficient value", caption = "") +
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
      legend.position = "none",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) 
  
  print(p_coeff_exp)
  
  # Generate predictions for all species
  pred_exp <- map_dfr(exp_species, function(sp) {
    sp_data <- filter(data_clean, species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    pred <- predict(mod, newdata = data.frame(x = x_seq))
    tibble(species = sp, x = x_seq, pred = pred)
  })
  
  # Combined final plot (Panel A)
  p_combined <- ggplot() +
    geom_point(data = data_clean,
               aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
               alpha = 0.7) +
    geom_line(data = pred_exp,
              aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 1) +
    annotate("text", x = 32, y = threshold,
             label = "median", fontface = "italic",
             hjust = 0.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette) +
    scale_shape_manual(values = c("Oak"=16,"Beech"=17,"Spruce"=15,"Pine"=18)) +
    scale_size_continuous(name = "Pixels per bin (%)",
                          range = c(1, 8),
                          labels = scales::percent_format(accuracy = 1)) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2)) +
    labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(color = "black", angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
    )
  
  # Stats for x50 and slope for all species
  stats_all <- map_dfr(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    coefs <- coef(mod)
    a_val <- coefs["a"]; b_val <- coefs["b"]; c_val <- coefs["c"]
    x50 <- ifelse((threshold - a_val)>0 & b_val>0, -log((threshold - a_val)/b_val)/c_val, NA)
    slope50 <- -c_val * (threshold - a_val)
    abs_slope <- abs(slope50)
    sp_data <- filter(data_clean, species == sp)
    fitted_vals <- predict(mod, newdata = sp_data)
    r2 <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) / sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    dm <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"), parameterNames=c("a","b","c"))
    se <- dm$SE; df_exp <- summary(mod)$df[2]; t_val <- slope50/se; p_val <- 2*(1-pt(abs(t_val), df_exp))
    tibble(species=sp, x50=x50, abs_slope=abs_slope, r_squared=r2, se=se, p_val=p_val)
  }) %>%
    mutate(species=factor(species, levels=exp_species))
  
  # Panel B: x50
  p_x50 <- ggplot(stats_all, aes(x=species, y=x50, fill=species)) +
    geom_col(width=0.7) + scale_fill_manual(values=cb_palette, guide=FALSE) +
    labs(y="transpiration deficit (mm)", x="") + ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(color = "black", angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
    )
  # Panel C: slope
  p_slope <- ggplot(stats_all, aes(x=species, y=abs_slope, fill=species)) +
    geom_col(width=0.7) +
    geom_errorbar(aes(ymin=pmax(0, abs_slope-se), ymax=abs_slope+se), width=0.2) +
    geom_text(aes(label=if_else(p_val<0.05, sprintf("%.2f*", r_squared), sprintf("%.2f", r_squared)),
                  y=abs_slope/2), size=5) +
    scale_fill_manual(values=cb_palette, guide=FALSE) + labs(y="absolute slope", x="") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(color = "black", angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
    )
  # Combine panels
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths=c(2,1), guides="collect") &
    theme(legend.position="bottom", legend.title=element_blank())
  print(final_plot)
  
  #### Plateau point calculation for exponential models (Spruce & Pine only)
  plateau_df <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    coefs <- coef(mod)
    c_val <- coefs["c"]
    x_plateau <- -log(0.05) / c_val
    data.frame(species = sp, x_plateau = x_plateau)
  }) %>% bind_rows()
  
  print("Estimated soil water potential (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  # Save outputs
  ggsave(combined_coef_fig, p_coeff_exp, width=10, height=8, dpi=300)
  ggsave(output_figure, final_plot, width=10, height=8, dpi=300)
}

plot_NDVI_TDiff_exp_slope_coeff(data,
                                "results_rootzone/Figures/NDVI_Q_TDiffbin_exp_coeff.png",
                                "results_rootzone/Figures/NDVI_Q_TDiffbin_exp_slope.png",
                                "results_rootzone/Figures/NDVI_Q_TDiffbin_exp_aic.png")
