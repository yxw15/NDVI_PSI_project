setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)

data <- combined %>% filter(month %in% c("July", "August"))
data <- data %>% filter(Quantiles > 0 & year < 2023)
data <- na.omit(data)

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
      avg_value  = mean(Quantiles, na.rm = TRUE),
      bin_median = median(soil_water_potential, na.rm = TRUE),  # x-value for plotting
      bin_mean   = mean(soil_water_potential, na.rm = TRUE),    # optional, keep if useful
      count      = n(),
      .groups    = "drop"
    ) %>%
    left_join(species_totals, by = "species") %>%
    filter(count >= 1000) %>%
    mutate(percentage = count / total_pixels) 
  return(meanNDVI_PSIbin_species)
}

## geom_smooth: linear + exponential
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
  data <- data %>% mutate(x = bin_mean)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
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
  # Combined Final Plot with Model Predictions + Confidence Bands (Panel A)
  #################################################
  
  p_combined <- ggplot() +
    
    # Raw binned observations
    geom_point(
      data = data_clean,
      aes(x = x, y = avg_value,
          color = species, shape = species, size = percentage),
      alpha = 0.7
    ) +
    
    # Linear species (Beech)
    geom_smooth(
      data = data_clean %>% filter(species %in% linear_species),
      aes(x = x, y = avg_value,
          color = species, fill = species),
      method = "lm",
      se = TRUE,            # <-- confidence band
      alpha = 0.25,
      linewidth = 1.2
    ) +
    
    # Exponential species (Oak, Spruce, Pine)
    geom_smooth(
      data = data_clean %>% filter(species %in% exp_species),
      aes(x = x, y = avg_value,
          color = species, fill = species),
      method = "nls",
      formula = y ~ a + b * exp(c * x),
      method.args = list(
        start = list(a = 5, b = 3, c = 0.001),
        control = nls.control(maxiter = 1200, minFactor = 1e-9)
      ),
      se = TRUE,            # <-- confidence band
      alpha = 0.25,
      linewidth = 1.2
    ) +
    
    # Threshold line
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    
    annotate(
      "text",
      x = min(data_clean$x, na.rm = TRUE),
      y = threshold,
      label = "median",
      hjust = -0.1,
      vjust = -0.3,
      fontface = "italic",
      size = 6
    ) +
    
    # Same palette for points, lines and ribbons
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, name = "") +
    
    # Shapes
    scale_shape_manual(
      values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
      guide = FALSE
    ) +
    
    # Point size = percentage
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    # Guides: avoid duplicate legend for ribbons
    guides(
      color = guide_legend(order = 1),
      fill  = "none",
      size  = guide_legend(order = 2)
    ) +
    
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
    ) +
    
    ggtitle("(a)") +
    
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.1,
        vjust = 1
      ),
      axis.text.x      = element_text(angle = 0, hjust = 0.5),
      axis.title       = element_text(face = "bold", size = 16),
      axis.text        = element_text(color = "black", size = 14),
      legend.position  = "bottom",
      legend.text      = element_text(size = 14),
      legend.title     = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid       = element_blank(),
      panel.border     = element_blank()
    )
  
  print(p_combined)
  
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

## geom_smooth: linear + gam 
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
  data <- data %>% mutate(x = bin_mean)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
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
  # Combined Final Plot with Model Predictions + Confidence Bands (Panel A)
  #################################################
  
  p_combined <- ggplot() +
    
    geom_point(
      data = data_clean,
      aes(x = x, y = avg_value,
          color = species, shape = species, size = percentage),
      alpha = 0.7
    ) +
    
    # Linear species (Beech)
    geom_smooth(
      data = data_clean %>% filter(species %in% linear_species),
      aes(x = x, y = avg_value, color = species, fill = species),
      method = "lm",
      se = TRUE,
      alpha = 0.25,
      linewidth = 1.2
    ) +
    
    # Nonlinear species (Oak, Spruce, Pine)
    geom_smooth(
      data = data_clean %>% filter(species %in% exp_species),
      aes(x = x, y = avg_value, color = species, fill = species),
      method = "gam",
      formula = y ~ s(x, k = 5),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1.2
    ) +
  
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    
    # Threshold line
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      color = "black",
      linewidth = 1
    ) +
    
    annotate(
      "text",
      x = min(data_clean$x, na.rm = TRUE),
      y = threshold,
      label = "median",
      hjust = -0.1,
      vjust = -0.3,
      fontface = "italic",
      size = 6
    ) +
    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, name = "") +
    
    scale_shape_manual(
      values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
      guide = FALSE
    ) +
    
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    guides(
      color = guide_legend(order = 1),
      fill  = "none",
      size  = guide_legend(order = 2)
    ) +
    
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
    ) +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        face  = "bold",
        size  = 18,
        hjust = 0.1,
        vjust = 1
      ),
      axis.text.x      = element_text(angle = 0, hjust = 0.5),
      axis.title       = element_text(face = "bold", size = 16),
      axis.text        = element_text(color = "black", size = 14),
      legend.position  = "bottom",
      legend.text      = element_text(size = 14),
      # legend.title     = element_text(size = 14, face = "bold"),
      legend.title     = element_text(size = 14),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid       = element_blank(),
      panel.border     = element_blank()
    )
  
   print(p_combined)
  

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

## MASS (bootstrapping 1000): linear + exponential
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
  library(MASS)
  
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
  data <- data %>% mutate(x = bin_mean)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
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
  
  pred_linear_ci <- lapply(linear_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE),
                 max(sp_data$x, na.rm = TRUE), length.out = 100)
    
    lm_mod <- models_linear[[sp]]
    pr <- predict(lm_mod, newdata = data.frame(x = x_seq), se.fit = TRUE)
    
    crit <- qt(0.975, df = df.residual(lm_mod))
    
    data.frame(
      species = sp,
      x = x_seq,
      fit = pr$fit,
      lwr = pr$fit - crit * pr$se.fit,
      upr = pr$fit + crit * pr$se.fit
    )
  })
  pred_linear_ci_df <- do.call(rbind, pred_linear_ci)
  
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

  pred_exp_ci <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)

    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE),
                 max(sp_data$x, na.rm = TRUE), length.out = 100)

    # Fitted coefficients and covariance
    beta_hat <- coef(mod)
    V <- vcov(mod)

    # Bootstrap
    nboot <- 1000
    beta_sim <- MASS::mvrnorm(nboot, beta_hat, V)

    # Prediction function
    pred_fun <- function(b, x) {
      a <- b[1]; b2 <- b[2]; c <- b[3]
      a + b2 * exp(c * x)
    }

    # Simulated predictions
    y_sim <- sapply(1:nboot, function(i) pred_fun(beta_sim[i, ], x_seq))
    # y_sim: rows = x, cols = bootstrap runs

    data.frame(
      species = sp,
      x = x_seq,
      fit = predict(mod, newdata = data.frame(x = x_seq)),
      lwr = apply(y_sim, 1, quantile, probs = 0.025, na.rm = TRUE),
      upr = apply(y_sim, 1, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
  pred_exp_ci_df <- bind_rows(pred_exp_ci)
  
  pred_exp_ci_df <- do.call(rbind, pred_exp_ci)
  
  pred_ci_all <- bind_rows(pred_linear_ci_df, pred_exp_ci_df)
  
  # Combine all predictions for Panel A
  pred_all <- bind_rows(pred_linear_df, pred_exp_df)
  
  p_combined <- ggplot() +
    # Incorporate point aesthetics from plot_Quantiles_PSI_exp_linear_slope_coeff_2
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species, shape = species, size = percentage), alpha = 0.7) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    geom_ribbon(
      data = pred_ci_all,
      aes(x = x, ymin = lwr, ymax = upr, fill = species),
      alpha = 0.2,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_fill_manual(values = cb_palette, name = "") +
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
  
  print(p_combined)
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

## MASS (bootstrapping 2000): linear + exponential  
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
  library(nlme)    
  library(car)     
  library(broom)
  library(MASS)
  library(scales) # Required for percent_format
  
  value_col <- "avg_value"
  linear_species <- c("Beech")
  exp_species <- c("Oak", "Spruce", "Pine")
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce"= "#009E73", "Pine"  = "#F0E442")
  
  data <- data %>% mutate(x = bin_mean)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  ##########################
  # AIC Calculation        #
  ##########################
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  aic_list <- lapply(unique(data_clean$species), function(sp) {
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
  
  aic_long <- aic_df %>%
    pivot_longer(cols = c(AIC_linear, AIC_exponential), names_to = "model", values_to = "AIC") %>%
    mutate(model = dplyr::recode(model, AIC_linear = "linear", AIC_exponential = "exponential"))
  
  p_aic <- ggplot(aic_long, aes(x = species, y = AIC, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("linear" = "dodgerblue", "exponential" = "orange"), name = "") +
    theme_minimal() + 
    theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA))
  
  ggsave(filename = aic_barplot_fig, plot = p_aic, width = 8, height = 6, dpi = 300)
  
  ##########################
  # Model Fitting          #
  ##########################
  models_linear <- list()
  for(sp in linear_species) {
    models_linear[[sp]] <- lm(avg_value ~ x, data = data_clean %>% filter(species == sp))
  }
  
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    if(nrow(sp_data) < 5) next
    models_exp[[sp]] <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), data = sp_data, start = start_list, control = ctrl)
    }, error = function(e) NULL)
  }
  
  ##############################################
  # Predictions and Confidence Intervals       #
  ##############################################
  
  # Linear CIs
  pred_linear_ci <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    x_seq <- seq(min(data_clean$x[data_clean$species==sp]), max(data_clean$x[data_clean$species==sp]), length.out = 100)
    pr <- predict(mod, newdata = data.frame(x = x_seq), se.fit = TRUE)
    crit <- qt(0.975, df = df.residual(mod))
    data.frame(species = sp, x = x_seq, fit = pr$fit, lwr = pr$fit - crit * pr$se.fit, upr = pr$fit + crit * pr$se.fit)
  })
  
  # Exponential CIs (Improved Monte Carlo Approach)
  pred_exp_ci <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    
    x_seq <- seq(min(data_clean$x[data_clean$species==sp]), max(data_clean$x[data_clean$species==sp]), length.out = 100)
    beta_hat <- coef(mod)
    V <- vcov(mod)
    
    nboot <- 2000
    beta_sim <- MASS::mvrnorm(nboot, beta_hat, V)
    
    pred_fun <- function(b, x) b[1] + b[2] * exp(b[3] * x)
    
    # Calculate predictions for every bootstrap sample across the x_seq
    y_sim <- apply(beta_sim, 1, function(b) pred_fun(b, x_seq)) 
    
    data.frame(
      species = sp,
      x = x_seq,
      fit = predict(mod, newdata = data.frame(x = x_seq)),
      lwr = apply(y_sim, 1, quantile, probs = 0.025, na.rm = TRUE),
      upr = apply(y_sim, 1, quantile, probs = 0.975, na.rm = TRUE)
    )
  })
  
  pred_ci_all <- bind_rows(c(pred_linear_ci, pred_exp_ci))
  
  #################################################
  # Panel A: Main Plot                            #
  #################################################
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species, shape = species, size = percentage), alpha = 0.5) +
    geom_ribbon(data = pred_ci_all, aes(x = x, ymin = lwr, ymax = upr, fill = species), alpha = 0.2, show.legend = FALSE) +
    geom_line(data = pred_ci_all, aes(x = x, y = fit, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed") +
    scale_fill_manual(values = cb_palette) +
    scale_color_manual(values = cb_palette) +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18)) +
    scale_size_continuous(range = c(1, 8), labels = scales::percent_format()) +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)", title = "(a)") +
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
  print(p_combined)
  
  #################################################
  # Panels B & C: x50 and Slopes                  #
  #################################################
  
  # Calculation Logic
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    summ <- summary(mod)$coefficients
    intercept <- summ[1,1]; slope <- summ[2,1]
    data.frame(species = sp, x50 = (threshold - intercept)/slope, abs_slope = abs(slope), 
               se = summ[2,2], p_val = summ[2,4], r_squared = summary(mod)$r.squared)
  })
  
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    cc <- coef(mod)
    x50_v <- log((threshold - cc["a"])/cc["b"]) / cc["c"]
    slope50 <- cc["c"] * (threshold - cc["a"])
    dm <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"))
    
    fitted_vals <- predict(mod)
    sp_data <- data_clean %>% filter(species == sp)
    r2 <- 1 - sum(residuals(mod)^2)/sum((sp_data$avg_value - mean(sp_data$avg_value))^2)
    
    data.frame(species = sp, x50 = as.numeric(x50_v), abs_slope = abs(slope50), 
               se = dm$SE, p_val = 2*(1-pt(abs(slope50/dm$SE), summary(mod)$df[2])), r_squared = r2)
  })
  
  stats_all <- bind_rows(c(lin_stats, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = levels(data$species))
  
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(y = "soil water potential (kPa)", title = "(b)", x="") + theme_minimal() + theme(legend.position="none", panel.grid = element_blank())
  
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")), y = abs_slope/2)) +
    scale_fill_manual(values = cb_palette) +
    labs(y = "absolute slope", title = "(c)", x="") + theme_minimal() + theme(legend.position="none", panel.grid = element_blank())
  
  # Final Layout
  final_plot <- (p_combined + (p_x50 / p_slope)) + plot_layout(widths = c(2, 1), guides = "collect") & theme(legend.position = "bottom")
  
  ggsave(filename = output_figure, plot = final_plot, width = 12, height = 8, dpi = 300)
  
  # Plateau point calculation
  plateau_df <- map_df(models_exp, function(mod) {
    if(is.null(mod)) return(NULL)
    data.frame(x_plateau = -log(0.05) / coef(mod)["c"])
  }, .id = "species")
  
  print(plateau_df)
}

# Delta method
plot_NDVI_PSI_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure, aic_barplot_fig) {
  # 1. Data Preparation
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(nlme)    
  library(car)     
  library(MASS)
  library(scales) 
  
  value_col <- "avg_value"
  linear_species <- c("Beech")
  exp_species <- c("Oak", "Spruce", "Pine")
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce"= "#009E73", "Pine"  = "#F0E442")
  data <- data %>% mutate(x = bin_mean)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  # 2. Model Fitting
  start_list <- list(a = 5, b = 3, c = 0.001)
  ctrl       <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  models_linear <- list()
  for(sp in linear_species) {
    models_linear[[sp]] <- lm(avg_value ~ x, data = data_clean %>% filter(species == sp))
  }
  
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    models_exp[[sp]] <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), data = sp_data, start = start_list, control = ctrl)
    }, error = function(e) NULL)
  }
  
  # 3. Prediction and SYMMETRICAL Confidence Intervals
  # Linear CIs (Standard)
  pred_linear_ci <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    x_seq <- seq(min(data_clean$x[data_clean$species==sp]), max(data_clean$x[data_clean$species==sp]), length.out = 100)
    pr <- predict(mod, newdata = data.frame(x = x_seq), se.fit = TRUE)
    crit <- qt(0.975, df = df.residual(mod))
    data.frame(species = sp, x = x_seq, fit = pr$fit, lwr = pr$fit - crit * pr$se.fit, upr = pr$fit + crit * pr$se.fit)
  })
  
  # Exponential CIs (Using Delta Method logic to ensure symmetry)
  pred_exp_ci <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    
    x_seq <- seq(min(data_clean$x[data_clean$species==sp]), max(data_clean$x[data_clean$species==sp]), length.out = 100)
    
    # We use the deltaMethod to get SE for each point on the curve
    # This ensures the ribbon is centered on the 'fit'
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
  
  # 4. Plotting Panel A
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species, shape = species, size = percentage), alpha = 0.4) +
    geom_ribbon(data = pred_ci_all, aes(x = x, ymin = lwr, ymax = upr, fill = species), alpha = 0.25) +
    geom_line(data = pred_ci_all, aes(x = x, y = fit, color = species), linewidth = 1.2) +
    geom_hline(yintercept = threshold, linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = cb_palette) +
    scale_color_manual(values = cb_palette) +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18)) +
    scale_size_continuous(range = c(1, 8), labels = scales::percent_format()) +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)", title = "(a)") +
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
  print(p_combined)
  # 5. Stats Analysis for B & C
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    summ <- summary(mod)$coefficients
    data.frame(species = sp, x50 = (threshold - summ[1,1])/summ[2,1], abs_slope = abs(summ[2,1]), 
               se = summ[2,2], p_val = summ[2,4], r_squared = summary(mod)$r.squared)
  })
  
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    cc <- coef(mod)
    x50_v <- log((threshold - cc["a"])/cc["b"]) / cc["c"]
    slope50 <- cc["c"] * (threshold - cc["a"])
    dm <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"))
    r2 <- 1 - sum(residuals(mod)^2)/sum((data_clean$avg_value[data_clean$species==sp] - mean(data_clean$avg_value[data_clean$species==sp]))^2)
    data.frame(species = sp, x50 = as.numeric(x50_v), abs_slope = abs(slope50), 
               se = dm$SE, p_val = 2*(1-pt(abs(slope50/dm$SE), summary(mod)$df[2])), r_squared = r2)
  })
  
  stats_all <- bind_rows(c(lin_stats, exp_stats))
  
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) + scale_fill_manual(values = cb_palette) +
    labs(y = "soil water potential (kPa)", title = "(b)", x="") + theme_minimal() + theme(legend.position="none")
  
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")), y = abs_slope/2)) +
    scale_fill_manual(values = cb_palette) +
    labs(y = "absolute slope", title = "(c)", x="") + theme_minimal() + theme(legend.position="none")
  
  # Final Layout and Save
  final_plot <- (p_combined + (p_x50 / p_slope)) + plot_layout(widths = c(1.5, 1), guides = "collect") & theme(legend.position = "bottom")
  
  ggsave(filename = output_figure, plot = final_plot, width = 12, height = 8, dpi = 300)
  
  print("Process Complete. Symmetrical ribbons applied via Delta Method.")
}

plot_NDVI_PSI_exp_linear_slope_coeff(data, 
                                     "results_rootzone/test/NDVI_Q_PSIbin_exp_linear_coeff.png",
                                     "results_rootzone/test/NDVI_Q_PSIbin_exp_linear_slope.png",
                                     "results_rootzone/test/NDVI_Q_PSIbin_exp_linear_aic.png")
