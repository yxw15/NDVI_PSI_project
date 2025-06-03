library(ncdf4)
library(reshape2)
library(ggplot2)
library(terra)
library(dplyr)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

df <- final_df %>% filter(month %in% c("July", "August"))

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

plot_Quantiles_PSI_exp_linear_slope_coeff <- function(data_input, combined_coef_fig, output_figure) {
  # Process data
  data <- NDVI_PSIbin(data_input)
  data <- na.omit(data)
  
  # Load libraries
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(patchwork); library(purrr); library(nlme); library(car); library(broom)
  
  # Define species and palette
  value_col      <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species    <- c("Spruce", "Pine")
  data$species   <- factor(data$species, levels = c(linear_species, exp_species))
  cb_palette     <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # Prepare data
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  models_linear <- list()
  for (sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    models_linear[[sp]] <- lm(avg_value ~ x, data = sp_data)
  }
  models_exp <- list()
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-9)
  for (sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    if (nrow(sp_data) < 5) { models_exp[[sp]] <- NULL; next }
    models_exp[[sp]] <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), data = sp_data,
          start = start_list, control = control_params)
    }, error = function(e) NULL)
  }
  
  ##############################################
  # Predictions for Panel A
  ##############################################
  pred_linear <- bind_rows(lapply(linear_species, function(sp) {
    df <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(df$x), max(df$x), length.out = 100)
    data.frame(species = sp, x = x_seq,
               pred = predict(models_linear[[sp]], newdata = data.frame(x = x_seq)))
  }))
  pred_exp <- bind_rows(lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]; if (is.null(mod)) return(NULL)
    df <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(df$x), max(df$x), length.out = 100)
    data.frame(species = sp, x = x_seq,
               pred = predict(mod, newdata = data.frame(x = x_seq)))
  }))
  pred_all <- bind_rows(pred_linear, pred_exp)
  
  # Panel A: combined
  p_combined <- ggplot() +
    geom_point(data = data_clean,
               aes(x = x, y = avg_value, color = species, shape = species, size = percentage), alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x), y = threshold,
             label = "median", fontface = "italic",
             hjust = -0.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette, name = "Species") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18), name = "Species") +
    # scale_size_continuous(name = "% pixels per bin", range = c(1, 8)) +
    scale_size_continuous(
      name = "Pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
    )
  
  
  #################################################
  # Stats for Panels B & C
  #################################################
  lin_stats <- bind_rows(lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]; summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]; slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]; p_val <- summ["x", "Pr(>|t|)"]
    x50 <- (threshold - intercept) / slope; abs_sl <- abs(slope)
    r2 <- summary(mod)$r.squared
    data.frame(species = sp, x50 = x50, abs_slope = abs_sl,
               se = se_slope, p_val = p_val, r_squared = r2)
  }))
  exp_stats <- bind_rows(lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]; if (is.null(mod)) return(NULL)
    co <- coef(mod); a <- co["a"]; b <- co["b"]; c <- co["c"]
    x50 <- ifelse((threshold - a) > 0 & b > 0, log((threshold - a)/b)/c, NA)
    slope50 <- c * (threshold - a); abs_sl <- abs(slope50)
    df_sp <- data_clean %>% filter(species == sp)
    r2 <- 1 - sum((df_sp[[value_col]] - predict(mod, df_sp))^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    dm <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"), parameterNames = c("a","b","c"))
    se <- dm$SE; df_res <- summary(mod)$df[2]
    t_val <- slope50 / se; p_val <- 2*(1-pt(abs(t_val), df_res))
    data.frame(species = sp, x50 = x50, abs_slope = abs_sl,
               se = se, p_val = p_val, r_squared = r2)
  }))
  stats_all <- bind_rows(lin_stats, exp_stats)
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: x50
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  # Panel C: absolute slope
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2), size = 5) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  # Combine panels
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position       = "bottom",
      legend.title          = element_blank(),
      legend.text           = element_text(size = 14),
      legend.key            = element_rect(fill = "white", color = NA),
      legend.background     = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(output_figure, final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_PSI_exp_linear_slope_coeff(df, 
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_coeff.png",
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_slope.png")

plot_Quantiles_TDiff_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data and remove missing rows
  data <- NDVI_TDiffbin(data) %>% na.omit()
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(nlme)    # for consistency
  library(car)     # for deltaMethod
  library(broom)
  
  # Define species groups and palette
  value_col       <- "avg_value"
  linear_species  <- c("Oak", "Beech")
  exp_species     <- c("Spruce", "Pine")
  data$species    <- factor(data$species, levels = c(linear_species, exp_species))
  cb_palette      <- c(
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
  # Fit Models by Species  #
  ##########################
  models_linear <- set_names(
    map(linear_species, ~ lm(avg_value ~ x, data = filter(data_clean, species == .x))),
    linear_species
  )
  models_exp <- set_names(
    map(exp_species, function(sp) {
      sp_data <- filter(data_clean, species == sp)
      tryCatch(
        nls(avg_value ~ a + b * exp(-c * x), data = sp_data,
            start = list(a = 5, b = 7, c = 0.04),
            control = nls.control(maxiter = 1200, minFactor = 1e-9)),
        error = function(e) NULL
      )
    }),
    exp_species
  )
  
  ##############################################
  # Linear Coefficients Extraction & Plotting  #
  ##############################################
  lin_coef_df <- map_dfr(linear_species, function(sp) {
    summ <- summary(models_linear[[sp]])$coefficients
    tibble(
      species     = sp,
      Coefficient = dplyr::recode(rownames(summ), `(Intercept)` = "a", x = "b"),
      Value        = summ[, "Estimate"],
      pvalue       = summ[, "Pr(>|t|)"]
    )
  }) %>%
    mutate(
      species = factor(species, levels = linear_species),
      label   = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue))
    )
  p_coeff_linear <- ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_col(position = position_dodge(0.9)) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(0.9), vjust = 0.5) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models", subtitle = expression(NDVI == a + b * x + epsilon), x = NULL, y = "coefficient value", caption = "(a)") +
    theme_minimal() + theme(legend.position = "top")
  
  #################################################
  # Exponential Coefficients Extraction & Plotting #
  #################################################
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
    geom_col(position = position_dodge(0.9)) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(0.9), vjust = 0.5) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models", subtitle = expression(NDVI == a + b * e^{-c * x} + epsilon), x = NULL, y = NULL, caption = "(b)") +
    theme_minimal() + theme(legend.position = "top")
  
  # Combine coefficient plots and save
  combined_coeff <- p_coeff_linear + p_coeff_exp
  print(combined_coeff)
  ggsave(combined_coef_fig, combined_coeff, width = 10, height = 8, dpi = 300)
  
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
    geom_point(data = data_clean,
               aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
               alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 32, y = threshold,
             label = "median", fontface = "italic",
             hjust = -0.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette, name = "Species") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       name = "Species") +
    scale_size_continuous(
      name = "Pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
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
    data.frame(species = sp, model = "exponential", x50 = x50_exp, abs_slope = abs_slope_exp,
               r_squared = r2_exp, se = se_exp, p_val = p_val_exp)
  })
  
  stats_all <- rbind(do.call(rbind, lin_stats), do.call(rbind, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: Bar plot of x50 (transpiration deficit) for all species
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "transpiration deficit (mm)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
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
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position       = "bottom",
      legend.title          = element_blank(),
      legend.text           = element_text(size = 14),
      legend.key            = element_rect(fill = "white", color = NA),
      legend.background     = element_blank(),
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
  
  print("Estimated transpiration deficit (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_exp_linear_slope_coeff(df,
                                            "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_linear_coeff.png",
                                            "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_linear_slope.png")

df_year <- final_df %>% filter(month %in% c("July", "August")) %>% filter(year < 2016)

plot_Quantiles_PSI_exp_linear_slope_coeff(df_year, 
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_coeff_2003_2015.png",
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_slope_2003_2015.png")

plot_Quantiles_TDiff_linear_all <- function(data, combined_coef_fig, output_figure) {
  # Process data and remove missing rows
  data <- NDVI_TDiffbin(data) %>% na.omit()
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(broom)
  
  # Define species groups and palette
  value_col       <- "avg_value"
  linear_species  <- c("Oak", "Beech")
  exp_species     <- c("Spruce", "Pine")
  all_species     <- c(linear_species, exp_species)
  data$species    <- factor(data$species, levels = all_species)
  cb_palette      <- c(
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
  # Fit Linear Models by Species
  ##########################
  models_linear <- set_names(
    map(all_species, ~ lm(avg_value ~ x, data = filter(data_clean, species == .x))),
    all_species
  )
  
  ##############################################
  # Coefficients Extraction & Plotting
  ##############################################
  coef_df <- map_dfr(all_species, function(sp) {
    summ <- summary(models_linear[[sp]])$coefficients
    tibble(
      species     = sp,
      Coefficient = dplyr::recode(rownames(summ), `(Intercept)` = "a", x = "b"),
      Value       = summ[, "Estimate"],
      pvalue      = summ[, "Pr(>|t|)"]
    )
  }) %>%
    mutate(
      species = factor(species, levels = all_species),
      label   = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue))
    )
  
  p_coeff <- ggplot(coef_df, aes(x = species, y = Value, fill = species)) +
    geom_col(position = position_dodge(0.9)) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(0.9), vjust = 0.5) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models", subtitle = expression(NDVI == a + b * x + epsilon),
         x = NULL, y = "coefficient value", caption = "(a)") +
    theme_minimal() +
    theme(
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 14),
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      strip.text        = element_text(size = 14, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  # Save coefficient figure
  ggsave(combined_coef_fig, p_coeff, width = 10, height = 8, dpi = 300)
  
  #################################################
  # Combined Final Plot with Model Predictions (Panel A)
  #################################################
  pred_all <- map_dfr(all_species, function(sp) {
    sp_data <- filter(data_clean, species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    pred <- predict(models_linear[[sp]], newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  
  p_combined <- ggplot() +
    geom_point(data = data_clean,
               aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
               alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = max(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", fontface = "italic", hjust = 1.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette, name = "Species") +
    scale_shape_manual(values = c(Oak=16, Beech=17, Spruce=15, Pine=18), name = "Species") +
    scale_size_continuous(
      name = "Pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 1), size = guide_legend(order = 2)) +
    labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
    )
  
  #################################################
  # Panel B & C: x50 and Absolute Slope
  #################################################
  stats_all <- map_dfr(all_species, function(sp) {
    mod <- models_linear[[sp]]
    summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]
    slope     <- summ["x", "Estimate"]
    se_slope  <- summ["x", "Std. Error"]
    p_val     <- summ["x", "Pr(>|t|)"]
    x50       <- (threshold - intercept) / slope
    abs_slope <- abs(slope)
    r2        <- summary(mod)$r.squared
    tibble(species = sp, x50 = x50, abs_slope = abs_slope, se = se_slope, p_val = p_val, r_squared = r2)
  }) %>%
    mutate(species = factor(species, levels = all_species))
  
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "transpiration deficit (mm)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title       = element_text(face = "bold", size = 16),
      axis.text        = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05, sprintf("%.2f*", r_squared), sprintf("%.2f", r_squared)), y = abs_slope/2),
              color = "black", size = 5) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title       = element_text(face = "bold", size = 16),
      axis.text        = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position       = "bottom",
      legend.title          = element_blank(),
      legend.text           = element_text(size = 14),
      legend.key            = element_rect(fill = "white", color = NA),
      legend.background     = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(output_figure, final_plot, width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_linear_all(df_year,
                               "results/key_displays_July_August/NDVI_Q_TDiffbin_linear_coeff_2003_2015.png",
                                "results/key_displays_July_August/NDVI_Q_TDiffbin_linear_slope_2003_2015.png")
