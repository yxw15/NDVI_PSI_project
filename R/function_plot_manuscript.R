setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(dplyr)

load("results/Data/AllSpecies_AllMonths_depth50.RData")
combined_df_50 <- combined %>% filter(Quantiles > 0)
combined_df_50 <- combined_df_50 %>% filter(month %in% c("July"))

load("results/Data/AllSpecies_AllMonths_depth100.RData")
combined_df_100 <- combined %>% filter(Quantiles > 0)
combined_df_100 <- combined_df_100 %>% filter(month %in% c("July"))

load("results/Data/AllSpecies_AllMonths_depth150.RData")
combined_df_150 <- combined %>% filter(Quantiles > 0)
combined_df_150 <- combined_df_150 %>% filter(month %in% c("July"))

### Figure 1 ###
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
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Create an x-axis scale that shows breaks every 2 years from 2003 to 2024
  x_scale <- scale_x_continuous(breaks = seq(2003, 2024, by = 2), limits = c(2003, 2024))
  
  ### TIME SERIES PLOTS ###
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), color = "species") +
    scale_color_manual(values = cb_palette, name = "", guide = guide_legend(nrow = 1)) +
    custom_theme +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    x_scale +
    annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold")
  
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("soil water potential"), bold("(kPa)"))), color = "species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    x_scale +
    annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold")
  
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("transpiration deficit"), bold("(mm)"))), color = "species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    x_scale +
    annotate("text", x = -Inf, y = Inf, label = "(c)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold")
  
  ts_plot <- p1 / p2 / p3
  
  ### CORRELATION PLOTS ###
  annotations_d <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_psi, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  print(annotations_d)
  
  annotations_e <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_tdiff, .data[[ndvi_column]], use = "complete.obs"), 2))
  print(annotations_e)
  
  # text_d <- paste0(species_order, ": r = ", 
  #                  annotations_d$corr[match(species_order, annotations_d$species)],
  #                  collapse = "\n")
  # 
  # text_e <- paste0(species_order, ": r = ", 
  #                  annotations_e$corr[match(species_order, annotations_e$species)],
  #                  collapse = "\n")
  
  p_d <- ggplot(df_species, aes(x = avg_psi, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "soil water potential (kPa)", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), color = "species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = -Inf, y = Inf, label = "(d)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold") 
  # annotate("text", x = -Inf, y = Inf, label = text_d, hjust = 0, vjust = 4.5, size = 5, color = "black")
  
  p_e <- ggplot(df_species, aes(x = avg_tdiff, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "transpiration deficit (mm)", y = "", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = Inf, y = Inf, label = "(e)", hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold") 
  # annotate("text", x = Inf, y = Inf, label = text_e, hjust = 1, vjust = 1.5, size = 5, color = "black")
  
  corr_plot <- p_d | p_e
  
  combined_main <- ts_plot / corr_plot
  combined_main <- combined_main & theme(legend.position = "top",
                                         plot.background = element_rect(fill = "white", color = NA))
                                         # plot.margin = margin(0,0,0,0, "cm"))
  
  legend <- suppressWarnings(cowplot::get_legend(p1))
  
  final_plot <- cowplot::plot_grid(combined_main, legend, ncol = 1, rel_heights = c(1, 0.05)) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  print(final_plot)
  ggsave(filename = output_path, plot = final_plot, width = 12, height = 14, dpi = 300)
}

plot_time_series_and_correlation_combined(combined_df_50,
                                          "results/Figures/time_series_Quantiles_PSI_TDiff_species_depth50.png")
### Figure 2 ###
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

plot_Quantiles_PSI_exp_linear_slope_coeff <- function(data_input, combined_coef_fig, output_figure) {
  # Process data
  data <- NDVI_PSIbin(data_input)
  data <- na.omit(data)
  
  # Load libraries
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(patchwork); library(purrr); library(nlme); library(car); library(broom)
  
  # Define species and palette
  value_col      <- "avg_value"
  linear_species <- c("Oak", "Beech", "Spruce")
  exp_species    <- c("Pine")
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

plot_Quantiles_PSI_exp_linear_slope_coeff(combined_df_50, 
                                          "results/Figures/NDVI_Q_PSIbin_exp_linear_coeff_depth50.png",
                                          "results/Figures/NDVI_Q_PSIbin_exp_linear_slope_depth50.png")

plot_Quantiles_PSI_exp_linear_slope_coeff(combined_df_100, 
                                          "results/Figures/NDVI_Q_PSIbin_exp_linear_coeff_depth100.png",
                                          "results/Figures/NDVI_Q_PSIbin_exp_linear_slope_depth100.png")

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

plot_Quantiles_PSI_exp_linear_slope_coeff(combined_df_150, 
                                          "results/Figures/NDVI_Q_PSIbin_exp_linear_coeff_depth150.png",
                                          "results/Figures/NDVI_Q_PSIbin_exp_linear_slope_depth150.png")

### Figure 3 ###
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

plot_Quantiles_TDiff_exp_linear_slope_coeff(combined_df_50,
                                            "results/Figures/NDVI_Q_TDiffbin_exp_linear_coeff.png",
                                            "results/Figures/NDVI_Q_TDiffbin_exp_linear_slope.png")


