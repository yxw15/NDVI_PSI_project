library(tidyverse)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

NDVI_TDiffbin_df <- NDVI_TDiffbin(all_results_df)

plot_Quantiles_TDiff_poly_3_slope_coeff <- function(data, coef_fig, output_figure) {
  require(ggplot2)
  require(nlme)
  require(dplyr)
  require(tibble)
  require(patchwork)
  require(purrr)
  require(car)      # for deltaMethod
  require(broom)
  require(tidyr)
  
  value_col <- "avg_value"
  # Order species as "Oak", "Beech", "Spruce", "Pine"
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce" = "#009E73", "Pine" = "#F0E442")
  
  # Prepare data using your custom NDVI_TDiffbin function and set x = bin_median
  data <- NDVI_TDiffbin(data)
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  data_clean$species <- factor(data_clean$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  threshold <- 11.5
  
  # Fit cubic model using poly(x, 3, raw=TRUE) so that the model is: a + b*x + c*x^2 + d*x^3
  start_list_poly <- list(a = 1, b = 0.1, c = 0.1, d = 0.01)
  nls_models_poly <- nlsList(
    avg_value ~ a + b * poly(x, 3, raw = TRUE)[,1] + c * poly(x, 3, raw = TRUE)[,2] + d * poly(x, 3, raw = TRUE)[,3] | species,
    data = data_clean,
    start = start_list_poly,
    control = nls.control()
  )
  print(summary(nls_models_poly))
  
  tidy_poly <- map_df(nls_models_poly, tidy, .id = "species")
  tidy_poly$species <- factor(tidy_poly$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  tidy_poly <- tidy_poly %>% mutate(label = if_else(p.value < 0.05, "*", sprintf("%.2f", p.value)))
  
  p_coeff <- ggplot(tidy_poly, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "Cubic Coefficients by Species for NDVI",
         subtitle = expression(NDVI == a + b*x + c*x^2 + d*x^3),
         x = "Species", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold", size = 12))
  
  print(p_coeff)
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # PANEL A: Observed Data and Fitted Curves
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models_poly[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list) %>% mutate(pred = as.numeric(pred))
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold, label = "median",
             hjust = -2.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI", caption = "(a)") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  calc_x50_exact <- function(a, b, c, d, threshold, x_range) {
    # Set up polynomial coefficients for: d*x^3 + c*x^2 + b*x + (a - threshold) = 0
    coeffs <- c(d, c, b, a - threshold)
    roots <- polyroot(coeffs)
    # Keep only roots with negligible imaginary part
    real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
    
    if (length(real_roots) == 0) {
      return(NA)
    }
    
    # Find roots that lie within the observed x range
    valid_roots <- real_roots[real_roots >= x_range[1] & real_roots <= x_range[2]]
    
    if (length(valid_roots) > 0) {
      # If multiple valid roots exist, choose the one closest to the median
      x50 <- valid_roots[which.min(abs(valid_roots - median(x_range)))]
    } else {
      # Otherwise, choose the real root closest to the median of x_range
      x50 <- real_roots[which.min(abs(real_roots - median(x_range)))]
    }
    
    return(x50)
  }
  
  
  stats_df <- do.call(rbind, stats_list)
  
  print(stats_df %>% select(species, x50, slope_abs, se, slope_p))
  
  rsq_df <- map_df(levels(data_clean$species), function(sp) {
    sp_data <- filter(data_clean, species == sp)
    mod <- nls_models_poly[[sp]]
    pred <- predict(mod, newdata = sp_data)
    r2_val <- 1 - sum((sp_data[[value_col]] - pred)^2) / sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    data.frame(species = sp, r2 = r2_val)
  })
  
  stats_df <- left_join(stats_df, rsq_df, by = "species")
  stats_df <- stats_df %>%
    mutate(slope_label = if_else(slope_p < 0.05,
                                 sprintf("%.2f*", r2),
                                 sprintf("p=%.2f\nR²=%.2f", slope_p, r2)))
  
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # PANEL B: Bar Plot of x50 Values
  p_x50 <- ggplot(stats_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Species", y = "Transpiration Deficit", caption = "(b)") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Labels
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs / 2, label = slope_label), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope", caption = "(c)") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    labs(caption = "(c)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  ggsave(filename = output_figure,
         plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}
