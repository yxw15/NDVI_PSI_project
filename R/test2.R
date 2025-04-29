library(dplyr)

load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

df <- final_df %>% filter(month %in% c("July", "August"))

NDVI_PSIbin_log <- function(df, bin_width_log = 0.5) {
  library(dplyr)
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Filter to valid, non-zero potentials
  df <- df %>% filter(!is.na(soil_water_potential), soil_water_potential != 0)
  
  # Log-transform magnitude of soil water potential
  df <- df %>% mutate(log_PSI = log(abs(soil_water_potential)))
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define log-scale bin breaks
  psi_min <- floor(min(df$log_PSI, na.rm = TRUE))
  psi_max <- ceiling(max(df$log_PSI, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width_log)
  
  # Bin log(soil_water_potential magnitude)
  df <- df %>%
    mutate(PSI_bin = cut(log_PSI, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute stats per species and bin
  result <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        edges <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(edges)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.01) %>%
    select(species, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(result)
}

plot_NDVI_Q_PSIbin_log <- function(data, save_coeff_fig, save_slope_fig) {
  
  data <- NDVI_PSIbin_log(data)
  data <- na.omit(data)
  
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)
  library(tidyr)
  
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce" = "#009E73", "Pine" = "#F0E442")
  
  data <- data %>% mutate(x = -bin_median + 8) %>% filter(x > 0)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x)) %>% mutate(percent_pixels = percentage * 100)
  threshold <- 11.5
  
  nls_models <- nlsList(avg_value ~ a * log(x) + b | species,
                        data = data_clean,
                        start = list(a = -1, b = 12),
                        control = nls.control(maxiter = 1000))
  
  coef_df <- coef(nls_models) %>% as.data.frame() %>% 
    rownames_to_column(var = "species") %>% 
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  coef_df <- coef_df %>% mutate(x50 = exp((threshold - b) / a), slope50 = a / x50)
  
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x), max(.$x), length.out = 100)
      pred <- predict(nls_models[[as.character(sp)]], newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- coefs["a"] / exp((threshold - coefs["b"]) / coefs["a"])
    dm_result <- deltaMethod(mod, paste0("a / exp((", threshold, " - b)/a)"), parameterNames = c("a", "b"))
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    p_val <- 2 * (1 - pt(abs(slope50_est / se), df_mod))
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) / sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    tibble(species = sp, slope50 = slope50_est, slope_abs = abs(slope50_est), se = se, p_val = p_val, r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  # Panel A
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species, size = percent_pixels)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_size_continuous(name = "Pixel %", range = c(1, 6)) +
    labs(x = "transformed soil water potential", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14)
    )
  
  # Panel B
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "transformed soil water potential") +
    ggtitle("(b)") +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Panel C
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05, sprintf("%.2f*", r_squared), sprintf("%.2f", r_squared)), y = slope_abs/2), color = "black", size = 5) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  final_slope_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  # Coefficient plot
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  
  coeffs_long <- coef_df %>%
    select(species, a, b) %>%
    pivot_longer(cols = c("a", "b"), names_to = "Coefficient", values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, `Std. Error`, `Pr(>|t|)`),
                           by = c("species", "Coefficient")) %>%
    mutate(label = if_else(`Pr(>|t|)` < 0.05, "*", sprintf("%.2f", `Pr(>|t|)`)))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = Value - `Std. Error`, ymax = Value + `Std. Error`), position = position_dodge(width = 0.9), width = 0.2) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(width = 0.9), color = "black", size = 5) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    scale_fill_manual(values = cb_palette, name = "") +
    labs(title = "NDVI quantiles ~ soil water potential",
         subtitle = expression(NDVI == a * log(x) + b + epsilon),
         x = "", y = "coefficient value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  ggsave(filename = save_coeff_fig, plot = p_coeffs, width = 10, height = 8, dpi = 300)
}
plot_NDVI_Q_PSIbin_log(df, 
                       "results/key_displays_July_August/Quantiles_PSI_bin_log_coeff.png",
                       "results/key_displays_July_August/Quantiles_PSI_bin_log_slope.png")
