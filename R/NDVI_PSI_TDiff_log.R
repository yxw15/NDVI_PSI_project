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

# assign one color per species
cols <- rainbow(length(unique(NDVI_PSIbin_df$species)))
names(cols) <- unique(NDVI_PSIbin_df$species)

with(NDVI_PSIbin_df, 
     plot(-bin_median+8, avg_value,
          col   = cols[species],
          pch   = 19,
          xlab  = "Bin Median",
          ylab  = "Average Value",
          main  = "Average Value vs. Bin Median"))

legend("topleft", legend = names(cols), col = cols, pch = 19, title = "Species")

plot_NDVI_Q_PSIbin_log <- function(data, save_coeff_fig, save_slope_fig) {
  
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
  
  data <- NDVI_PSIbin_log(data)
  data <- na.omit(data)
  
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

plot_NDVI_Q_PSIbin_GAM <- function(data, save_slope_fig) {
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(tidyr)
  
  # Preprocess
  data_clean <- data %>%
    NDVI_PSIbin_log() %>%
    na.omit() %>%
    mutate(
      PSI = -bin_median + 8,
      percent_pixels = percentage * 100,
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine"))
    ) %>%
    filter(PSI > 0, !is.na(avg_value))
  
  # Settings
  threshold_psi <- 11.5
  cb_palette <- c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442")
  
  # Fit GAMs
  gam_models <- data_clean %>%
    group_by(species) %>%
    summarise(mod = list(gam(avg_value ~ s(PSI, bs = "cs", k = 5), data = cur_data())))
  
  # Compute metrics: PSI50, slope50, slope_se, r_squared, p_value
  stats_df <- gam_models %>% rowwise() %>% do({
    sp <- .$species; m <- .$mod
    # grid
    xs <- seq(min(data_clean$PSI[data_clean$species == sp]),
              max(data_clean$PSI[data_clean$species == sp]), length.out = 200)
    pr <- predict(m, newdata = tibble(PSI = xs), se.fit = TRUE)
    ys <- pr$fit; ses <- pr$se.fit
    # first crossing
    above <- which(ys >= threshold_psi)
    PSI50 <- if(length(above) > 0) xs[above[1]] else NA_real_
    # slope
    if(!is.na(PSI50)) {
      h <- xs[2] - xs[1]
      idx <- which.min(abs(xs - PSI50))
      if(idx > 1 && idx < length(xs)) {
        y1 <- ys[idx + 1]; y2 <- ys[idx - 1]
        s1 <- ses[idx + 1]; s2 <- ses[idx - 1]
        slope50 <- (y1 - y2) / (2 * h)
        slope_se <- sqrt(s1^2 + s2^2) / (2 * h)
      } else {
        slope50 <- NA_real_; slope_se <- NA_real_
      }
    } else {
      slope50 <- NA_real_; slope_se <- NA_real_
    }
    r_squared <- summary(m)$r.sq
    p_value <- if(!is.na(slope_se) && slope_se > 0) 2 * (1 - pnorm(abs(slope50 / slope_se))) else NA_real_
    tibble(species = sp, PSI50 = PSI50, slope50 = slope50,
           se = slope_se, p_val = p_value, r_squared = r_squared)
  })
  
  # Prediction lines
  pred_all <- gam_models %>% rowwise() %>% do({
    sp <- .$species; m <- .$mod
    xs <- seq(min(data_clean$PSI[data_clean$species == sp]),
              max(data_clean$PSI[data_clean$species == sp]), length.out = 200)
    tibble(species = sp, PSI = xs, pred = predict(m, newdata = tibble(PSI = xs)))
  })
  
  # Panel A: GAM fit
  p_combined <- ggplot(data_clean, aes(x = PSI, y = avg_value, color = species, size = percent_pixels)) +
    geom_point() +
    geom_line(data = pred_all, inherit.aes = FALSE,
              aes(x = PSI, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold_psi, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 0, y = threshold_psi, label = "median",
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
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
  
  # Panel B: PSI50
  p_x50 <- ggplot(stats_df, aes(x = species, y = PSI50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "transformed soil water potential") +
    ggtitle("(b)") + expand_limits(y = 0) +
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
  
  # Panel C: slope magnitude with error bars and r2/p-value
  p_slope <- ggplot(stats_df, aes(x = species, y = abs(slope50), fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs(slope50) - se), ymax = abs(slope50) + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  paste0(sprintf("R²=%.2f", r_squared), "*\n",
                                         sprintf("p=%.3f", p_val)),
                                  paste0("R²=", sprintf("%.2f", r_squared), "\n", sprintf("p=%.2f", p_val))),
                  y = abs(slope50) / 2), color = "black", size = 5) +
    scale_fill_manual(values = cb_palette, name = "") + guides(fill = "none") +
    labs(x = "", y = "absolute slope") + ggtitle("(c)") + expand_limits(y = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
  
  # Combine and save
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
  ggsave(save_slope_fig, final_slope_plot, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_GAM(df, "results/key_displays_July_August/Quantiles_PSI_bin_GAM_slope.png")

NDVI_TDiffbin_log <- function(df, bin_width_log = 0.5) {
  library(dplyr)
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Filter to positive deficits
  df <- df %>% filter(!is.na(transpiration_deficit), transpiration_deficit > 0)
  
  # Log-transform transpiration deficit
  df <- df %>% mutate(log_TDiff = log(transpiration_deficit))
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define log-scale bin breaks
  tdiff_min <- floor(min(df$log_TDiff, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$log_TDiff, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width_log)
  
  # Bin
  df <- df %>%
    mutate(TDiff_bin = cut(log_TDiff, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute stats
  result <- df %>%
    group_by(species, TDiff_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(TDiff_bin), function(bin_label) {
        edges <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(edges)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.01) %>%
    select(species, TDiff_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(result)
}

# assign one color per species
cols <- rainbow(length(unique(NDVI_TDiffbin_df$species)))
names(cols) <- unique(NDVI_TDiffbin_df$species)

with(NDVI_TDiffbin_df, 
     plot(-bin_median+8, avg_value,
          col   = cols[species],
          pch   = 19,
          xlab  = "Bin Median",
          ylab  = "Average Value",
          main  = "Average Value vs. Bin Median"))

legend("topleft", legend = names(cols), col = cols, pch = 19, title = "Species")

plot_NDVI_Q_TDiffbin_log <- function(data, save_coeff_fig, save_slope_fig) {
  
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
  
  data <- NDVI_TDiffbin_log(data)  # upstream calculation still based on PSI bins
  data <- na.omit(data)
  
  # compute transpiration deficit (TDiff) from PSI bin median
  data <- data %>% mutate(TDiff = -bin_median + 8) %>% filter(TDiff > 0)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(TDiff)) %>% mutate(percent_pixels = percentage * 100)
  threshold_tdiff <- 11.5  # keeps same cutoff for comparison
  
  # fit log models by species
  nls_models <- nlsList(avg_value ~ a * log(TDiff) + b | species,
                        data = data_clean,
                        start = list(a = -1, b = 12),
                        control = nls.control(maxiter = 1000))
  
  # extract coefficients and compute TDiff50 and slope at that point
  coef_df <- coef(nls_models) %>% as.data.frame() %>% rownames_to_column(var = "species") %>% filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  coef_df <- coef_df %>% mutate(
    TDiff50 = exp((threshold_tdiff - b) / a),
    slope50 = a / TDiff50
  )
  
  # generate predictions for plotting
  pred_list <- data_clean %>% group_by(species) %>% do({
    sp <- unique(.$species)
    x_seq <- seq(min(.$TDiff), max(.$TDiff), length.out = 100)
    pred <- predict(nls_models[[as.character(sp)]], newdata = data.frame(TDiff = x_seq))
    data.frame(TDiff = x_seq, pred = pred)
  })
  pred_all <- bind_rows(pred_list)
  
  # compute summary statistics
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- coefs["a"] / exp((threshold_tdiff - coefs["b"]) / coefs["a"])
    dm_result <- car::deltaMethod(mod,
                                  paste0("a/exp((", threshold_tdiff, "-b)/a)"),
                                  parameterNames = c("a", "b")
    )
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    p_val <- 2 * (1 - pt(abs(slope50_est / se), df_mod))
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) / sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    tibble(
      species = sp,
      slope50 = slope50_est,
      slope_abs = abs(slope50_est),
      se = se,
      p_val = p_val,
      r_squared = r_squared
    )
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  # Panel A: NDVI quantiles vs TDiff
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = TDiff, y = avg_value, color = species, size = percent_pixels)) +
    geom_line(data = pred_all, aes(x = TDiff, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold_tdiff, linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_size_continuous(name = "Pixel %", range = c(1, 6)) +
    annotate("text", x = 4, y = threshold_tdiff, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    labs(x = "transformed transpiration deficit", y = "NDVI quantiles (rank)") +
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
  
  # Panel B: TDiff50 by species
  p_x50 <- ggplot(coef_df, aes(x = species, y = TDiff50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "transformed transpiration deficit") +
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
  
  # Panel C: slope magnitude
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
  
  # Coefficient panel
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
  
  coeffs_long <- left_join(
    coeffs_long,
    coeff_stats %>% select(species, Coefficient, `Std. Error`, `Pr(>|t|)`),
    by = c("species", "Coefficient")
  ) %>%
    mutate(label = if_else(`Pr(>|t|)` < 0.05, "*", sprintf("%.2f", `Pr(>|t|)`)))
  
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = Value - `Std. Error`, ymax = Value + `Std. Error`), position = position_dodge(width = 0.9), width = 0.2) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(width = 0.9), color = "black", size = 5) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    scale_fill_manual(values = cb_palette, name = "") +
    labs(
      title = "NDVI quantiles ~ transformed transpiration deficit",
      subtitle = expression(NDVI == a * log(T[d]) + b + epsilon),
      x = "", y = "coefficient value"
    ) +
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
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  ggsave(filename = save_coeff_fig, plot = p_coeffs, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_TDiffbin_log(df, 
                       "results/key_displays_July_August/Quantiles_TDiff_bin_log_coeff.png",
                       "results/key_displays_July_August/Quantiles_TDiff_bin_log_slope.png")

plot_NDVI_Q_TDiffbin_GAM <- function(data, save_slope_fig) {
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(tidyr)
  
  # Preprocess
  data_clean <- data %>%
    NDVI_TDiffbin_log() %>%
    na.omit() %>%
    mutate(
      TDiff = -bin_median + 8,
      percent_pixels = percentage * 100,
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine"))
    ) %>%
    filter(TDiff > 0, !is.na(avg_value))
  
  # Settings
  threshold_tdiff <- 11.5
  cb_palette <- c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442")
  
  # Fit GAMs
  gam_models <- data_clean %>%
    group_by(species) %>%
    summarise(mod = list(gam(avg_value ~ s(TDiff, bs = "cs", k = 5), data = cur_data())))
  
  # Compute metrics: TDiff50, slope50, slope_se, r_squared, p_value
  stats_df <- gam_models %>% rowwise() %>% do({
    sp <- .$species; m <- .$mod
    # grid
    xs <- seq(min(data_clean$TDiff[data_clean$species == sp]),
              max(data_clean$TDiff[data_clean$species == sp]), length.out = 200)
    pr <- predict(m, newdata = tibble(TDiff = xs), se.fit = TRUE)
    ys <- pr$fit; ses <- pr$se.fit
    # first crossing
    above <- which(ys >= threshold_tdiff)
    TDiff50 <- if(length(above)>0) xs[above[1]] else NA_real_
    # slope
    if(!is.na(TDiff50)) {
      h <- xs[2] - xs[1]
      idx <- which.min(abs(xs - TDiff50))
      if(idx>1 && idx<length(xs)) {
        y1 <- ys[idx+1]; y2 <- ys[idx-1]
        s1 <- ses[idx+1]; s2 <- ses[idx-1]
        slope50 <- (y1 - y2)/(2*h)
        slope_se <- sqrt(s1^2 + s2^2)/(2*h)
      } else {
        slope50 <- NA_real_; slope_se <- NA_real_
      }
    } else {
      slope50 <- NA_real_; slope_se <- NA_real_
    }
    r_squared <- summary(m)$r.sq
    p_value <- if(!is.na(slope_se) && slope_se>0) 2*(1 - pnorm(abs(slope50/slope_se))) else NA_real_
    tibble(species = sp, TDiff50 = TDiff50, slope50 = slope50,
           se = slope_se, p_val = p_value, r_squared = r_squared)
  })
  
  # Prediction lines
  pred_all <- gam_models %>% rowwise() %>% do({
    sp <- .$species; m <- .$mod
    xs <- seq(min(data_clean$TDiff[data_clean$species == sp]),
              max(data_clean$TDiff[data_clean$species == sp]), length.out = 200)
    tibble(species = sp, TDiff = xs, pred = predict(m, newdata = tibble(TDiff = xs)))
  })
  
  # Panel A: GAM fit
  p_combined <- ggplot(data_clean, aes(x = TDiff, y = avg_value, color = species, size = percent_pixels)) +
    geom_point() +
    geom_line(data = pred_all, inherit.aes = FALSE,
              aes(x = TDiff, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold_tdiff, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 4, y = threshold_tdiff, label = "median",
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_size_continuous(name = "Pixel %", range = c(1, 6)) +
    labs(x = "transformed transpiration deficit", y = "NDVI quantiles (rank)") +
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
  
  # Panel B: TDiff50
  p_x50 <- ggplot(stats_df, aes(x = species, y = TDiff50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "transformed transpiration deficit") +
    ggtitle("(b)") + expand_limits(y = 0) +
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
  
  # Panel C: slope magnitude with error bars and r2/p-value
  p_slope <- ggplot(stats_df, aes(x = species, y = abs(slope50), fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs(slope50) - se), ymax = abs(slope50) + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  paste0(sprintf("R²=%.2f", r_squared), "*\n",
                                         sprintf("p=%.3f", p_val)),
                                  paste0("R²=", sprintf("%.2f", r_squared), "\n", sprintf("p=%.2f", p_val))),
                  y = abs(slope50)/2), color = "black", size = 5) +
    scale_fill_manual(values = cb_palette, name = "") + guides(fill = "none") +
    labs(x = "", y = "absolute slope") + ggtitle("(c)") + expand_limits(y = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
  
  # Combine and save
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
  ggsave(save_slope_fig, final_slope_plot, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_TDiffbin_GAM(df, "results/key_displays_July_August/Quantiles_TDiff_bin_GAM_slope.png")

plot_log_vs_linear_AIC <- function(data, save_combined_fig) {
  
  library(dplyr);   library(purrr)
  library(broom);   library(tidyr)
  library(ggplot2); library(patchwork)
  
  species_order   <- c("Oak","Beech","Spruce","Pine")
  model_palette   <- c(linear="orange", logarithm="dodgerblue")
  base_size       <- 14
  
  #### Panel A: NDVI ~ PSI ####
  data_a <- NDVI_PSIbin_log(data) %>%
    na.omit() %>%
    mutate(
      species = factor(species, levels=species_order),
      x_lin   = -bin_median + 8,
      x_log   = -bin_median + 8
    ) %>%
    filter(x_log > 0)
  
  aic_df_a <- data_a %>%
    group_by(species) %>%
    nest() %>%
    mutate(
      fit_lin = map(data, ~ lm(avg_value ~ x_lin,      data = .x)),
      fit_log = map(data, ~ lm(avg_value ~ log(x_log), data = .x)),
      aic_lin = map_dbl(fit_lin, AIC),
      aic_log = map_dbl(fit_log, AIC),
      r2_lin  = map_dbl(fit_lin, ~ summary(.x)$r.squared),
      r2_log  = map_dbl(fit_log, ~ summary(.x)$r.squared)
    ) %>%
    select(species, aic_lin, aic_log, r2_lin, r2_log) %>%
    pivot_longer(
      cols      = starts_with(c("aic_","r2_")),
      names_to  = c(".value","model"),
      names_sep = "_"
    ) %>%
    mutate(
      model = factor(model, levels = c("lin","log"),
                     labels = c("linear","logarithm")),
      y_pos  = aic / 2
    )
  
  dodge_width <- 0.8
  # bar_width   <- dodge_width / 2 - 0.05
  bar_width   <- dodge_width
  
  p_a <- ggplot(aic_df_a, aes(x = species, y = aic, fill = model)) +
    geom_col(position = position_dodge(width = dodge_width), width = bar_width) +
    geom_text(aes(label = round(r2,2), y = y_pos),
              position = position_dodge(width = dodge_width), vjust = 0.5, size = 5) +
    scale_fill_manual(values = model_palette, name = "") +
    labs(
      x   = NULL,
      y   = "AIC",
      title = "",
      tag = "(a)"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.tag = element_text(face = "bold", size = 16),
      plot.tag.position = c(0,1),
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
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  #### Panel B: NDVI ~ TDiff ####
  data_b <- NDVI_TDiffbin_log(data) %>%
    na.omit() %>%
    mutate(
      species = factor(species, levels = species_order),
      x       = -bin_median + 8
    ) %>%
    filter(x > 0)
  
  aic_df_b <- data_b %>%
    group_by(species) %>%
    nest() %>%
    mutate(
      fit_lin = map(data, ~ lm(avg_value ~ x,      data = .x)),
      fit_log = map(data, ~ lm(avg_value ~ log(x), data = .x)),
      aic_lin = map_dbl(fit_lin, AIC),
      aic_log = map_dbl(fit_log, AIC),
      r2_lin  = map_dbl(fit_lin, ~ summary(.x)$r.squared),
      r2_log  = map_dbl(fit_log, ~ summary(.x)$r.squared)
    ) %>%
    select(species, aic_lin, aic_log, r2_lin, r2_log) %>%
    pivot_longer(
      cols      = starts_with(c("aic_","r2_")),
      names_to  = c(".value","model"),
      names_sep = "_"
    ) %>%
    mutate(
      model = factor(model, levels = c("lin","log"),
                     labels = c("linear","logarithm")),
      y_pos  = aic / 2
    )
  
  p_b <- ggplot(aic_df_b, aes(x = species, y = aic, fill = model)) +
    geom_col(position = position_dodge(width = dodge_width), width = bar_width) +
    geom_text(aes(label = round(r2,2), y = y_pos),
              position = position_dodge(width = dodge_width), vjust = 0.5, size = 5) +
    scale_fill_manual(values = model_palette, name = "") +
    labs(
      x   = NULL,
      y   = "AIC",
      title = "",
      tag = "(b)"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.tag = element_text(face = "bold", size = 16),
      plot.tag.position = c(0,1),
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
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Combine, save, print
  combined <- (p_a | p_b) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  dir.create(dirname(save_combined_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_combined_fig, combined, width = 14, height = 6, dpi = 300)
  print(combined)
}

plot_log_vs_linear_AIC(df, "results/key_displays_July_August/log_linear_AIC.png")
