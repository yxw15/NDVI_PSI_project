setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(dplyr)
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")
df <- final_df %>% filter(month %in% c("July", "August"))

plot_distribution <- function(input_data, output_figure) {
  
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  # custom theme
  custom_theme <- theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14)
  )
  
  # (a) Quantiles: 1-unit bins from 1 to 22, % on y-axis
  p1 <- ggplot(input_data, aes(x = Quantiles)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count) * 100)),
      binwidth = 1, boundary = 0.5,
      fill   = "#009E73", colour = "black"
    ) +
    scale_x_continuous(
      breaks = 1:22,
      limits = c(0.5, 22.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(scale = 1)
    ) +
    labs(
      x   = "NDVI quantiles",
      y   = "",
      tag = "(a)"
    ) +
    custom_theme +
    theme(
      plot.tag.position = c(0.02, 0.98),
      plot.tag          = element_text(face = "bold", size = 16)
    )
  
  # (b) Soil water potential: 50‐unit bins, % on y-axis
  p2 <- ggplot(input_data, aes(x = soil_water_potential)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count) * 100)),
      binwidth = 50,
      fill     = "#E69F00",
      colour   = "black"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(scale = 1)
    ) +
    labs(
      x   = "soil water potential (kPa)",
      y   = "",
      tag = "(b)"
    ) +
    custom_theme +
    theme(
      plot.tag.position = c(0.02, 0.98),
      plot.tag          = element_text(face = "bold", size = 16)
    )
  
  # (c) Transpiration deficit: 3‐unit bins, % on y-axis
  p3 <- ggplot(input_data, aes(x = transpiration_deficit)) +
    geom_histogram(
      aes(y = after_stat(count / sum(count) * 100)),
      binwidth = 3,
      fill     = "#0072B2",
      colour   = "black"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(scale = 1)
    ) +
    labs(
      x   = "transpiration deficit (mm)",
      y   = "",
      tag = "(c)"
    ) +
    custom_theme +
    theme(
      plot.tag.position = c(0.02, 0.98),
      plot.tag          = element_text(face = "bold", size = 16)
    )
  
  # stack vertically
  combined <- (p1 / p2 / p3) + plot_layout(ncol = 1)
  
  # save to file
  ggsave(
    filename = output_figure,
    plot     = combined,
    width    = 8,      # in inches
    height   = 12,     # in inches
    dpi      = 300
  )
  
  invisible(combined)
}

plot_distribution(df, "results/key_displays_July_August/distribution.png")

plot_distribution_log <- function(input_data, output_figure) {
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  # copy & augment
  data <- input_data
  data$pF <- log10(abs(data$soil_water_potential))
  data$pT <- log10(abs(data$transpiration_deficit))
  data <- subset(data, is.finite(pF) & is.finite(pT))
  
  # full 0.3‐step sequences
  pf_breaks <- seq(floor(min(data$pF)/0.3)*0.3,
                   ceiling(max(data$pF)/0.3)*0.3,
                   by = 0.3)
  pt_breaks <- seq(floor(min(data$pT)/0.3)*0.3,
                   ceiling(max(data$pT)/0.3)*0.3,
                   by = 0.3)
  # every‐3rd, then drop the first element
  pf_ticks <- pf_breaks[seq(1, length(pf_breaks), by = 3)][-1]
  pt_ticks <- pt_breaks[seq(1, length(pt_breaks), by = 3)][-1]
  
  custom_theme <- theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", colour = "white"),
    panel.background  = element_rect(fill = "white"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(colour = "black", size = 14),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold")
  )
  
  # (a) Quantiles every 3rd
  p1 <- ggplot(data, aes(Quantiles)) +
    geom_histogram(aes(y = after_stat(count/sum(count)*100)),
                   binwidth = 1, boundary = 0.5,
                   fill = "#009E73", colour = "black") +
    scale_x_continuous(breaks = seq(1, 22, by = 3),
                       limits = c(0.5, 22.5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = "NDVI quantiles", y = "", tag = "(a)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (b) pF dropping smallest tick
  p2 <- ggplot(data, aes(pF)) +
    geom_histogram(aes(y = after_stat(count/sum(count)*100)),
                   binwidth = 0.3, boundary = pf_breaks[1],
                   fill = "#E69F00", colour = "black") +
    scale_x_continuous(breaks = pf_ticks,
                       limits = range(pf_breaks), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = expression(atop("(log"[10]~"|"~Psi[~soil]~"|)")),
         y = "", tag = "(b)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (c) pT dropping smallest tick
  p3 <- ggplot(data, aes(pT)) +
    geom_histogram(aes(y = after_stat(count/sum(count)*100)),
                   binwidth = 0.3, boundary = pt_breaks[1],
                   fill = "#0072B2", colour = "black") +
    scale_x_continuous(breaks = pt_ticks,
                       limits = range(pt_breaks), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = expression(log[10](abs(T[d]))),
         y = "", tag = "(c)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  combined <- (p1 / p2 / p3) + plot_layout(ncol = 1)
  ggsave(output_figure, combined, width = 8, height = 12, dpi = 300)
  invisible(combined)
}

plot_distribution_log(df, "results/key_displays_July_August/distribution_log.png")

plot_distribution_combined <- function(input_data, output_figure) {
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  # custom theme
  custom_theme <- theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14)
  )
  
  # (a) NDVI quantiles
  p1 <- ggplot(input_data, aes(x = Quantiles)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 1, boundary = 0.5,
                   fill   = "#009E73", colour = "black") +
    scale_x_continuous(breaks = 1:22,
                       limits = c(0.5, 22.5),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = "NDVI quantiles", y = NULL, tag = "(a)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (b) Soil water potential original
  p2 <- ggplot(input_data, aes(x = soil_water_potential)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 50, fill = "#E69F00", colour = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = "soil water potential (kPa)", y = NULL, tag = "(b)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (c) Soil water potential log-transformed
  data_log <- input_data
  data_log$pF <- log10(abs(data_log$soil_water_potential))
  data_log <- subset(data_log, is.finite(pF))
  pf_breaks <- seq(floor(min(data_log$pF) / 0.3) * 0.3,
                   ceiling(max(data_log$pF) / 0.3) * 0.3,
                   by = 0.3)
  pf_ticks  <- pf_breaks[seq(1, length(pf_breaks), by = 3)][-1]
  p3 <- ggplot(data_log, aes(x = pF)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 0.3, boundary = pf_breaks[1],
                   fill = "#E69F00", colour = "black") +
    scale_x_continuous(breaks = pf_ticks,
                       limits = range(pf_breaks), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(
      x = expression(paste(log[10], " | ", Psi[soil], " |")),
      y = NULL,
      tag = "(c)"
    ) +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (d) Transpiration deficit original
  p4 <- ggplot(input_data, aes(x = transpiration_deficit)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 3, fill = "#0072B2", colour = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = "transpiration deficit (mm)", y = NULL, tag = "(d)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (e) Transpiration deficit log-transformed
  data_log$pT <- log10(abs(data_log$transpiration_deficit))
  data_log <- subset(data_log, is.finite(pT))
  pt_breaks <- seq(floor(min(data_log$pT) / 0.3) * 0.3,
                   ceiling(max(data_log$pT) / 0.3) * 0.3,
                   by = 0.3)
  pt_ticks  <- pt_breaks[seq(1, length(pt_breaks), by = 3)][-1]
  p5 <- ggplot(data_log, aes(x = pT)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 0.3, boundary = pt_breaks[1],
                   fill = "#0072B2", colour = "black") +
    scale_x_continuous(breaks = pt_ticks,
                       limits = range(pt_breaks), expand = c(0, 0)) +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
    #                    labels = percent_format(scale = 1)) +
    scale_y_continuous(
      breaks = c(0, 5, 10, 15, 20),  
      labels = scales::percent_format(scale = 1),
      expand = expansion(mult = c(0, 0.1))
    ) +
    labs(x = expression(log[10](abs(T[d]))), y = NULL, tag = "(e)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # combine layout: a on top, b/c second row, d/e third row
  combined <- p1 / (p2 | p3) / (p4 | p5) + plot_layout(ncol = 1)
  
  # save to file
  ggsave(
    filename = output_figure,
    plot     = combined,
    width    = 10,
    height   = 14,
    dpi      = 300
  )
  invisible(combined)
}

plot_distribution_combined(df, "results/key_displays_July_August/distribution_combined.png")

plot_distribution_combined_PSI <- function(input_data, output_figure) {
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  # custom theme
  custom_theme <- theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14)
  )
  
  # (b) Soil water potential original
  p2 <- ggplot(input_data, aes(x = soil_water_potential)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 50, fill = "#E69F00", colour = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(x = "soil water potential (kPa)", y = NULL, tag = "(a)") +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # (c) Soil water potential log-transformed
  data_log <- input_data
  data_log$pF <- log10(abs(data_log$soil_water_potential))
  data_log <- subset(data_log, is.finite(pF))
  pf_breaks <- seq(floor(min(data_log$pF) / 0.3) * 0.3,
                   ceiling(max(data_log$pF) / 0.3) * 0.3,
                   by = 0.3)
  pf_ticks  <- pf_breaks[seq(1, length(pf_breaks), by = 3)][-1]
  p3 <- ggplot(data_log, aes(x = pF)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   binwidth = 0.3, boundary = pf_breaks[1],
                   fill = "#E69F00", colour = "black") +
    scale_x_continuous(breaks = pf_ticks,
                       limits = range(pf_breaks), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = percent_format(scale = 1)) +
    labs(
      x = expression(paste(log[10], " | ", Psi[soil], " |")),
      y = NULL,
      tag = "(b)"
    ) +
    custom_theme +
    theme(plot.tag.position = c(0.02, 0.98),
          plot.tag          = element_text(face = "bold", size = 16))
  
  # combine layout: a on top, b/c second row, d/e third row
  combined <- p2 | p3
  
  # save to file
  ggsave(
    filename = output_figure,
    plot     = combined,
    width    = 10,
    height   = 5,
    dpi      = 300
  )
  invisible(combined)
}

plot_distribution_combined_PSI(df, "results/key_displays_July_August/distribution_combined_psi.png")

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
      Model = c("linear", "exponential"),
      AIC = c(aic_linear, aic_exp),
      R2 = c(r2_linear, r2_exp)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_a[[sp]] <- sp_results
  }
  aic_df_a <- do.call(rbind, aic_results_a)
  aic_df_a$species <- factor(aic_df_a$species, levels = species_order)
  
  # Shared color palette for Panels A and B
  model_palette_shared <- c("linear" = "orange",
                            "exponential" = "dodgerblue")
  
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
      Model = c("linear", "exponential"),
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
      Model = c("linear", "poly2", "poly3"),
      AIC = c(aic_linear, aic_poly2, aic_poly3),
      R2 = c(r2_linear, r2_poly2, r2_poly3)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_c[[sp]] <- sp_results
  }
  aic_df_c <- do.call(rbind, aic_results_c)
  aic_df_c$species <- factor(aic_df_c$species, levels = species_order)
  
  model_palette_c <- c("linear" = "orange",
                       "poly2"  = "dodgerblue",
                       "poly3"  = "green4")
  
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
    annotate("text", x = 0, y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
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
  
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
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
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(filename = save_coeff_fig, plot = p_coeffs, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_log(df, 
                       "results/key_displays_July_August/Quantiles_PSI_bin_log_coeff.png",
                       "results/key_displays_July_August/Quantiles_PSI_bin_log_slope.png")
