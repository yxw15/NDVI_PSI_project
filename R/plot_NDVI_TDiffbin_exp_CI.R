setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)

data <- combined %>% filter(month %in% c("July", "August")) %>% filter(year < 2023)
data <- data %>% filter(Quantiles > 0)
data <- na.omit(data)

# ---- theme & axes ----
# 1. Visual Style only
base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# 2. Legend Ordering only
base_guides <- guides(
  size  = guide_legend(order = 1, override.aes = list(alpha = 1)),
  color = guide_legend(order = 2, override.aes = list(linewidth = 1.2, alpha = 1)),
  shape = guide_legend(order = 2)
)

data <- combined %>% filter(month %in% c("July", "August")) %>% filter(year < 2023)
data <- data %>% filter(Quantiles > 0)
data <- na.omit(data)

NDVI_TDiffbin <- function(df, bin_width = 3, min_count = 1000) {
  
  suppressPackageStartupMessages({
    library(tidyverse)
  })
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  # Bin the transpiration deficit values
  df <- df %>%
    mutate(
      TDiff_bin = cut(
        transpiration_deficit,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # Compute statistics for each species and bin
  meanNDVI_TDiffbin_species <- df %>%
    group_by(species, TDiff_bin) %>%
    summarise(
      avg_value = mean(Quantiles, na.rm = TRUE),
      bin_mean  = mean(transpiration_deficit, na.rm = TRUE),
      count     = n(),
      .groups   = "drop"
    ) %>%
    left_join(species_totals, by = "species") %>%
    filter(count >= min_count) %>%
    mutate(percentage = count / total_pixels)
  
  return(meanNDVI_TDiffbin_species)
}

NDVI_TDiffbin_df <- NDVI_TDiffbin(data)

# --- helper for ribbon+line legend glyph (same as F2) ---
ribbon_line_key <- function(data, params, size) {
  grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(
        col  = NA,
        fill = scales::alpha(data$colour %||% data$fill, 0.25)
      )
    ),
    grid::linesGrob(
      x = c(0.08, 0.92), y = c(0.5, 0.5),
      gp = grid::gpar(
        col   = data$colour %||% data$color,
        lwd   = (data$linewidth %||% data$size %||% 0.5) * .pt,
        lty   = data$linetype %||% 1,
        alpha = 1
      )
    )
  )
}

# ---- NDVI_TDiff_exp ----
plot_NDVI_TDiff_exp_slope_coeff <- function(data,
                                            output_figure,
                                            bin_width = 3,
                                            conf_level = 0.95) {
  
  # 1) bin + clean
  data_b <- NDVI_TDiffbin(data, bin_width = bin_width) %>%
    tidyr::drop_na()
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(car)
  library(scales)
  library(latex2exp)
  
  value_col <- "avg_value"
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  
  data_clean <- data_b %>%
    mutate(
      species = factor(species, levels = species_levels),
      x = bin_mean
    ) %>%
    filter(!is.na(.data[[value_col]]), is.finite(x))
  
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  threshold <- 11.5
  alpha_level <- 1 - conf_level
  
  # -------------------------
  # 2) Fit linear + exponential, choose best by AIC
  # exponential: f(x)=a + b*exp(-c*x)
  # -------------------------
  start_nls <- list(a = 5, b = 7, c = 0.04)
  ctrl_nls  <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  fits <- purrr::set_names(purrr::map(species_levels, \(sp) {
    df_sp <- filter(data_clean, species == sp)
    
    res <- list(lm = NULL, nls = NULL,
                AIC_linear = NA_real_, AIC_exponential = NA_real_,
                best = NA_character_, species = sp)
    
    if (nrow(df_sp) >= 2) {
      res$lm <- tryCatch(lm(avg_value ~ x, data = df_sp), error = \(e) NULL)
      if (!is.null(res$lm)) res$AIC_linear <- AIC(res$lm)
    }
    
    if (nrow(df_sp) >= 5) {
      res$nls <- tryCatch(
        nls(avg_value ~ a + b * exp(-c * x),
            data = df_sp,
            start = start_nls,
            control = ctrl_nls),
        error = \(e) NULL
      )
      if (!is.null(res$nls)) res$AIC_exponential <- AIC(res$nls)
    }
    
    aics <- c(linear = res$AIC_linear, exponential = res$AIC_exponential)
    res$best <- if (all(is.na(aics))) NA_character_ else names(which.min(aics))
    res
  }), species_levels)
  
  # -------------------------
  # 3) Panel (a): predictions + SYMMETRICAL CI (F2 style)
  #    IMPORTANT: use observed x values (unique bin_mean), not a seq grid
  # -------------------------
  pred_ci <- purrr::map_dfr(species_levels, \(sp) {
    sp_data <- filter(data_clean, species == sp)
    if (nrow(sp_data) == 0) return(tibble())
    
    x_seq <- sp_data %>%
      pull(x) %>%
      unique() %>%
      sort()
    
    f <- fits[[sp]]
    if (is.na(f$best)) return(tibble())
    
    if (f$best == "linear") {
      pr   <- predict(f$lm, newdata = data.frame(x = x_seq), se.fit = TRUE)
      crit <- qt(1 - alpha_level / 2, df = df.residual(f$lm))
      
      tibble(
        species = sp,
        x = x_seq,
        fit = as.numeric(pr$fit),
        lwr = as.numeric(pr$fit - crit * pr$se.fit),
        upr = as.numeric(pr$fit + crit * pr$se.fit)
      )
    } else {
      if (is.null(f$nls)) return(tibble())
      
      ci_points <- lapply(x_seq, function(val) {
        dm <- car::deltaMethod(f$nls, paste0("a + b * exp(-c * ", val, ")"),
                               parameterNames = c("a","b","c"))
        data.frame(x = val, fit = dm$Estimate, se = dm$SE)
      }) %>% dplyr::bind_rows()
      
      crit <- qt(1 - alpha_level / 2, df = summary(f$nls)$df[2])
      
      tibble(
        species = sp,
        x = ci_points$x,
        fit = ci_points$fit,
        lwr = ci_points$fit - (crit * ci_points$se),
        upr = ci_points$fit + (crit * ci_points$se)
      )
    }
  }) %>%
    mutate(species = factor(species, levels = species_levels))
  
  # Panel (a) plot (match F2)
  p_line <- ggplot() +
    geom_point(
      data = data_clean,
      aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
      alpha = 0.4
    ) +
    geom_ribbon(
      data = pred_ci,
      aes(x = x, ymin = lwr, ymax = upr, fill = species),
      alpha = 0.25,
      show.legend = FALSE
    ) +
    geom_line(
      data = pred_ci,
      aes(x = x, y = fit, color = species),
      linewidth = 1.2,
      key_glyph = ribbon_line_key
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", alpha = 0.5) +
    annotate(
      "text",
      x = min(data_clean$x, na.rm = TRUE),
      y = threshold,
      label = "median",
      hjust = -7,
      vjust = -0.3,
      size = 6
    ) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette,  name = "") +
    scale_shape_manual(values = c("Oak"=16,"Beech"=17,"Spruce"=15,"Pine"=18), name = "") +
    scale_size_continuous(range = c(1, 8), labels = scales::percent_format(), name = "") +
    labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(color = "grey80", shape = 16, alpha = 1)
      ),
      color = guide_legend(
        order = 2,
        override.aes = list(shape = NA, linewidth = 1.2)
      ),
      fill  = "none",
      shape = "none"
    ) +
    base_theme +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.1, vjust = 1)
    )
  
  # -------------------------
  # 4) Panels (b) & (c): stats from chosen best model (F2 style)
  # -------------------------
  one_species_stats <- function(sp) {
    f <- fits[[sp]]
    if (is.na(f$best)) return(NULL)
    sp_data <- filter(data_clean, species == sp)
    
    if (f$best == "linear") {
      summ <- summary(f$lm)$coefficients
      a <- summ["(Intercept)", "Estimate"]
      b <- summ["x", "Estimate"]
      se_b <- summ["x", "Std. Error"]
      p_b  <- summ["x", "Pr(>|t|)"]
      r2   <- summary(f$lm)$r.squared
      x50  <- (threshold - a) / b
      
      tibble(species = sp, x50 = x50, abs_slope = abs(b),
             r_squared = r2, se = se_b, p_val = p_b)
    } else {
      coefs <- coef(f$nls); a <- coefs["a"]; b <- coefs["b"]; c <- coefs["c"]
      
      x50 <- ifelse((threshold - a) > 0 & b > 0, -log((threshold - a)/b)/c, NA_real_)
      slope50 <- c * (a - threshold) # == -c*(threshold-a)
      abs_slope <- abs(slope50)
      
      fitted_vals <- predict(f$nls, newdata = sp_data)
      r2 <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) /
        sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
      
      dm <- tryCatch(
        car::deltaMethod(f$nls, paste0("c*(a - ", threshold, ")"),
                         parameterNames = c("a","b","c")),
        error = \(e) NULL
      )
      se <- if (is.null(dm)) NA_real_ else dm$SE
      df_exp <- tryCatch(summary(f$nls)$df[2], error = \(e) NA_real_)
      p_val <- if (is.finite(se) && !is.na(df_exp))
        2 * (1 - stats::pt(abs(slope50 / se), df_exp)) else NA_real_
      
      tibble(species = sp, x50 = x50, abs_slope = abs_slope,
             r_squared = r2, se = se, p_val = p_val)
    }
  }
  
  stats_all <- bind_rows(lapply(species_levels, one_species_stats)) %>%
    mutate(species = factor(species, levels = species_levels))
  
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = TeX("$T_{d{,}Q_{50}}$ (mm)")) +
    base_theme + theme(legend.position = "none") +
    ggtitle("(b)") +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.25, vjust = 1)
    )
  
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")),
                  y = abs_slope/2)) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = TeX("$|\\lambda_{\\, \\, Q_{50}}(T_{d})|$")) +
    ggtitle("(c)") +
    base_theme +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(size = 18, hjust = 0.25, vjust = 1)
    )
  
  # -------------------------
  # 5) Final layout + save (match F2)
  # -------------------------
  final_plot <- (p_line + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(filename = output_figure, plot = final_plot, width = 12, height = 8, dpi = 300)
  
  invisible(final_plot)
}


plot_NDVI_TDiff_exp_slope_coeff(
  data,
  "results_rootzone/Figures_till2022/SI_PSI/SI_NDVI_Q_TDiffbin_exp_slope.png"
)

