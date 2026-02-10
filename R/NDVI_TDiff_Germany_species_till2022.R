setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)

########### only species located data #######################
data <- combined %>% 
  filter(month %in% c("July", "August"),
         Quantiles > 0,
         year <= 2022) %>%
  drop_na()

# -------------------------
# 1) BINNING: use bin_mean = mean(x in each bin)
# -------------------------
NDVI_TDiffbin <- function(df, bin_width = 3, min_count = 1000) {
  
  value_col <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df %>%
    as_tibble() %>%
    mutate(
      TDiff_bin = cut(transpiration_deficit,
                      breaks = bin_breaks,
                      include.lowest = TRUE,
                      right = FALSE)
    ) %>%
    group_by(species, TDiff_bin) %>%
    summarise(
      bin_mean  = mean(transpiration_deficit, na.rm = TRUE),  # ✅ NEW
      avg_value = mean(.data[[value_col]], na.rm = TRUE),
      count     = n(),
      .groups   = "drop"
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    select(species, TDiff_bin, bin_mean, avg_value, count, total_pixels, percentage)
}

# -------------------------
# 2) Delta-method CI for exponential model: f(x)=a + b*exp(-c*x)
# -------------------------
pred_exp_delta_TDiff <- function(fit, xseq, conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  zcrit <- qnorm(1 - alpha/2)
  
  co <- coef(fit)
  V  <- vcov(fit)
  
  a <- unname(co["a"]); b <- unname(co["b"]); c <- unname(co["c"])
  ex <- exp(-c * xseq)
  
  pred <- a + b * ex
  
  # gradient wrt (a,b,c):
  # df/da = 1
  # df/db = exp(-c x)
  # df/dc = b * d/dc exp(-c x) = b * (-x) * exp(-c x)
  se <- map_dbl(seq_along(xseq), \(i) {
    g <- c(1, ex[i], b * (-xseq[i]) * ex[i])
    as.numeric(sqrt(t(g) %*% V %*% g))
  })
  
  tibble(
    x = xseq,
    pred = pred,
    lower = pred - zcrit * se,
    upper = pred + zcrit * se
  )
}

# -------------------------
# 3) Main plot function (bin_mean + CI ribbon in panel a)
# -------------------------
plot_NDVI_TDiff_exp_slope_coeff <- function(data,
                                            combined_coef_fig,
                                            output_figure,
                                            aic_barplot_fig,
                                            bin_width = 3,
                                            conf_level = 0.95) {
  
  # bin + clean
  data_b <- NDVI_TDiffbin(data, bin_width = bin_width) %>%
    drop_na()
  
  library(ggplot2)
  library(patchwork)
  library(car)
  
  value_col <- "avg_value"
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  
  data_clean <- data_b %>%
    mutate(
      species = factor(species, levels = species_levels),
      x = bin_mean                     # ✅ NEW: use bin_mean for x
    ) %>%
    filter(!is.na(.data[[value_col]]), is.finite(x))
  
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  threshold <- 11.5
  
  # -------------------------
  # AIC comparison per species
  # -------------------------
  start_nls <- list(a = 5, b = 7, c = 0.04)  # f(x)=a + b*exp(-c*x)
  ctrl_nls  <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  fits <- set_names(map(species_levels, \(sp) {
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
  
  aic_df <- map_dfr(fits, \(f) {
    tibble(species = f$species,
           AIC_linear = f$AIC_linear,
           AIC_exponential = f$AIC_exponential)
  }) %>% mutate(species = factor(species, levels = species_levels))
  
  aic_long <- aic_df %>%
    pivot_longer(c(AIC_linear, AIC_exponential), names_to = "model", values_to = "AIC") %>%
    mutate(model = dplyr::recode(model, AIC_linear = "linear", AIC_exponential = "exponential"))
  
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
  
  # -------------------------
  # Coefficient plots (unchanged style)
  # -------------------------
  chosen_linear <- names(Filter(\(f) identical(f$best, "linear"), fits))
  chosen_exp    <- names(Filter(\(f) identical(f$best, "exponential"), fits))
  
  lin_coef_df <- if (length(chosen_linear) > 0) {
    map_dfr(chosen_linear, \(sp) {
      mod <- fits[[sp]]$lm
      if (is.null(mod)) return(tibble())
      summ <- summary(mod)$coefficients
      tibble(
        species     = sp,
        Coefficient = dplyr::recode(rownames(summ), "(Intercept)" = "a", "x" = "b"),
        Value       = summ[, "Estimate"],
        pvalue      = summ[, "Pr(>|t|)"]
      )
    }) %>%
      mutate(species = factor(species, levels = species_levels),
             label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  } else tibble()
  
  p_coeff_linear <- if (nrow(lin_coef_df) > 0) {
    ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
      geom_col(position = position_dodge(0.9)) +
      geom_text(aes(label = label, y = Value/2), position = position_dodge(0.9), vjust = 0.5) +
      scale_fill_manual(values = cb_palette) +
      facet_wrap(~ Coefficient, scales = "free_y") +
      labs(title = "linear models",
           subtitle = expression(NDVI == a + b * x + epsilon),
           x = NULL, y = "coefficient value") +
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
  } else ggplot() + theme_void()
  
  exp_coef_df <- if (length(chosen_exp) > 0) {
    map_dfr(chosen_exp, \(sp) {
      mod <- fits[[sp]]$nls
      if (is.null(mod)) return(tibble())
      summ <- summary(mod)$coefficients
      tibble(
        species     = sp,
        Coefficient = tolower(rownames(summ)),
        Value       = summ[, "Estimate"],
        pvalue      = summ[, "Pr(>|t|)"]
      )
    }) %>%
      mutate(species = factor(species, levels = species_levels),
             label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  } else tibble()
  
  p_coeff_exp <- if (nrow(exp_coef_df) > 0) {
    ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
      geom_col(position = position_dodge(0.9)) +
      geom_text(aes(label = label, y = Value/2), position = position_dodge(0.9), vjust = 0.5) +
      scale_fill_manual(values = cb_palette) +
      facet_wrap(~ Coefficient, scales = "free_y") +
      labs(title = "exponential models",
           subtitle = expression(NDVI == a + b * e^{-~c * x} + epsilon),
           x = NULL, y = "coefficient value") +
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
  } else ggplot() + theme_void()
  
  # combined_coeff <- if (inherits(p_coeff_linear, "ggplot") && inherits(p_coeff_exp, "ggplot")) {
  #   p_coeff_linear + p_coeff_exp
  # } else p_coeff_exp
  
  combined_coeff <- p_coeff_exp
   
  print(combined_coeff)
  ggsave(combined_coef_fig, combined_coeff, width = 10, height = 8, dpi = 300)
  
  # -------------------------
  # Panel (a): predictions + confidence bands for BEST model
  # -------------------------
  pred_ci <- map_dfr(species_levels, \(sp) {
    sp_data <- filter(data_clean, species == sp)
    if (nrow(sp_data) == 0) return(tibble())
    
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 200)
    f <- fits[[sp]]
    
    if (is.na(f$best)) return(tibble())
    
    if (f$best == "linear") {
      pr <- predict(f$lm, newdata = data.frame(x = x_seq),
                    interval = "confidence", level = conf_level)
      
      tibble(species = sp, x = x_seq,
             pred = pr[, "fit"], lower = pr[, "lwr"], upper = pr[, "upr"])
    } else {
      if (is.null(f$nls)) return(tibble())
      pred_exp_delta_TDiff(f$nls, x_seq, conf_level = conf_level) %>%
        mutate(species = sp)
    }
  }) %>%
    mutate(species = factor(species, levels = species_levels))
  
  p_combined <- ggplot() +
    # ✅ confidence ribbon behind
    geom_ribbon(
      data = pred_ci,
      aes(x = x, ymin = lower, ymax = upper, fill = species),
      alpha = 0.18,
      color = NA
    ) +
    geom_point(
      data = data_clean,
      aes(x = x, y = avg_value, color = species, shape = species, size = percentage),
      alpha = 0.7
    ) +
    geom_line(
      data = pred_ci,
      aes(x = x, y = pred, color = species),
      linewidth = 1
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", linewidth = 1) +
    annotate("text", x = 32, y = threshold,
             label = "median", fontface = "italic",
             hjust = 0.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    scale_shape_manual(values = c("Oak"=16,"Beech"=17,"Spruce"=15,"Pine"=18)) +
    scale_size_continuous(name = "Pixels per bin (%)",
                          range = c(1, 8),
                          labels = scales::percent_format(accuracy = 1)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 1),
           size = guide_legend(order = 2)) +
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
  
  # -------------------------
  # Panels (b) and (c): stats from chosen best model (same as your current logic)
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
      slope50 <- -c * (threshold - a)
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
        2*(1-pt(abs(slope50/se), df_exp)) else NA_real_
      
      tibble(species = sp, x50 = x50, abs_slope = abs_slope,
             r_squared = r2, se = se, p_val = p_val)
    }
  }
  
  stats_all <- bind_rows(lapply(species_levels, one_species_stats)) %>%
    mutate(species = factor(species, levels = species_levels))
  
  p_x50 <- ggplot(stats_all, aes(x=species, y=x50, fill=species)) +
    geom_col(width=0.7) +
    scale_fill_manual(values=cb_palette, guide=FALSE) +
    labs(y="transpiration deficit (mm)", x="") +
    ggtitle("(b)") +
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
  
  p_slope <- ggplot(stats_all, aes(x=species, y=abs_slope, fill=species)) +
    geom_col(width=0.7) +
    geom_errorbar(aes(ymin=pmax(0, abs_slope-se), ymax=abs_slope+se), width=0.2) +
    geom_text(aes(label=if_else(p_val<0.05, sprintf("%.2f*", r_squared), sprintf("%.2f", r_squared)),
                  y=abs_slope/2), size=5) +
    scale_fill_manual(values=cb_palette, guide=FALSE) +
    labs(y="absolute slope", x="") +
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
  
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths=c(2,1), guides="collect") &
    theme(legend.position="bottom", legend.title=element_blank())
  
  print(final_plot)
  
  # plateau for best exponential only
  plateau_df <- bind_rows(lapply(chosen_exp, function(sp) {
    mod <- fits[[sp]]$nls
    if (is.null(mod)) return(NULL)
    c_val <- coef(mod)["c"]
    tibble(species = sp, x_plateau = -log(0.05) / c_val)
  }))
  
  print("Estimated transpiration deficit (x_plateau) where NDVI is ~95% of its asymptotic minimum (chosen exponential models only):")
  print(plateau_df)
  
  ggsave(combined_coef_fig, combined_coeff, width=10, height=8, dpi=300)
  ggsave(output_figure, final_plot, width=10, height=8, dpi=300)
}

# Run
plot_NDVI_TDiff_exp_slope_coeff(
  data,
  "results_rootzone/test/NDVI_Q_TDiffbin_exp_coeff.png",
  "results_rootzone/test/NDVI_Q_TDiffbin_exp_slope.png",
  "results_rootzone/test/NDVI_Q_TDiffbin_exp_aic.png",
  bin_width = 3,
  conf_level = 0.95
)

NDVI_TDiffbin <- function(df, bin_width = 3, min_count = 1000) {
  
  # totals per (area, species) for percentage (same logic as PSI version)
  species_totals <- df %>%
    dplyr::group_by(area, species) %>%
    dplyr::summarise(total_pixels = dplyr::n(), .groups = "drop")
  
  # global bin breaks so bins consistent across everything
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df %>%
    dplyr::mutate(
      TDiff_bin = cut(
        transpiration_deficit,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    dplyr::group_by(area, species, TDiff_bin) %>%
    dplyr::summarise(
      bin_mean  = mean(transpiration_deficit, na.rm = TRUE),
      avg_value = mean(NDVI, na.rm = TRUE),
      count     = dplyr::n(),
      .groups   = "drop"
    ) %>%
    dplyr::left_join(species_totals, by = c("area", "species")) %>%
    dplyr::mutate(percentage = count / total_pixels) %>%
    dplyr::filter(count >= min_count) %>%
    dplyr::select(area, species, TDiff_bin, bin_mean, avg_value, count, total_pixels, percentage)
}

### loess ###
plot_NDVI_Q_TDiffbin_area <- function(
    df,
    out_file,
    bin_width = 3,
    min_count = 1000,
    width = 12,
    height = 5,
    dpi = 300
) {
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  facet_labs <- c(
    "Germany" = "(a) Germany",
    "species" = "(b) species"
  )
  
  p <- NDVI_TDiffbin(df, bin_width = bin_width, min_count = min_count) %>%
    dplyr::mutate(area = factor(area, levels = c("Germany", "species"))) %>%
    ggplot2::ggplot(ggplot2::aes(bin_mean, avg_value)) +
    ggplot2::geom_point(ggplot2::aes(color = species, size = percentage), alpha = 0.8) +
    ggplot2::geom_smooth(
      ggplot2::aes(color = species, fill = species),
      method = "loess",
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~ area, nrow = 1, labeller = ggplot2::as_labeller(facet_labs)) +
    ggplot2::scale_color_manual(values = cb_palette, name = "") +
    ggplot2::scale_fill_manual(values = cb_palette, guide = "none") +
    ggplot2::scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title  = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      strip.text = ggplot2::element_text(size = 16, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
}


# ---- run ----
plot_NDVI_Q_TDiffbin_area(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_Q_TDiffbin_area_loess.png"
)

### GAM ###
plot_NDVI_Q_TDiffbin_area <- function(
    df,
    out_file,
    bin_width = 3,
    min_count = 1000,
    k = 6,                 # GAM flexibility (increase for more wiggle)
    bs = "cs",             # shrinkage spline
    width = 12,
    height = 5,
    dpi = 300
) {
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  facet_labs <- c(
    "Germany" = "(a) Germany",
    "species" = "(b) species"
  )
  
  p <- NDVI_TDiffbin(df, bin_width = bin_width, min_count = min_count) %>%
    dplyr::mutate(area = factor(area, levels = c("Germany", "species"))) %>%
    ggplot2::ggplot(ggplot2::aes(bin_mean, avg_value)) +
    ggplot2::geom_point(ggplot2::aes(color = species, size = percentage), alpha = 0.8) +
    ggplot2::geom_smooth(
      ggplot2::aes(color = species, fill = species),
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~ area, nrow = 1, labeller = ggplot2::as_labeller(facet_labs)) +
    ggplot2::scale_color_manual(values = cb_palette, name = "") +
    ggplot2::scale_fill_manual(values = cb_palette, guide = "none") +
    ggplot2::scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title  = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      strip.text = ggplot2::element_text(size = 16, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
}


# ---- run ----
plot_NDVI_Q_TDiffbin_area(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_Q_TDiffbin_area_gam.png",
  k = 6,
  bs = "cs"
)

# --- NDVI vs TDiff: 4 species facets + 2 area curves (Germany vs species), GAM ---

plot_NDVI_Q_TDiffbin_species4_area2 <- function(
    df,
    out_file,
    bin_width = 3,
    min_count = 1000,
    k = 6,
    bs = "cs",
    width = 12,
    height = 8,
    dpi = 300
) {
  
  # species order + colors for the two "area" curves
  sp_levels <- c("Oak", "Beech", "Spruce", "Pine")
  area_cols <- c("Germany" = "blue", "species" = "red")
  
  p <- NDVI_TDiffbin(df, bin_width = bin_width, min_count = min_count) %>%
    dplyr::mutate(
      area    = factor(area, levels = c("Germany", "species")),
      species = factor(species, levels = sp_levels)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(bin_mean, avg_value)) +
    
    # points (colored + shaped by area)
    ggplot2::geom_point(
      ggplot2::aes(color = area, shape = area, size = percentage),
      alpha = 0.8
    ) +
    
    # GAM smooths (two lines per facet), ribbon follows area color
    ggplot2::geom_smooth(
      ggplot2::aes(color = area, fill = area),
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.20,
      linewidth = 1
    ) +
    
    # 4 panels, ordered Oak -> Beech -> Spruce -> Pine
    ggplot2::facet_wrap(~ species, nrow = 2) +
    
    # blue/red for areas
    ggplot2::scale_color_manual(values = area_cols, name = "") +
    ggplot2::scale_fill_manual(values = area_cols, guide = "none") +
    
    # shapes for the two lines/points
    ggplot2::scale_shape_manual(values = c("Germany" = 16, "species" = 17), name = "") +
    
    ggplot2::scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    ggplot2::labs(
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)"
    ) +
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title  = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      strip.text = ggplot2::element_text(size = 16, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
  invisible(p)
}

# ---- run ----
plot_NDVI_Q_TDiffbin_species4_area2(
  df = DT_combined,
  out_file = "results_rootzone/test/NDVI_Q_TDiffbin_4species_2areas_gam.png",
  bin_width = 3,
  min_count = 1000,
  k = 5,
  bs = "cs"
)

# ---------- 2) NDVI difference along TDiff bins ----------
NDVI_diff_species_minus_Germany_TDiffbin <- function(df, bin_width = 3, min_count = 1000) {
  
  # common TDiff bin breaks (global so Germany/species use identical bins)
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df_binned <- df %>%
    dplyr::mutate(
      TDiff_bin = cut(
        transpiration_deficit,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # mean NDVI for species-mask (area == "species")
  sp_tbl <- df_binned %>%
    dplyr::filter(area == "species") %>%
    dplyr::group_by(species, TDiff_bin) %>%
    dplyr::summarise(
      bin_mean = mean(transpiration_deficit, na.rm = TRUE),
      sp_mean  = mean(NDVI, na.rm = TRUE),
      sp_n     = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::filter(sp_n >= min_count)
  
  # mean NDVI for Germany baseline (area == "Germany")
  de_tbl <- df_binned %>%
    dplyr::filter(area == "Germany") %>%
    dplyr::group_by(species, TDiff_bin) %>%
    dplyr::summarise(
      de_mean = mean(NDVI, na.rm = TRUE),
      de_n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(de_n >= min_count)
  
  # join + difference
  sp_tbl %>%
    dplyr::inner_join(de_tbl, by = c("species", "TDiff_bin")) %>%
    dplyr::mutate(diff_NDVI = sp_mean - de_mean) %>%
    dplyr::select(species, TDiff_bin, bin_mean, sp_mean, de_mean, sp_n, de_n, diff_NDVI)
}

plot_NDVI_diff_TDiff <- function(data, figure_output = NULL,
                                 bin_width = 3, min_count = 1000,
                                 k = 3, bs = "cs",
                                 width = 8, height = 6, dpi = 300) {
  
  df_diff <- NDVI_diff_species_minus_Germany_TDiffbin(data, bin_width, min_count) %>%
    dplyr::mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  p <- ggplot2::ggplot(df_diff,
                       ggplot2::aes(x = bin_mean, y = diff_NDVI, color = species)) +
    
    ggplot2::geom_hline(yintercept = 0,
                        linewidth = 0.6,
                        linetype = "dashed",
                        color = "grey40") +
    
    ggplot2::geom_point(alpha = 0.5, size = 2) +
    
    ggplot2::geom_smooth(
      ggplot2::aes(fill = species),
      method = "gam",
      formula = y ~ s(x, k = k, bs = bs),
      se = TRUE,
      alpha = 0.25,
      linewidth = 1
    ) +
    
    ggplot2::scale_color_manual(values = cb_palette, name = "") +
    ggplot2::scale_fill_manual(values = cb_palette, guide = "none") +
    
    ggplot2::labs(
      x = "transpiration deficit (TDiff)",
      y = expression(Delta[NDVI]~"(species - Germany)")
    ) +
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title  = ggplot2::element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  print(p)
  
  if (!is.null(figure_output)) {
    dir.create(dirname(figure_output), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(figure_output, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  invisible(p)
}

plot_NDVI_diff_TDiff(
  data = DT_combined,
  figure_output = "results_rootzone/test/NDVI_TDiffbin_area_diff_till2022.png",
  bin_width = 3,
  min_count = 1000,
  k = 3
)

### difference bar plot:
plot_NDVI_diff_TDiff_bar <- function(data, figure_output = NULL,
                                     bin_width = 3, min_count = 1000,
                                     width = 6, height = 5, dpi = 300) {
  
  library(dplyr)
  library(ggplot2)
  
  df_bar <- NDVI_diff_species_minus_Germany_TDiffbin(data, bin_width, min_count) %>%
    mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine"))) %>%
    group_by(species) %>%
    summarise(
      mean_diff = mean(diff_NDVI, na.rm = TRUE),
      se_diff   = sd(diff_NDVI, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  p <- ggplot(df_bar, aes(x = species, y = mean_diff, fill = species)) +
    
    geom_hline(yintercept = 0,
               linewidth = 0.6,
               linetype = "dashed",
               color = "grey40") +
    
    geom_col(width = 0.7, alpha = 0.9) +
    
    geom_errorbar(
      aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
      width = 0.2, linewidth = 0.6
    ) +
    
    scale_fill_manual(values = cb_palette, guide = "none") +
    
    labs(
      x = "",
      y = expression(mean~Delta[NDVI]~"(species - Germany)")
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  if (!is.null(figure_output)) {
    dir.create(dirname(figure_output), recursive = TRUE, showWarnings = FALSE)
    ggsave(figure_output, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  invisible(p)
}

plot_NDVI_diff_TDiff_bar(DT_combined, figure_output = "results_rootzone/test/NDVI_TDiffbin_area_diff_bar_till2022.png")


