setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

data_all <- combined %>% filter(Quantiles > 0) %>% filter(month %in% c("May", "June", "July", "August"))
data_all <- na.omit(data_all)
data <- data_all %>% filter(month %in% c("July", "August")) %>% filter(year < 2023)

### S11 TDiff PSIbin ###
TDiff_PSIbin <- function(df, bin_width = 50) {
  
  library(dplyr)
  
  # Identify the correct column
  value_column <- if ("transpiration_deficit" %in% names(df)) "transpiration_deficit"
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks dynamically
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean transpiration_deficit and count per bin per species
  meanTDiff_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_transpiration_deficit = mean(.data[[value_column]], na.rm = TRUE),
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
    filter(percentage > 0.001) %>%
    select(species, PSI_bin, bin_median, avg_transpiration_deficit, count, total_pixels, percentage)
  
  return(meanTDiff_PSIbin_species)
}


plot_TDiff_PSIbin_slope <- function(data,
                                    coef_output,
                                    figure_output,
                                    aic_barplot_fig) {
  # -------------------------
  # Libraries & preprocessing
  # -------------------------
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(purrr)
  library(scales)  # percent_format
  
  # Prepare data
  df <- TDiff_PSIbin(data)
  df <- na.omit(df)
  
  # Species & palette (match your first function)
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_levels)
  
  # shorthand
  xcol <- "bin_median"
  ycol <- "avg_transpiration_deficit"
  
  # -------------------------
  # Fit models & choose by AIC
  # -------------------------
  degrees <- c(1, 2, 3)
  fit_one_species <- function(sp) {
    sdf <- df %>% filter(species == sp) %>% select(all_of(c(xcol, ycol)))
    if (nrow(sdf) < 5 || any(!is.finite(sdf[[xcol]])) || any(!is.finite(sdf[[ycol]]))) {
      return(list(species = sp,
                  aic_tbl = tibble(model = c("linear","poly2","poly3"), AIC = NA_real_),
                  best_deg = NA_integer_, best_fit = NULL))
    }
    fits <- lapply(degrees, function(d) {
      fm <- reformulate(termlabels = paste0("poly(", xcol, ", ", d, ")"), response = ycol)
      tryCatch(lm(fm, data = sdf), error = function(e) NULL)
    })
    aics <- sapply(fits, function(m) if (is.null(m)) NA_real_ else AIC(m))
    aic_tbl <- tibble(model = c("linear","poly2","poly3"), AIC = aics)
    best_idx <- suppressWarnings(which.min(aics))
    best_deg <- if (length(best_idx) == 0 || all(is.na(aics))) NA_integer_ else degrees[best_idx]
    best_fit <- if (is.na(best_deg)) NULL else fits[[best_idx]]
    list(species = sp, aic_tbl = aic_tbl, best_deg = best_deg, best_fit = best_fit)
  }
  fits_all <- lapply(species_levels, fit_one_species)
  
  # AIC dataframe (long) for barplot
  aic_df <- bind_rows(lapply(fits_all, function(x) {
    x$aic_tbl %>% mutate(species = x$species)
  }))
  aic_df$species <- factor(aic_df$species, levels = species_levels)
  aic_df$model <- factor(aic_df$model, levels = c("linear","poly2","poly3"))
  
  message("AIC values per species:")
  print(aic_df)
  
  # -------------------------
  # AIC barplot (saved)
  # -------------------------
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("linear" = "#0072B2",
                                 "poly2"  = "#E69F00",
                                 "poly3"  = "#009E73"),
                      name = "") +
    labs(x = "", y = "AIC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
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
  
  # ---------------------------------
  # Predictions from the best models
  # ---------------------------------
  pred_list <- lapply(fits_all, function(ff) {
    sp <- ff$species
    if (is.null(ff$best_fit) || is.na(ff$best_deg)) return(NULL)
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) return(NULL)
    xseq <- seq(xr[1], xr[2], length.out = 200)
    nd <- data.frame(bin_median = xseq)
    nd$species <- factor(sp, levels = species_levels)
    nd$predicted <- predict(ff$best_fit, newdata = nd)
    nd$degree <- ff$best_deg
    nd
  })
  pred_df <- bind_rows(pred_list)
  
  # Lines for mean / median based on observed y
  line_mean   <- mean(df[[ycol]], na.rm = TRUE)
  line_median <- median(df[[ycol]], na.rm = TRUE)
  
  print(line_median)
  
  # -------------
  # Panel (a): scatter + fitted lines + median reference
  # -------------
  plot_mixed <- ggplot() +
    geom_point(data = df,
               aes(x = .data[[xcol]],
                   y = .data[[ycol]],
                   color = species,
                   shape = species,
                   size = percentage),
               alpha = 0.7) +
    geom_line(data = pred_df,
              aes(x = bin_median, y = predicted, color = species),
              linewidth = 1) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text",
             x = -2000,
             y = line_median,
             label = "median",
             vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "Pixels per bin (%)",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "soil water potential (kPa)",
         y = "transpiration deficit") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold", size = 18, hjust = 0.1, vjust = 1),
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
  
  # -------------
  # Panel (b): soil water potential at median transpiration deficit
  #   -> **root finding only** with automatic bracketing (no interpolation fallback)
  # -------------
  solve_x_at_y <- function(fit, y_target, x_range, n_grid = 400) {
    # function for root finding
    f <- function(x) predict(fit, newdata = data.frame(bin_median = x)) - y_target
    
    # quick check at the endpoints
    fx_lo <- f(x_range[1]); fx_hi <- f(x_range[2])
    
    # build a grid to locate a sign change bracket
    xs <- seq(x_range[1], x_range[2], length.out = n_grid)
    fs <- suppressWarnings(f(xs))
    
    # find first sign change interval
    sc_idx <- which(diff(sign(fs)) != 0 & is.finite(fs[-length(fs)]) & is.finite(fs[-1]))
    if (length(sc_idx) == 0) {
      # no sign change -> no root in range
      return(NA_real_)
    }
    # choose the bracket where |f| near zero is smallest
    pick <- sc_idx[which.min(pmin(abs(fs[sc_idx]), abs(fs[sc_idx+1])))]
    lo <- xs[pick]; hi <- xs[pick + 1]
    
    # refine with uniroot
    tryCatch(uniroot(function(z) f(z), interval = c(lo, hi))$root,
             error = function(e) NA_real_)
  }
  
  # compute per species (root finding only)
  x_at_median <- lapply(fits_all, function(ff) {
    sp <- ff$species
    fit <- ff$best_fit
    if (is.null(fit)) return(tibble(species = sp, x_at_median = NA_real_))
    
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    xm <- solve_x_at_y(fit, line_median, xr, n_grid = 600)
    tibble(species = sp, x_at_median = xm)
  }) %>% bind_rows() %>% mutate(species = factor(species, levels = species_levels))
  
  p_bar_x <- ggplot(x_at_median, aes(x = species, y = x_at_median, fill = species)) +
    geom_bar(stat = "identity", width = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold", size = 18, hjust = 0.25, vjust = 1),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "none",
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # -------------
  # Panel (c): local slope near 1% of max(pred)
  # -------------
  local_slope_data <- pred_df %>%
    group_by(species) %>%
    group_modify(~{
      d <- .x
      maxp <- max(d$predicted, na.rm = TRUE)
      thresh <- 0.01 * maxp
      idx <- which.min(abs(d$predicted - thresh))
      w0 <- max(1, idx - 5)
      w1 <- min(nrow(d), idx + 5)
      dw <- d[w0:w1, , drop = FALSE]
      lm_loc <- lm(predicted ~ bin_median, data = dw)
      sm <- summary(lm_loc)
      data.frame(
        slope = coef(lm_loc)[["bin_median"]],
        slope_se = sm$coefficients["bin_median","Std. Error"],
        p_value = sm$coefficients["bin_median","Pr(>|t|)"],
        r2 = sm$r.squared
      )
    }) %>% ungroup() %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               sprintf("%.2f*", r2),
                               sprintf("%.2f", r2)))
  
  p_bar <- ggplot(local_slope_data, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - slope_se), ymax = slope_abs + slope_se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 5, color = "black") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold", size = 18, hjust = 0.25, vjust = 1),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "none",
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # -------------
  # Combine panels
  # -------------
  final_plot <- (plot_mixed + (p_bar_x / p_bar)) +
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
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  # -------------------------
  # Coefficients (best model per species)
  # -------------------------
  coef_rows <- lapply(fits_all, function(ff) {
    sp <- ff$species
    fit <- ff$best_fit
    deg <- ff$best_deg
    if (is.null(fit) || is.na(deg)) return(NULL)
    sm <- summary(fit)$coefficients
    out <- as.data.frame(sm)
    out$term_raw <- rownames(out)
    out$species <- sp
    out$degree <- deg
    out
  }) %>% bind_rows()
  
  if (!is.null(coef_rows) && nrow(coef_rows) > 0) {
    coef_data <- coef_rows %>%
      mutate(
        term = dplyr::case_when(
          term_raw == "(Intercept)" ~ "a",
          grepl("poly\\(bin_median,\\s*1\\)1", term_raw) ~ "b",
          grepl("poly\\(bin_median,\\s*2\\)1", term_raw) ~ "b",
          grepl("poly\\(bin_median,\\s*2\\)2", term_raw) ~ "c",
          grepl("poly\\(bin_median,\\s*3\\)1", term_raw) ~ "b",
          grepl("poly\\(bin_median,\\s*3\\)2", term_raw) ~ "c",
          grepl("poly\\(bin_median,\\s*3\\)3", term_raw) ~ "d",
          TRUE ~ term_raw
        ),
        species = factor(species, levels = species_levels)
      ) %>%
      dplyr::rename(value = Estimate, p_value = `Pr(>|t|)`) %>%
      mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)),
             term = factor(term, levels = c("a","b","c","d")))
    
    # Keep terms up to selected degree for each species
    coef_data <- coef_data %>%
      group_by(species) %>%
      filter(
        (term %in% c("a","b") & degree >= 1) |
          (term %in% c("c")    & degree >= 2) |
          (term %in% c("d")    & degree >= 3)
      ) %>% ungroup()
    
    plot_coeff <- ggplot(coef_data, aes(x = species, y = value, fill = species)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = label_text, y = value/2),
                color = "black", size = 5,
                position = position_dodge(width = 0.8)) +
      facet_wrap(~ term, scales = "free_y") +
      scale_fill_manual(values = cb_palette, name = "") +
      labs(title = "transpiration deficit (mm) ~ soil water potential (kPa)",
           x = "", y = "coefficient value") +
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
    
    print(plot_coeff)
    ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
    
    best_summary <- tibble(
      species = sapply(fits_all, `[[`, "species"),
      best_degree = sapply(fits_all, `[[`, "best_deg")
    )
    message("Best polynomial degree by AIC:")
    print(best_summary)
  } else {
    warning("No coefficients available to plot (fits failed or insufficient data).")
    plot_coeff <- NULL
  }
  
  print(list(combined_plot = final_plot,
               coeff_plot = plot_coeff,
               aic_plot = p_aic))

}

plot_TDiff_PSIbin_slope(
  data = data,
  coef_output = "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_coeff_till2022.png",
  figure_output = "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_slope_till2022.png",
  aic_barplot_fig = "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_AIC_till2022.png")
