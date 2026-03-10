setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)
library(patchwork)
library(mgcv)
library(scales)
library(purrr)

data <- combined %>%
  filter(month %in% c("July", "August")) %>%
  filter(year < 2023) %>%
  filter(Quantiles > 0) %>%
  na.omit()

# ---- theme ----
base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 14),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title        = element_text(size = 16),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_blank(),
    strip.text        = element_text(size = 14)
  )

TDiff_PSIbin <- function(df, bin_width = 50, min_count = 1000) {
  
  # identify transpiration deficit column
  value_column <- if ("transpiration_deficit" %in% names(df)) {
    "transpiration_deficit"
  } else {
    stop("Column 'transpiration_deficit' not found in the input data frame.")
  }
  
  # total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # define PSI bin breaks
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE) / bin_width) * bin_width
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE) / bin_width) * bin_width
  bin_breaks <- seq(psi_min, psi_max + bin_width, by = bin_width)
  
  # assign bins
  df_binned <- df %>%
    mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # calculate mean PSI and mean TDiff per species × bin
  meanTDiff_PSIbin_species <- df_binned %>%
    group_by(species, PSI_bin) %>%
    summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      avg_transpiration_deficit = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = "drop"
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    select(
      species,
      PSI_bin,
      bin_mean,
      avg_transpiration_deficit,
      count,
      total_pixels,
      percentage
    )
  
  return(meanTDiff_PSIbin_species)
}

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

plot_TDiff_PSIbin_onlyA <- function(data,
                                    figure_output,
                                    bin_width = 50,
                                    min_count = 1000) {
  
  library(dplyr)
  library(ggplot2)
  library(mgcv)
  library(scales)
  library(purrr)
  library(tibble)
  
  # prepare binned data
  df <- TDiff_PSIbin(data, bin_width = bin_width, min_count = min_count) %>%
    drop_na()
  
  # species order and palette
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  
  df <- df %>%
    mutate(species = factor(species, levels = species_levels))
  
  xcol <- "bin_mean"
  ycol <- "avg_transpiration_deficit"
  
  # fit GAM for each species and print chosen k
  fit_one_species <- function(sp) {
    sdf <- df %>%
      filter(species == sp) %>%
      arrange(.data[[xcol]])
    
    n_bins <- nrow(sdf)
    
    if (n_bins < 5) {
      cat("Species:", sp, "| bins:", n_bins, "| skipped (not enough bins for GAM)\n")
      return(NULL)
    }
    
    k_use <- min(6, max(3, floor(n_bins / 2)))
    
    cat("Species:", sp, "| bins:", n_bins, "| k used:", k_use, "\n")
    
    gam_fit <- tryCatch(
      mgcv::gam(
        reformulate(sprintf("s(%s, k = %d)", xcol, k_use), response = ycol),
        data = sdf,
        method = "REML"
      ),
      error = function(e) {
        cat("Species:", sp, "| GAM fitting failed:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    
    if (is.null(gam_fit)) return(NULL)
    
    xseq <- seq(
      min(sdf[[xcol]], na.rm = TRUE),
      max(sdf[[xcol]], na.rm = TRUE),
      length.out = 200
    )
    
    pred_df <- tibble(
      species = factor(sp, levels = species_levels),
      bin_mean = xseq
    )
    
    pred <- predict(gam_fit, newdata = pred_df, se.fit = TRUE)
    
    pred_df %>%
      mutate(
        predicted = pred$fit,
        se = pred$se.fit,
        lower = predicted - 1.96 * se,
        upper = predicted + 1.96 * se
      )
  }
  
  pred_df <- purrr::map_dfr(species_levels, fit_one_species)
  
  # final plot
  plot_clean <- ggplot() +
    geom_point(
      data = df,
      aes(
        x = .data[[xcol]],
        y = .data[[ycol]],
        color = species,
        shape = species,
        size = percentage
      ),
      alpha = 0.7
    ) +
    geom_ribbon(
      data = pred_df,
      aes(
        x = bin_mean,
        ymin = lower,
        ymax = upper,
        fill = species,
        group = species
      ),
      alpha = 0.25,
      show.legend = FALSE
    ) +
    geom_line(
      data = pred_df,
      aes(
        x = bin_mean,
        y = predicted,
        color = species
      ),
      linewidth = 1.2,
      key_glyph = ribbon_line_key
    ) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, name = "") +
    scale_shape_manual(
      values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
      name = ""
    ) +
    scale_size_continuous(
      range = c(1, 8),
      labels = scales::percent_format(accuracy = 1),
      name = ""
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit (mm)"
    ) +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          color = "grey80",
          shape = 16,
          alpha = 1
        )
      ),
      color = guide_legend(
        order = 2,
        override.aes = list(
          shape = NA,
          linewidth = 1.2
        )
      ),
      fill  = "none",
      shape = "none"
    ) +
    base_theme +
    theme(
      legend.position = "bottom"
    )
  
  print(plot_clean)
  
  if (!missing(figure_output) && !is.null(figure_output)) {
    ggsave(
      filename = figure_output,
      plot = plot_clean,
      width = 12,
      height = 8,
      dpi = 300
    )
  }
  
  invisible(plot_clean)
}

plot_TDiff_PSIbin_onlyA(
  data = data,
  figure_output = "results_rootzone/Figures_till2022/SI_PSI/SI_TDiff_PSIbin_onlyA_till2022.png",
  bin_width = 50,
  min_count = 1000
)
