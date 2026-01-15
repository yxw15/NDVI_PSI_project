library(tidyverse)
library(patchwork)
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)
library(grid) 

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

plot_TDiff_PSIbin_onlyA <- function(data, figure_output) {
  
  # -------------------------
  # Libraries & preprocessing
  # -------------------------
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(purrr)
  library(scales)
  
  # Prepare data
  df <- TDiff_PSIbin(data)
  df <- na.omit(df)
  
  # Species & palette
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_levels)
  
  xcol <- "bin_median"
  ycol <- "avg_transpiration_deficit"
  
  # -------------------------
  # Fit models & choose best degree by AIC
  # -------------------------
  degrees <- c(1, 2, 3)
  fit_one_species <- function(sp) {
    sdf <- df %>% filter(species == sp) %>% select(all_of(c(xcol, ycol)))
    if (nrow(sdf) < 5) {
      return(list(species = sp, best_deg = NA, best_fit = NULL))
    }
    fits <- lapply(degrees, function(d) {
      fm <- reformulate(termlabels = paste0("poly(", xcol, ", ", d, ")"), response = ycol)
      tryCatch(lm(fm, data = sdf), error = function(e) NULL)
    })
    aics <- sapply(fits, function(m) if (is.null(m)) Inf else AIC(m))
    best_idx <- which.min(aics)
    list(species = sp, best_deg = degrees[best_idx], best_fit = fits[[best_idx]])
  }
  fits_all <- lapply(species_levels, fit_one_species)
  
  # ---------------------------------
  # Predictions
  # ---------------------------------
  pred_list <- lapply(fits_all, function(ff) {
    if (is.null(ff$best_fit)) return(NULL)
    sp <- ff$species
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    xseq <- seq(xr[1], xr[2], length.out = 200)
    nd <- data.frame(bin_median = xseq)
    nd$species <- factor(sp, levels = species_levels)
    nd$predicted <- predict(ff$best_fit, newdata = nd)
    nd
  })
  pred_df <- bind_rows(pred_list)
  
  # -------------
  # Final plot: ONLY points + fitted lines
  # -------------
  plot_clean <- ggplot() +
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
    scale_color_manual(values = cb_palette, name = "") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "soil water potential (kPa)",
         y = "transpiration deficit") +
    theme_minimal() +
    theme(
      plot.title = element_blank(),        # remove title "(a)"
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid       = element_blank()
    )
  
  print(plot_clean)
  
  if (!missing(figure_output)) {
    ggsave(figure_output, plot = plot_clean, width = 8, height = 6, dpi = 300)
  }
  
  invisible(plot_clean)
}

plot_TDiff_PSIbin_onlyA(data, "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_onlyA_till2022.png")
