# ─── 0. Setup ────────────────────────────────────────────────────────────────
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0001/yixuan/NDVI_PSI_project")

# Libraries
library(terra)
library(tidyverse)

prepare_tdiff_summary_data <- function(
    years       = 2003:2022,
    months      = c("July", "August"),
    species     = c("Oak", "Beech", "Spruce", "Pine"),
    base_dir    = "results_monthly_rootzone",
    file_prefix = "tdiff_",
    bin_width   = 1,        # mm bins; change as needed
    min_count   = 1000,
    output_path = "results_rootzone/Data/tdiff_summary_pairwise.RData"
) {
  suppressPackageStartupMessages({
    library(terra)
    library(tidyverse)
  })
  
  # 1) Helper: Read rasters
  read_tdiff <- function(year, month, species_name) {
    file_path <- file.path(base_dir, month, species_name, paste0(file_prefix, year, ".tif"))
    if (!file.exists(file_path)) return(NULL)
    
    r <- rast(file_path)
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    colnames(df)[3] <- "tdiff_value"
    
    df %>%
      mutate(
        species = species_name,
        month   = month,
        year    = year
      )
  }
  
  message("Reading rasters...")
  tdiff_all <- expand.grid(
    year = years, month = months, species = species, stringsAsFactors = FALSE
  ) %>%
    pmap_dfr(~ read_tdiff(..1, ..2, ..3))
  
  if (nrow(tdiff_all) == 0) {
    stop("No tdiff rasters were found. Check base_dir/month/species/subdir_name/file_prefix.")
  }
  
  # 2) Pivot Wide for Pairwise Processing
  message("Processing pairwise joins...")
  tdiff_wide <- tdiff_all %>%
    pivot_wider(names_from = species, values_from = tdiff_value)
  
  # 3) Long + Pairwise Join
  tdiff_long <- tdiff_wide %>%
    pivot_longer(cols = all_of(species), names_to = "species", values_to = "tdiff_value")
  
  tdiff_pairs <- tdiff_long %>%
    inner_join(
      tdiff_long,
      by = c("x", "y", "month", "year"),
      suffix = c("_x", "_y"),
      relationship = "many-to-many"
    ) %>%
    filter(species_x != species_y) %>%
    mutate(
      species_x = factor(species_x, levels = species),
      species_y = factor(species_y, levels = species)
    )
  
  # 4) Binning and Summarising
  message("Binning data...")
  tdiff_summary <- tdiff_pairs %>%
    group_by(species_x) %>%
    mutate(
      TDIFF_bin = cut(
        tdiff_value_x,
        breaks = seq(
          floor(min(tdiff_value_x, na.rm = TRUE)),
          ceiling(max(tdiff_value_x, na.rm = TRUE)),
          by = bin_width
        ),
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    group_by(species_x, species_y, TDIFF_bin) %>%
    summarise(
      bin_mean      = mean(tdiff_value_x, na.rm = TRUE),
      mean_tdiff_y  = mean(tdiff_value_y, na.rm = TRUE),
      count         = n(),
      .groups       = "drop"
    ) %>%
    filter(count >= min_count) %>%
    # Match your PSI definition: percentage within each (species_x, species_y)
    group_by(species_x, species_y) %>%
    mutate(percentage = count / sum(count, na.rm = TRUE)) %>%
    ungroup()
  
  # Save Results
  if (!dir.exists(dirname(output_path))) {
    dir.create(dirname(output_path), recursive = TRUE)
  }
  
  save(tdiff_summary, file = output_path)
  message(paste("Process Complete. Data saved to:", output_path))
  
  return(tdiff_summary)
}

prepare_tdiff_summary_data(years = 2003:2022)

load("results_rootzone/Data/tdiff_summary_pairwise.RData") 

# ─── 7. Plot: four-panel tdiff–tdiff comparison ──────────────────────────────────
plot_tdiff_pairwise_LM <- function(
    tdiff_summary,
    point_alpha = 0.4,
    band_alpha = 0.25,
    line_width = 1,
    abline_lwd = 0.8,
    cb_palette = c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442"),
    cb_shapes  = c(Oak=16, Beech=17, Spruce=15, Pine=18),
    save_path = NULL,
    save_width = 12,
    save_height = 8,
    save_dpi = 300
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
    library(dplyr)
  })
  
  # Ensure consistent order (same as PSI function)
  tdiff_summary <- tdiff_summary %>%
    mutate(
      species_x = factor(species_x, levels = names(cb_palette)),
      species_y = factor(species_y, levels = names(cb_palette))
    )
  
  p <- ggplot(
    tdiff_summary,
    aes(x = bin_mean, y = mean_tdiff_y)
  ) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", alpha = 0.5, linewidth = abline_lwd
    ) +
    
    # Points: show ONLY size legend (no species legend from points)
    geom_point(
      aes(color = species_y, shape = species_y, size = percentage),
      alpha = point_alpha,
      show.legend = c(color = FALSE, shape = FALSE, size = TRUE)
    ) +
    
    # Linear regression with confidence band (per species_y, within each facet)
    geom_smooth(
      aes(color = species_y, fill = species_y, group = species_y),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      alpha = band_alpha,
      linewidth = line_width
    ) +
    
    facet_wrap(~ species_x, ncol = 2, scales = "free") +
    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette,  name = "") +
    scale_shape_manual(values = cb_shapes,  name = "") +
    
    scale_size_continuous(
      name   = "",
      range  = c(1, 8),
      labels = percent_format(accuracy = 1)
    ) +
    
    labs(
      x = "transpiration deficit of panel species (mm)",
      y = "transpiration deficit of line species (mm)"
    ) +
    
    base_theme +
    
    # Legend control: size first (grey), then fill (species), no color/shape legends
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          colour = "grey80",
          shape = 16,
          fill = NA,
          linewidth = 0,
          alpha = 1
        )
      ),
      fill = guide_legend(
        order = 2,
        override.aes = list(
          linewidth = line_width,
          alpha = 0.4,
          colour = unname(cb_palette),
          shape = NA,
          size = 3
        )
      ),
      color = "none",
      shape = "none"
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.background = element_blank()
    )
  
  print(p)
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(save_path, p, width = save_width, height = save_height, dpi = save_dpi, bg = "white")
  }
  
  invisible(p)
}

plot_tdiff_pairwise_LM(
  tdiff_summary, # Input the data you created/loaded
  save_path = "results_rootzone/Figures_till2022/SI_PSI/SI_TDiff_TDiff_pairwise_4panel_LM.png"
)