setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Define
species <- c("Oak", "Beech", "Spruce", "Pine")
months <- c("July", "August")
years <- 2003:2024

# Load species masks
Oak    <- rast("species_map_MODIS/Oak.tif")
Beech  <- rast("species_map_MODIS/Beech.tif")
Spruce <- rast("species_map_MODIS/Spruce.tif")
Pine   <- rast("species_map_MODIS/Pine.tif")
species_masks <- list(Oak=Oak, Beech=Beech, Spruce=Spruce, Pine=Pine)

get_psi_path <- function(sp, month) {
  file.path("results_whole_Germany", paste0(sp, "_PSI_", month, ".tif"))
}

get_ndvi_path <- function(sp, month) {
  file.path("results_whole_Germany", paste0("NDVI_", month, ".tif"))
}

# Function to process one NDVI_species/PSI_species pair
process_ndvi_psi_pair <- function(ndvi_species, psi_species, months, years) {
  mask_r <- species_masks[[ndvi_species]]
  out_list <- list()
  for (month in months) {
    ndvi_file <- file.path("results_whole_Germany", paste0("NDVI_", month, ".tif"))
    ndvi_r <- rast(ndvi_file)
    # Project mask to NDVI if needed
    if (!compareGeom(ndvi_r, mask_r, stopOnError=FALSE)) {
      mask_r_proj <- project(mask_r, ndvi_r, method="near")
    } else {
      mask_r_proj <- mask_r
    }
    ndvi_r_masked <- mask(ndvi_r, mask_r_proj)
    psi_file <- get_psi_path(psi_species, month)
    psi_r <- rast(psi_file)
    # Project PSI to NDVI grid if needed
    if (!compareGeom(ndvi_r, psi_r, stopOnError=FALSE)) {
      psi_r_proj <- project(psi_r, ndvi_r)
    } else {
      psi_r_proj <- psi_r
    }
    psi_r_masked <- mask(psi_r_proj, mask_r_proj)
    n_layers <- min(nlyr(ndvi_r_masked), nlyr(psi_r_masked))
    for (lyr in seq_len(n_layers)) {
      ndvi_lyr <- ndvi_r_masked[[lyr]]
      psi_lyr  <- psi_r_masked[[lyr]]
      year     <- years[lyr]
      names(ndvi_lyr) <- "NDVI"
      names(psi_lyr)  <- "PSI"
      comb_r <- c(ndvi_lyr, psi_lyr)
      df <- as.data.frame(comb_r, xy=TRUE, na.rm=TRUE)
      df$NDVI_species <- ndvi_species
      df$PSI_species  <- psi_species
      df$Year         <- year
      df$Month        <- month
      out_list[[paste0(month, "_", year)]] <- df
    }
  }
  dplyr::bind_rows(out_list)
}

# Generate all pairs (exclude self-pairs if needed)
species_pairs <- expand.grid(
  NDVI_species = species,
  PSI_species  = species,
  stringsAsFactors = FALSE
)

# Process all pairs and combine
all_results <- list()
for (i in seq_len(nrow(species_pairs))) {
  ndvi_sp <- species_pairs$NDVI_species[i]
  psi_sp  <- species_pairs$PSI_species[i]
  cat("Processing NDVI:", ndvi_sp, "PSI:", psi_sp, "\n")
  all_results[[paste(ndvi_sp, psi_sp, sep = "_")]] <- process_ndvi_psi_pair(ndvi_sp, psi_sp, months, years)
}
big_df <- bind_rows(all_results)

write.csv(big_df, "results_whole_Germany/NDVI_PSI_crossspecies.csv", row.names = FALSE)


# Binning and summary - updated!
NDVI_PSIbin_multi <- function(df, bin_width = 50) {
  library(dplyr)
  
  df <- na.omit(df)
  psi_min <- floor(min(df$PSI, na.rm = TRUE))
  psi_max <- ceiling(max(df$PSI, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  df <- df %>%
    mutate(PSI_bin = cut(PSI, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # 1. Calculate total pixel count per NDVI_species and PSI_species
  total_pixel_df <- df %>%
    group_by(NDVI_species, PSI_species) %>%
    summarise(total_pixel = n(), .groups = "drop")
  
  # 2. Calculate mean NDVI, count, and median of each bin
  meanNDVI_PSIbin_species <- df %>%
    group_by(NDVI_species, PSI_species, PSI_bin) %>%
    summarise(
      avg_value = mean(NDVI, na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    # 3. Add the total pixel count
    left_join(total_pixel_df, by = c("NDVI_species", "PSI_species")) %>%
    # 4. Calculate pixel percentage
    mutate(
      pixel_percentage = count / total_pixel,
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    # 5. Filter rows with pixel_percentage >= 0.001
    filter(pixel_percentage >= 0.001)
  
  return(meanNDVI_PSIbin_species)
}


big_bin_df <- NDVI_PSIbin_multi(big_df)
write.csv(big_bin_df, "results_whole_Germany/NDVI_PSI_bin_summary_crossspecies.csv", row.names = FALSE)

cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)


big_bin_df$NDVI_species <- factor(big_bin_df$NDVI_species, levels = c("Oak", "Beech", "Spruce", "Pine"))
big_bin_df$PSI_species  <- factor(big_bin_df$PSI_species,  levels = c("Oak", "Beech", "Spruce", "Pine"))

big_bin_df$pixel_percent <- big_bin_df$pixel_percentage * 100

library(ggplot2)

ggplot(big_bin_df, aes(
  x = bin_median,
  y = avg_value,
  color = PSI_species,
  shape = PSI_species,
  size  = pixel_percent
)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = cb_palette, name = "PSI species") +
  scale_shape_manual(
    values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
    name = "species of soil water potential"
  ) +
  scale_size_continuous(
    name = "pixel percentage (%)",
    range = c(1, 8)
  ) +
  facet_wrap(~ NDVI_species, nrow = 2, ncol = 2) +
  labs(
    x = "soil water potential (kPa)",
    y = "NDVI quantiles"
  ) +
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
    legend.position = "right",
    legend.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 6)  # Make legend points larger for color/shape legend
    ),
    shape = guide_legend(
      override.aes = list(size = 6)  # Make legend points larger for color/shape legend
    ),
    size = guide_legend()            # Keeps the size legend as is
  )

# Save outputs
ggsave("results_whole_Germany/NDVI_vs_PSI_crossspecies_facet_by_NDVI.png", width = 12, height = 8, dpi = 300)

big_bin_df <- read.csv("results_whole_Germany/NDVI_PSI_bin_summary_crossspecies.csv")


library(dplyr)
library(broom)

fit_best_curve <- function(df) {
  # Fit linear
  fit_lm <- lm(avg_value ~ bin_median, data = df)
  aic_lm <- AIC(fit_lm)
  
  # Fit exponential
  fit_exp <- tryCatch({
    nls(avg_value ~ a + b * exp(c * bin_median),
        data = df,
        start = list(a = min(df$avg_value), 
                     b = (max(df$avg_value)-min(df$avg_value))/2, 
                     c = 0.001),
        control = nls.control(maxiter = 1000))
  }, error = function(e) NULL)
  
  # Model selection logic
  use_exp <- FALSE
  if (!is.null(fit_exp)) {
    aic_exp <- AIC(fit_exp)
    # significance check for 'c'
    tidy_exp <- broom::tidy(fit_exp)
    c_row <- tidy_exp[tidy_exp$term == "c", ]
    if (aic_exp < aic_lm && !is.null(c_row$p.value) && c_row$p.value < 0.05) {
      use_exp <- TRUE
    }
  }
  
  # Prediction grid (finer than original)
  pred_grid <- data.frame(bin_median = seq(min(df$bin_median), max(df$bin_median), length.out = 100))
  if (use_exp) {
    pred_grid$fit <- predict(fit_exp, newdata = pred_grid)
    pred_grid$model <- "exponential"
  } else {
    pred_grid$fit <- predict(fit_lm, newdata = pred_grid)
    pred_grid$model <- "linear"
  }
  pred_grid
}

library(tidyr)
library(purrr)

big_bin_df <- na.omit(big_bin_df)

fitted_df <- big_bin_df %>%
  group_by(NDVI_species, PSI_species) %>%
  group_modify(~ fit_best_curve(.x)) %>%
  ungroup()


cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

ggplot(big_bin_df, aes(
  x = bin_median,
  y = avg_value,
  color = PSI_species,
  shape = PSI_species,
  size  = pixel_percent
)) +
  geom_point(alpha = 0.85) +
  # Fitted lines:
  geom_line(data = fitted_df, 
            aes(x = bin_median, y = fit, color = PSI_species, linetype = model, group = interaction(PSI_species, model)),
            linewidth = 1, inherit.aes = FALSE) +
  scale_color_manual(values = cb_palette, name = "species of soil water potential") +
  scale_shape_manual(
    values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
    name = "species of soil water potential"
  ) +
  scale_size_continuous(
    name = "pixel percentage (%)",
    range = c(1, 8)
  ) +
  facet_wrap(~ NDVI_species, nrow = 2, ncol = 2) +
  labs(
    x = "soil water potential (kPa)",
    y = "NDVI quantiles"
  ) +
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
    legend.position = "right",
    legend.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 6)
    ),
    # shape = guide_legend(
    #   override.aes = list(size = 6)
    # ),
    size = guide_legend()
  )
ggsave("results_whole_Germany/NDVI_vs_PSI_crossspecies_facet_by_NDVI_fitted_line.png", width = 12, height = 8, dpi = 300)
