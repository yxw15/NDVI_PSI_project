# ─── 0. Setup ────────────────────────────────────────────────────────────────
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Libraries
library(terra)
library(dplyr)
library(tidyverse)

# Parameters
years   <- 2003:2022
months  <- c("July", "August")
species <- c("Oak", "Beech", "Spruce", "Pine")
base_dir <- "results_monthly_rootzone"

# ─── 1. Read all PSI rasters into one dataframe ──────────────────────────────
read_psi <- function(year, month, species_name) {
  file_path <- file.path(base_dir, month, species_name, paste0("psi_", year, ".tif"))
  r <- rast(file_path)
  
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "soil_water_potential"
  
  df %>%
    mutate(
      species = species_name,
      month   = month,
      year    = year
    )
}

psi_all <- expand.grid(year = years, month = months, species = species, stringsAsFactors = FALSE) %>%
  pmap_dfr(~ read_psi(..1, ..2, ..3))


# ─── 1. Read all PSI rasters into one dataframe ──────────────────────────────
read_tdiff <- function(year, month, species_name) {
  file_path <- file.path(base_dir, month, species_name, paste0("tdiff_", year, ".tif"))
  r <- rast(file_path)
  
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "transpiration_deficit"
  
  df %>%
    mutate(
      species = species_name,
      month   = month,
      year    = year
    )
}

tdiff_all <- expand.grid(year = years, month = months, species = species, stringsAsFactors = FALSE) %>%
  pmap_dfr(~ read_tdiff(..1, ..2, ..3))

combined_df <- psi_all %>%
  left_join(tdiff_all, by = c("x", "y", "species", "month", "year"))

combined_df$area <- "Germany"
head(combined_df)
save(combined_df, file = "results_rootzone/Data/PSI_TDiff_species_month_year_Germany.RData")

load("results/Data/AllSpecies_AllMonths_rootzone.RData")
data <- combined %>% filter(month %in% c("July", "August"))
data <- data %>% filter(Quantiles > 0 & year <=2022)
data <- na.omit(data)
data$area <- "species"

data_fixed <- data %>%
  select(any_of(names(combined_df))) %>%      # keep only shared columns
  select(all_of(names(combined_df)))  

data_fixed <- data_fixed %>%
  mutate(
    year = as.integer(year),
    species = as.character(species),
    month = as.character(month),
    area = as.character(area)
  )

combined_df2 <- combined_df %>%
  mutate(
    year = as.integer(year),
    species = as.character(species),
    month = as.character(month),
    area = as.character(area)
  )

combined_all <- bind_rows(combined_df2, data_fixed)
save(combined_all, file = "results_rootzone/Data/PSI_TDiff_species_month_year_area.RData")

##################### Germany and species-located data #################
load("results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")

library(tidyverse)

TDiff_PSIbin_area <- function(df, bin_width = 50, min_count = 1000) {
  library(dplyr)
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df_binned <- df %>%
    mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    filter(!is.na(PSI_bin))  # <-- important
  
  totals <- df_binned %>%
    group_by(area, species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  df_binned %>%
    group_by(area, species, PSI_bin) %>%
    summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      avg_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE),
      count = n(),
      .groups = "drop"
    ) %>%
    left_join(totals, by = c("area", "species")) %>%
    mutate(percentage = count / total_pixels) %>%
    filter(count >= min_count) %>%
    tidyr::drop_na()
}

### one panel
plot_TDiff_PSIbin_byArea <- function(data, figure_output = NULL,
                                     bin_width = 50, min_count = 1000) {
  
  library(tidyverse)
  library(scales)
  
  df <- TDiff_PSIbin_area(data,
                          bin_width = bin_width,
                          min_count = min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      area    = factor(area, levels = c("Germany", "species"))
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  xcol <- "bin_mean"
  ycol <- "avg_transpiration_deficit"
  
  degrees <- c(1, 2, 3)
  
  fits <- df %>%
    group_by(area, species) %>%
    group_modify(~{
      sdf <- .x %>% select(all_of(c(xcol, ycol)))
      if (nrow(sdf) < 5) return(tibble(best_fit = list(NULL)))
      
      models <- map(degrees, \(d) {
        tryCatch(
          lm(reformulate(paste0("poly(", xcol, ", ", d, ")"), ycol), data = sdf),
          error = \(e) NULL
        )
      })
      
      aics <- map_dbl(models, \(m) if (is.null(m)) Inf else AIC(m))
      tibble(best_fit = list(models[[which.min(aics)]]))
    }) %>%
    ungroup()
  
  pred_df <- fits %>%
    left_join(
      df %>%
        group_by(area, species) %>%
        summarise(
          xmin = min(.data[[xcol]], na.rm = TRUE),
          xmax = max(.data[[xcol]], na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("area", "species")
    ) %>%
    mutate(
      xseq = map2(xmin, xmax, \(a, b) seq(a, b, length.out = 200)),
      pred_tbl = pmap(list(best_fit, xseq, area, species), \(fit, xs, ar, sp) {
        if (is.null(fit)) return(tibble())
        nd <- tibble(bin_mean = xs, area = ar, species = sp)
        nd %>% mutate(predicted = predict(fit, newdata = nd))
      })
    ) %>%
    select(pred_tbl) %>%
    unnest(pred_tbl)
  
  plot_clean <- ggplot() +
    geom_point(
      data = df,
      aes(x = .data[[xcol]],
          y = .data[[ycol]],
          color = species,
          shape = species,
          size = percentage),
      alpha = 0.7
    ) +
    geom_line(
      data = pred_df,
      aes(x = bin_mean,
          y = predicted,
          color = species,
          linetype = area),
      linewidth = 1
    ) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_linetype_manual(values = c("Germany" = "solid", "species" = "dashed"),
                          name = "") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color    = guide_legend(order = 1),
      linetype = guide_legend(order = 2),
      size     = guide_legend(order = 3)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
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
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = plot_clean,
           width = 8, height = 6, dpi = 300)
  }
  
  invisible(plot_clean)
}

### two panels
plot_TDiff_PSIbin_byArea <- function(data, figure_output = NULL,
                                     bin_width = 50, min_count = 1000) {
  
  library(tidyverse)
  library(scales)
  
  df <- TDiff_PSIbin_area(data,
                          bin_width = bin_width,
                          min_count = min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      area    = factor(area, levels = c("Germany", "species"))
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  xcol <- "bin_mean"
  ycol <- "avg_transpiration_deficit"
  
  degrees <- c(1, 2, 3)
  
  # -------------------------
  # Model fitting (per area × species)
  # -------------------------
  fits <- df %>%
    group_by(area, species) %>%
    group_modify(~{
      sdf <- .x %>% select(all_of(c(xcol, ycol)))
      if (nrow(sdf) < 5) return(tibble(best_fit = list(NULL)))
      
      models <- map(degrees, \(d) {
        tryCatch(
          lm(reformulate(paste0("poly(", xcol, ", ", d, ")"), ycol), data = sdf),
          error = \(e) NULL
        )
      })
      
      aics <- map_dbl(models, \(m) if (is.null(m)) Inf else AIC(m))
      tibble(best_fit = list(models[[which.min(aics)]]))
    }) %>%
    ungroup()
  
  # -------------------------
  # Predictions
  # -------------------------
  pred_df <- fits %>%
    left_join(
      df %>%
        group_by(area, species) %>%
        summarise(
          xmin = min(.data[[xcol]], na.rm = TRUE),
          xmax = max(.data[[xcol]], na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("area", "species")
    ) %>%
    mutate(
      xseq = map2(xmin, xmax, \(a, b) seq(a, b, length.out = 200)),
      pred_tbl = pmap(list(best_fit, xseq, area, species), \(fit, xs, ar, sp) {
        if (is.null(fit)) return(tibble())
        nd <- tibble(bin_mean = xs, area = ar, species = sp)
        nd %>% mutate(predicted = predict(fit, newdata = nd))
      })
    ) %>%
    select(pred_tbl) %>%
    unnest(pred_tbl)
  
  # -------------------------
  # Plot: TWO PANELS by area
  # -------------------------
  plot_clean <- ggplot() +
    geom_point(
      data = df,
      aes(x = .data[[xcol]],
          y = .data[[ycol]],
          color = species,
          shape = species,
          size = percentage),
      alpha = 0.7
    ) +
    geom_line(
      data = pred_df,
      aes(x = bin_mean,
          y = predicted,
          color = species),
      linewidth = 1
    ) +
    facet_wrap(~ area, ncol = 2) +   # ← THIS IS THE KEY LINE
    scale_color_manual(values = cb_palette, name = "") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17,
                                  "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
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
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = plot_clean,
           width = 10, height = 5, dpi = 300)
  }
  
  invisible(plot_clean)
}

### two panels with confidenct bands
plot_TDiff_PSIbin_byArea <- function(data, figure_output = NULL,
                                     bin_width = 50, min_count = 1000,
                                     conf_level = 0.95) {
  
  library(tidyverse)
  library(scales)
  
  df <- TDiff_PSIbin_area(data,
                          bin_width = bin_width,
                          min_count = min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      area    = factor(area, levels = c("Germany", "species"))
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  xcol <- "bin_mean"
  ycol <- "avg_transpiration_deficit"
  
  degrees <- c(1, 2, 3)
  
  # -------------------------
  # Model fitting (per area × species)
  # -------------------------
  fits <- df %>%
    group_by(area, species) %>%
    group_modify(~{
      sdf <- .x %>% select(all_of(c(xcol, ycol)))
      if (nrow(sdf) < 5) return(tibble(best_fit = list(NULL)))
      
      models <- map(degrees, \(d) {
        tryCatch(
          lm(reformulate(paste0("poly(", xcol, ", ", d, ")"), ycol), data = sdf),
          error = \(e) NULL
        )
      })
      
      aics <- map_dbl(models, \(m) if (is.null(m)) Inf else AIC(m))
      tibble(best_fit = list(models[[which.min(aics)]]))
    }) %>%
    ungroup()
  
  # -------------------------
  # Predictions + confidence bands
  # -------------------------
  alpha <- 1 - conf_level
  
  pred_df <- fits %>%
    left_join(
      df %>%
        group_by(area, species) %>%
        summarise(
          xmin = min(.data[[xcol]], na.rm = TRUE),
          xmax = max(.data[[xcol]], na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("area", "species")
    ) %>%
    mutate(
      xseq = map2(xmin, xmax, \(a, b) seq(a, b, length.out = 200)),
      pred_tbl = pmap(list(best_fit, xseq, area, species), \(fit, xs, ar, sp) {
        if (is.null(fit)) return(tibble())
        
        nd <- tibble(bin_mean = xs, area = ar, species = sp)
        
        pr <- predict(fit, newdata = nd, se.fit = TRUE)
        df_res <- df.residual(fit)
        tcrit <- qt(1 - alpha/2, df = df_res)
        
        nd %>%
          mutate(
            predicted = as.numeric(pr$fit),
            se = as.numeric(pr$se.fit),
            lower = predicted - tcrit * se,
            upper = predicted + tcrit * se
          )
      })
    ) %>%
    select(pred_tbl) %>%
    unnest(pred_tbl)
  
  # -------------------------
  # Plot: TWO PANELS by area + confidence bands
  # -------------------------
  plot_clean <- ggplot() +
    
    # confidence ribbon (behind lines/points)
    geom_ribbon(
      data = pred_df,
      aes(x = bin_mean, ymin = lower, ymax = upper, fill = species, group = interaction(area, species)),
      alpha = 0.18,
      color = NA
    ) +
    
    geom_point(
      data = df,
      aes(x = .data[[xcol]],
          y = .data[[ycol]],
          color = species,
          shape = species,
          size = percentage),
      alpha = 0.7
    ) +
    
    geom_line(
      data = pred_df,
      aes(x = bin_mean,
          y = predicted,
          color = species),
      linewidth = 1
    ) +
    facet_wrap(
      ~ area,
      ncol = 2,
      labeller = as_labeller(c(
        Germany = "(a) Germany",
        species = "(b) species"
      ))
    ) + 
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +   # ribbon uses fill
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17,
                                  "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit (mm)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
  
  
  print(plot_clean)
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = plot_clean,
           width = 10, height = 5, dpi = 300)
  }
  
  invisible(plot_clean)
}

# plot_TDiff_PSIbin_byArea(combined_all, "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_area_till2022.png")
plot_TDiff_PSIbin_byArea(DT_combined, "results_rootzone/test/TDiff_PSIbin_area_till2022.png")


### two panels GAM (k=3)
plot_TDiff_PSIbin_byArea <- function(data, figure_output = NULL,
                                     bin_width = 50, min_count = 1000,
                                     conf_level = 0.95) {
  
  library(tidyverse)
  library(scales)
  library(mgcv)
  
  df <- TDiff_PSIbin_area(data,
                          bin_width = bin_width,
                          min_count = min_count) %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      area    = factor(area, levels = c("Germany", "species"))
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  xcol <- "bin_mean"
  ycol <- "avg_transpiration_deficit"
  
  alpha <- 1 - conf_level
  
  # -------------------------
  # Fit GAM + predict per area × species
  # -------------------------
  pred_df <- df %>%
    group_by(area, species) %>%
    group_modify(~{
      sdf <- .x %>% select(all_of(c(xcol, ycol)))
      
      if (nrow(sdf) < 5) return(tibble())
      
      # Fit GAM with k = 3
      fit <- mgcv::gam(
        reformulate(termlabels = "s(bin_mean, k = 3)", response = ycol),
        data = sdf,
        method = "REML"
      )
      
      xr <- range(sdf[[xcol]], na.rm = TRUE)
      xseq <- seq(xr[1], xr[2], length.out = 200)
      nd <- tibble(bin_mean = xseq)
      
      pr <- predict(fit, newdata = nd, se.fit = TRUE)
      
      # For mgcv GAMs, using normal approx is standard
      zcrit <- qnorm(1 - alpha/2)
      
      nd %>%
        mutate(
          predicted = as.numeric(pr$fit),
          se = as.numeric(pr$se.fit),
          lower = predicted - zcrit * se,
          upper = predicted + zcrit * se
        )
    }) %>%
    ungroup()
  
  # -------------------------
  # Plot: TWO PANELS by area + confidence bands
  # -------------------------
  plot_clean <- ggplot() +
    
    geom_ribbon(
      data = pred_df,
      aes(x = bin_mean, ymin = lower, ymax = upper,
          fill = species, group = interaction(area, species)),
      alpha = 0.18,
      color = NA
    ) +
    
    geom_point(
      data = df,
      aes(x = .data[[xcol]],
          y = .data[[ycol]],
          color = species,
          shape = species,
          size = percentage),
      alpha = 0.7
    ) +
    
    geom_line(
      data = pred_df,
      aes(x = bin_mean,
          y = predicted,
          color = species),
      linewidth = 1
    ) +
    
    facet_wrap(
      ~ area,
      ncol = 2,
      labeller = as_labeller(c(
        Germany = "(a) Germany",
        species = "(b) species"
      ))
    ) +
    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17,
                                  "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit (mm)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
  
  print(plot_clean)
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = plot_clean,
           width = 10, height = 5, dpi = 300)
  }
  
  invisible(plot_clean)
}

# plot_TDiff_PSIbin_byArea(combined_all, "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_area_gam3_till2022.png")
plot_TDiff_PSIbin_byArea(DT_combined, "results_rootzone/test/TDiff_PSIbin_area_gam3_till2022.png")


### different
TDiff_diff_species_minus_Germany <- function(df, bin_width = 50, min_count = 1000) {
  library(tidyverse)
  
  # define common bin breaks (global so Germany/species use identical bins)
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df_binned <- df %>%
    mutate(
      PSI_bin = cut(
        soil_water_potential,
        breaks = bin_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # mean TDiff for species-mask (area == "species")
  sp_tbl <- df_binned %>%
    filter(area == "species") %>%
    group_by(species, PSI_bin) %>%
    summarise(
      bin_mean = mean(soil_water_potential, na.rm = TRUE),
      sp_mean  = mean(transpiration_deficit, na.rm = TRUE),
      sp_n     = n(),
      .groups  = "drop"
    ) %>%
    filter(sp_n >= min_count)
  
  # mean TDiff for Germany baseline (area == "Germany")
  de_tbl <- df_binned %>%
    filter(area == "Germany") %>%
    group_by(species, PSI_bin) %>%
    summarise(
      de_mean = mean(transpiration_deficit, na.rm = TRUE),
      de_n    = n(),
      .groups = "drop"
    ) %>%
    filter(de_n >= min_count)
  
  # join + difference
  diff_df <- sp_tbl %>%
    inner_join(de_tbl, by = c("species", "PSI_bin")) %>%
    mutate(diff_transpiration_deficit = sp_mean - de_mean) %>%
    select(species, PSI_bin, bin_mean, sp_mean, de_mean, sp_n, de_n, diff_transpiration_deficit)
  
  diff_df
}

plot_TDiff_diff <- function(data, figure_output = NULL,
                            bin_width = 50, min_count = 1000) {
  
  library(tidyverse)
  
  df_diff <- TDiff_diff_species_minus_Germany(data, bin_width, min_count) %>%
    mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  p <- ggplot(df_diff,
              aes(x = bin_mean,
                  y = diff_transpiration_deficit,
                  color = species)) +
    
    geom_hline(yintercept = 0,
               linewidth = 0.6,
               linetype = "dashed",
               color = "grey40") +
    
    geom_point(alpha = 0.5, size = 2) +
    
    # geom_smooth(
    #   aes(fill = species),          # ← key line
    #   se = TRUE,
    #   method = "loess",
    #   span = 0.7,
    #   linewidth = 1,
    #   alpha = 0.25
    # ) +
    
    geom_smooth(
      aes(fill = species),
      method = "gam",
      formula = y ~ s(x, k = 3),
      se = TRUE
    ) +

    
    scale_color_manual(values = cb_palette, name = "") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    
    labs(
      x = "soil water potential (kPa)",
      y = expression(Delta[T[d](species - Germany)]~"(mm)")
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = p, width = 8, height = 6, dpi = 300)
  }
  
  invisible(p)
}

# plot_TDiff_diff(combined_all, "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_area_diff_till2022.png")
plot_TDiff_diff(DT_combined, "results_rootzone/test/TDiff_PSIbin_area_diff_till2022.png")

### different bar plot 
plot_TDiff_diff_bar <- function(data, figure_output = NULL,
                                bin_width = 50, min_count = 1000) {
  
  library(tidyverse)
  
  df_bar <- TDiff_diff_species_minus_Germany(data, bin_width, min_count) %>%
    mutate(
      species = factor(species,
                       levels = c("Oak", "Beech", "Spruce", "Pine"))
    ) %>%
    group_by(species) %>%
    summarise(
      mean_diff = mean(diff_transpiration_deficit, na.rm = TRUE),
      se_diff   = sd(diff_transpiration_deficit, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  p <- ggplot(df_bar,
              aes(x = species,
                  y = mean_diff,
                  fill = species)) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed",
               linewidth = 0.6,
               color = "grey40") +
    
    geom_col(width = 0.7, alpha = 0.9) +
    
    geom_errorbar(
      aes(ymin = mean_diff - se_diff,
          ymax = mean_diff + se_diff),
      width = 0.2,
      linewidth = 0.6
    ) +
    
    scale_fill_manual(values = cb_palette, guide = "none") +
    
    labs(
      x = "",
      y = expression(mean~Delta[T[d](species - Germany)]~"(mm)")
    ) +
    
    theme_minimal() +
    theme(
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
    ggsave(figure_output, plot = p, width = 6, height = 5, dpi = 300)
  }
  
  invisible(p)
}

plot_TDiff_diff_bar(DT_combined, "results_rootzone/test/TDiff_PSIbin_area_diff_bar_till2022.png")

### Germany and species-located for species panels ###
plot_TDiff_PSIbin_species4_area2 <- function(
    df,
    out_file,
    bin_width = 50,
    min_count = 1000,
    k = 6,
    bs = "cs",
    width = 12,
    height = 8,
    dpi = 300,
    reverse_x = FALSE   # optional: TRUE to reverse PSI axis
) {
  # species order + colors for the two "area" curves
  sp_levels <- c("Oak", "Beech", "Spruce", "Pine")
  area_cols <- c("Germany" = "blue", "species" = "red")
  
  p <- TDiff_PSIbin_area(df, bin_width = bin_width, min_count = min_count) %>%
    dplyr::mutate(
      area    = factor(area, levels = c("Germany", "species")),
      species = factor(species, levels = sp_levels)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = bin_mean, y = avg_transpiration_deficit)) +
    
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
      x = "soil water potential (PSI)",
      y = "transpiration deficit (mm)"
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
  
  # optional PSI reverse axis
  if (isTRUE(reverse_x)) {
    p <- p + ggplot2::scale_x_reverse()
  }
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  
  print(p)
  invisible(p)
}

plot_TDiff_PSIbin_species4_area2(
  df = DT_combined,
  out_file = "results_rootzone/test/TDiff_vs_PSI_species4_area2.png",
  bin_width = 50,
  min_count = 1000,
  reverse_x = FALSE
)

### poly 3 ###
plot_TDiff_PSIbin_poly3 <- function(
    df,
    out_file = NULL,
    bin_width = 50,
    min_count = 1000,
    width = 12,
    height = 8,
    dpi = 300,
    reverse_x = FALSE
) {
  library(tidyverse)
  library(scales)
  
  # 1. Prepare Data
  sp_levels <- c("Oak", "Beech", "Spruce", "Pine")
  area_cols <- c("Germany" = "blue", "species" = "red") # Blue and Red-ish
  
  plot_df <- TDiff_PSIbin_area(df, bin_width = bin_width, min_count = min_count) %>%
    mutate(
      area    = factor(area, levels = c("Germany", "species")),
      species = factor(species, levels = sp_levels)
    )
  
  # 2. Build Plot
  p <- ggplot(plot_df, aes(x = bin_mean, y = avg_transpiration_deficit, 
                           color = area, fill = area, group = area)) +
    
    # Points sized by percentage
    geom_point(
      aes(shape = area, size = percentage),
      alpha = 0.6
    ) +
    
    # 3rd Degree Polynomial Smooth with Confidence Band
    # formula = y ~ poly(x, 3) provides the cubic fit you requested
    geom_smooth(
      method = "lm",
      formula = y ~ poly(x, 3),
      se = TRUE,      # Enables the confidence band
      alpha = 0.2,    # Transparency of the ribbon
      linewidth = 1.2
    ) +
    
    # 4 Panels by species
    facet_wrap(~ species, nrow = 2) +
    
    # Scales and Aesthetics
    scale_color_manual(values = area_cols, name = "Area Scope") +
    scale_fill_manual(values = area_cols, guide = "none") +
    scale_shape_manual(values = c("Germany" = 16, "species" = 17), name = "Area Scope") +
    scale_size_continuous(
      name = "pixels per bin (%)",
      range = c(1, 8),
      labels = percent_format(accuracy = 1)
    ) +
    
    labs(
      x = "soil water potential (kPa)",
      y = "transpiration deficit (mm)"
    ) +
    
    # Clean Theme
    theme_minimal() +
    ggplot2::theme_minimal() +
    theme(
      axis.text    = element_text(size = 14),
      axis.title   = element_text(size = 16, face = "bold"),
      strip.text   = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
    # ggplot2::theme(
    #   plot.title = ggplot2::element_blank(),
    #   axis.text.x = ggplot2::element_text(size = 14),
    #   axis.text.y = ggplot2::element_text(size = 14),
    #   axis.title  = ggplot2::element_text(size = 16, face = "bold"),
    #   legend.position = "bottom",
    #   legend.text = ggplot2::element_text(size = 14),
    #   legend.title = ggplot2::element_text(size = 14, face = "bold"),
    #   strip.text = ggplot2::element_text(size = 16, face = "bold"),
    #   panel.background = ggplot2::element_rect(fill = "white"),
    #   plot.background  = ggplot2::element_rect(fill = "white", color = "white"),
    #   panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.4),
    #   panel.grid.minor = ggplot2::element_blank()
    # )
  
  # Optional: Reverse X-axis (common for PSI/water potential)
  if (isTRUE(reverse_x)) {
    p <- p + scale_x_reverse()
  }
  
  # 3. Save and Return
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
  }
  
  return(p)
}

plot_TDiff_PSIbin_poly3(
  df = DT_combined,
  out_file = "results_rootzone/test/TDiff_vs_PSI_species4_area2_poly3.png",
  bin_width = 50,
  min_count = 1000,
  reverse_x = TRUE
)
