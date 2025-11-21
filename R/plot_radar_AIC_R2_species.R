# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggradar)
library(cowplot)
library(scales)

# ---- Global white background everywhere ----
theme_set(
  theme_minimal() +
    theme(
      plot.background        = element_rect(fill = "white", colour = NA),
      panel.background       = element_rect(fill = "white", colour = NA),
      legend.background      = element_rect(fill = "white", colour = NA),
      legend.box.background  = element_rect(fill = "white", colour = NA)
    )
)

# ---- Species order ----
species_levels <- c("Oak","Beech","Spruce","Pine")

# ---- Model palette (consistent across all panels) ----
model_palette <- c(
  linear      = "#E69F00",
  exponential = "#0072B2",
  poly2       = "#009E73",
  poly3       = "#F0E442"
)

# ---- Scaling helper ----
range01 <- function(x) {
  ux <- unique(x[is.finite(x)])
  if (length(ux) <= 1) return(ifelse(is.na(x), NA_real_, 0.5))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# =========================
# Data (from your LaTeX tables)
# =========================
panelA_long <- tibble::tribble(
  ~species, ~model,       ~AIC,   ~R2,
  "Oak",    "linear",      47.87, 0.97,
  "Oak",    "exponential", 39.20, 0.97,
  "Beech",  "linear",      66.81, 0.95,
  "Beech",  "exponential", NA,    NA,
  "Spruce", "linear",       2.74, 0.99,
  "Spruce", "exponential", -9.00, 0.99,
  "Pine",   "linear",      98.31, 0.73,
  "Pine",   "exponential",  8.97, 0.99
)

panelB_long <- tibble::tribble(
  ~species, ~model,       ~AIC,   ~R2,
  "Oak",    "linear",      29.92, 0.93,
  "Oak",    "exponential", 20.71, 0.97,
  "Beech",  "linear",      35.29, 0.85,
  "Beech",  "exponential", 22.14, 0.95,
  "Spruce", "linear",      37.00, 0.48,
  "Spruce", "exponential", 25.80, 0.84,
  "Pine",   "linear",      44.84, 0.76,
  "Pine",   "exponential", 19.46, 0.98
)

panelC_long <- tibble::tribble(
  ~species, ~model, ~AIC,   ~R2,
  "Oak",    "linear", 149.46, 0.96,
  "Oak",    "poly2",   99.34, 0.99,
  "Oak",    "poly3",   77.22, 0.99,
  "Beech",  "linear",  51.53, 0.99,
  "Beech",  "poly2",   53.47, 0.99,
  "Beech",  "poly3",   42.87, 1.00,
  "Spruce", "linear",  66.00, 0.98,
  "Spruce", "poly2",   50.27, 0.99,
  "Spruce", "poly3",    7.95, 1.00,
  "Pine",   "linear", 118.62, 0.85,
  "Pine",   "poly2",   94.60, 0.94,
  "Pine",   "poly3",   44.34, 0.99
)

# =========================
# Core builder: species = axes, lines = MODELS; metric = "AIC" or "R2"
# =========================
make_radar_models_by_metric <- function(
    long_tbl, model_levels, metric = c("AIC","R2"),
    panel_title_expr, 
    r_label = 1.1,
    r_label_offsets = c(Oak = 0.00, Beech = 0.12, Spruce = 0.00, Pine = 0.12),
    draw_lines = TRUE, point_size = 4.2, show_legend = FALSE
) {
  metric <- match.arg(metric)
  
  long_tbl <- long_tbl |>
    mutate(
      model   = factor(model, levels = model_levels),
      species = factor(species, levels = species_levels)
    ) |>
    arrange(species, model)
  
  # Select and scale per species
  scaled <- if (metric == "AIC") {
    long_tbl |>
      select(species, model, value = AIC) |>
      group_by(species) |>
      mutate(score = 1 - range01(value)) |>
      ungroup()
  } else {
    long_tbl |>
      select(species, model, value = R2) |>
      group_by(species) |>
      mutate(score = range01(value)) |>
      ungroup()
  }
  
  radar_df <- scaled |>
    mutate(model = as.character(model)) |>
    select(species, model, score) |>
    tidyr::pivot_wider(names_from = species, values_from = score) |>
    arrange(factor(model, levels = model_levels)) |>
    rename(group = model) |>
    mutate(across(-group, ~ ifelse(is.na(.), 0, .)))
  
  # ---- Axis labels (species) with per-species offsets + side-aware justification ----
  axis_labels <- species_levels
  n_vars <- length(axis_labels)
  angles <- seq(0, 2*pi, length.out = n_vars + 1)[-(n_vars + 1)]
  
  r_each <- r_label + as.numeric(r_label_offsets[axis_labels])
  axis_df <- data.frame(
    lab   = axis_labels,
    ang   = angles,
    r     = r_each
  ) |>
    mutate(
      x = r * sin(ang),
      y = r * cos(ang),
      hjust = dplyr::case_when(
        sin(ang) >  1e-8 ~ 0,   # right side
        sin(ang) < -1e-8 ~ 1,   # left side
        TRUE             ~ 0.5  # top/bottom
      ),
      vjust = dplyr::case_when(
        cos(ang) >  1e-8 ~ 1.1,  # top
        cos(ang) < -1e-8 ~ -0.1, # bottom
        TRUE             ~ 0.5
      )
    )
  
  # --- MORE breathing room so top/left labels never get cropped
  expand_r <- max(r_each) + 0.35   # was +0.18
  
  # Colors for present models, in model_levels order
  cols <- unname(model_palette[as.character(radar_df$group)])
  
  p <- ggradar(
    radar_df,
    axis.labels               = rep("", n_vars),  # we draw custom labels
    group.colours             = cols,
    grid.min = 0, grid.mid = 0.5, grid.max = 1,
    background.circle.colour  = "white",
    gridline.min.colour       = "grey80",
    gridline.mid.colour       = "grey70",
    gridline.max.colour       = "grey60",
    axis.label.size           = 0,
    group.line.width          = if (draw_lines) 1.2 else 0,  # points only by default
    group.point.size          = point_size
  ) +
    labs(title = panel_title_expr) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = if (show_legend) "top" else "none",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(22, 30, 32, 30)   # a touch more outside margin
    ) +
    coord_equal() +
    coord_cartesian(clip = "off") +
    expand_limits(x = c(-expand_r, expand_r), y = c(-expand_r, expand_r)) +
    geom_text(
      data = axis_df,
      aes(x, y, label = lab, hjust = hjust, vjust = vjust),
      parse = FALSE,
      size = 4.8,
      lineheight = 0.95
    ) 
  
  if (show_legend) {
    p <- p + guides(colour = guide_legend(
      title = NULL, nrow = 1, byrow = TRUE,
      override.aes = list(linetype = 0, size = point_size)  # point-only legend keys
    ))
  }
  p
  
  
}

# ---- helper: add under-panel tag ----
add_under_tag <- function(plot, tag_text) {
  tag <- ggdraw() + draw_label(tag_text, x = 0.5, y = 0.95, hjust = 0.5, vjust = 0)
  plot_grid(plot, tag, ncol = 1, rel_heights = c(1, 0.06))
}

# =========================
# Build the six panels
# =========================
title_A <- expression(NDVI~quantiles~" "~"~"~Psi[" "~soil])
title_B <- expression(NDVI~quantiles~" "~"~"~T[d])
title_C <- expression(T[d]~" "~"~"~Psi[" "~soil])

# Increased base radius and offsets to ensure all labels are visible
base_r <- 1.1   

# More generous offsets for all labels, especially Oak (top) and Spruce (left)
offs <- c(Oak = 0.23, Beech = 0.15, Spruce = 0.23, Pine = 0.15)

# Row 1 (AIC)
pA_AIC <- make_radar_models_by_metric(panelA_long, c("linear","exponential"), "AIC",
                                      title_A, r_label = base_r, r_label_offsets = offs)
pB_AIC <- make_radar_models_by_metric(panelB_long, c("linear","exponential"), "AIC",
                                      title_B, r_label = base_r, r_label_offsets = offs)
pC_AIC <- make_radar_models_by_metric(panelC_long, c("linear","poly2","poly3"), "AIC",
                                      title_C, r_label = base_r, r_label_offsets = offs)

# Row 2 (R²)
pA_R2  <- make_radar_models_by_metric(panelA_long, c("linear","exponential"), "R2",
                                      title_A, r_label = base_r, r_label_offsets = offs)
pB_R2  <- make_radar_models_by_metric(panelB_long, c("linear","exponential"), "R2",
                                      title_B, r_label = base_r, r_label_offsets = offs)
pC_R2  <- make_radar_models_by_metric(panelC_long, c("linear","poly2","poly3"), "R2",
                                      title_C, r_label = base_r, r_label_offsets = offs)

# ---- Tags
pA_AIC_tag <- add_under_tag(pA_AIC, "(a)")
pB_AIC_tag <- add_under_tag(pB_AIC, "(b)")
pC_AIC_tag <- add_under_tag(pC_AIC, "(c)")
pA_R2_tag  <- add_under_tag(pA_R2,  "(d)")
pB_R2_tag  <- add_under_tag(pB_R2,  "(e)")
pC_R2_tag  <- add_under_tag(pC_R2,  "(f)")

# =========================
# Shared legend (all models), point-only keys
# =========================
dummy_df <- data.frame(
  group = c("linear","exponential","poly2","poly3"),
  Oak   = runif(4, 0.4, 0.7),
  Beech = runif(4, 0.4, 0.7),
  Spruce= runif(4, 0.4, 0.7),
  Pine  = runif(4, 0.4, 0.7)
)
legend_plot <- ggradar(
  dummy_df,
  axis.labels   = rep("", 4),
  group.colours = unname(model_palette[dummy_df$group]),
  grid.min = 0, grid.mid = 0.5, grid.max = 1,
  background.circle.colour  = "white",
  axis.label.size = 0,
  group.line.width = 0,
  group.point.size = 4.2
) + guides(colour = guide_legend(
  title = NULL, nrow = 1, byrow = TRUE,
  override.aes = list(linetype = 0, size = 4.2)
))
legend_shared <- cowplot::get_legend(legend_plot)

# =========================
# Arrange 2 rows x 3 columns, legend on top, with left-side row labels
# =========================
row1 <- plot_grid(pA_AIC_tag, pB_AIC_tag, pC_AIC_tag, nrow = 1, rel_widths = c(1,1,1))
row2 <- plot_grid(pA_R2_tag,  pB_R2_tag,  pC_R2_tag,  nrow = 1, rel_widths = c(1,1,1))

label_theme <- theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )
label_AIC <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "AIC", angle = 0, size = 8) + label_theme
label_R2  <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "R^2", parse = TRUE, angle = 0, size = 8) + label_theme

row1_labeled <- plot_grid(label_AIC, row1, ncol = 2, rel_widths = c(0.06, 1))
row2_labeled <- plot_grid(label_R2,  row2, ncol = 2, rel_widths = c(0.06, 1))

six_panel <- plot_grid(
  legend_shared,
  row1_labeled,
  row2_labeled,
  ncol = 1,
  rel_heights = c(0.12, 1, 1)
) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

print(six_panel)

# ---- Save ----
out_dir <- "results_rootzone/Figures/supplementary"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(
  filename = file.path(out_dir, "radar_6panels_models_by_metric.png"),
  plot = six_panel,
  width = 16, height = 10, dpi = 300, bg = "white"
)