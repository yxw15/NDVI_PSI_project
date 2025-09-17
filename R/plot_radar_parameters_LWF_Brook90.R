setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# ---- Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggradar)

# ---- Color palette ----
cb_palette <- c(
  Oak    = "#E69F00",
  Beech  = "#0072B2",
  Spruce = "#009E73",
  Pine   = "#F0E442"
)

# ---- Parameters table (includes mean_vpd) ----
params_tbl <- tibble::tribble(
  ~species, ~LAI, ~gmax,  ~G,      ~Psi_CR, ~root_depth, ~beta_root, ~mean_vpd,
  "Oak",     4.5, 0.0055, 0.02475, -2.5,    -1.8,        0.976,      0.605,
  "Beech",   6.0, 0.0042, 0.02520, -2.0,    -1.4,        0.976,      0.512,
  "Spruce",  5.5, 0.0035, 0.01925, -2.0,    -1.2,        0.966,      0.455,
  "Pine",    3.5, 0.0055, 0.01925, -2.5,    -1.8,        0.966,      0.600
)

# ---- Prepare magnitudes for signed variables ----
params_aug <- params_tbl %>%
  mutate(
    Psi_CR_abs   = abs(Psi_CR),
    root_depth_m = abs(root_depth)
  ) %>%
  select(species, LAI, gmax, G, Psi_CR_abs, root_depth_m, beta_root, mean_vpd)

# ---- Min–max scale to 0–1 across species (per variable) ----
range01 <- function(x) {
  if (length(unique(x[is.finite(x)])) <= 1) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

params_scaled <- params_aug %>%
  mutate(across(where(is.numeric), range01))

# ggradar expects first column to be the group label
radar_df <- params_scaled %>%
  rename(group = species)

# ---- 1) Give ggradar simple (character) labels to avoid expression->data.frame error ----
# (Order must match columns after 'group')
axis_labels_plain <- c(
  "LAI",
  "g[max]",
  "G",
  "|Psi[CR]|",
  "|root[depth]|",
  "beta[root]",
  "bar(VPD)"
)

# ---- Radar base (hide built-in labels so we can add parsed ones) ----
radar_base <- ggradar(
  radar_df,
  axis.labels = axis_labels_plain,  # character vector OK
  group.colours = cb_palette,
  grid.min = 0, grid.mid = 0.5, grid.max = 1,
  background.circle.colour = "white",
  gridline.mid.colour = "grey80",
  gridline.min.colour = "grey90",
  gridline.max.colour = "grey90",
  axis.label.size = 0,              # hide the default labels
  group.line.width = 1.2,
  group.point.size = 3
) +
  theme(legend.position = "top")

# ---- Labels further outside the circle ----
n_vars <- ncol(radar_df) - 1
angles <- seq(0, 2*pi, length.out = n_vars + 1)[- (n_vars + 1)]

r_label <- 1.35   # increase to move labels further away (try 1.35–1.5)

axis_df <- data.frame(
  x = r_label * sin(angles),
  y = r_label * cos(angles),
  lab <- c(
    "LAI",
    "italic(g)[max]",
    "G",
    '"|"*italic(Psi)*phantom(0)[CR]*"|"',   # |Ψ_CR| with extra space
    '"|"*root~depth*"|"',                   # |root depth|
    "italic(beta)*phantom(0)[root]",        # β_root with extra space
    "bar(VPD)"
  )
)

radar_plot <- radar_base +
  geom_text(
    data = axis_df,
    aes(x, y, label = lab),
    parse = TRUE,
    size = 5
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    legend.position = "top",
  )

print(radar_plot)

out_dir <- "results_rootzone/Figures/supplementary"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(
  filename = file.path(out_dir, "radar_parameters.png"),
  plot = radar_plot,
  width = 7.5, height = 7.5, dpi = 300
)