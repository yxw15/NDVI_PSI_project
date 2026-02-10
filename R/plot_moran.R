setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(tidyverse)
library(ggplot2)

# ----------------------------
# Theme + colors
# ----------------------------
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

sp_cols <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")

# ----------------------------
# Read all csv in data/
# ----------------------------
files <- list.files("results_moran/data", pattern = "\\.csv$", full.names = TRUE)

combined <- purrr::map_dfr(files, function(f) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  # infer "Panel" from filename
  panel <- dplyr::case_when(
    str_detect(basename(f), "Quantiles") ~ "NDVI quantiles (rank)",
    str_detect(basename(f), "soil_water_potential") ~ "soil water potential (kPa)",
    str_detect(basename(f), "transpiration_deficit") ~ "transpiration deficit (mm)",
    TRUE ~ "Other"
  )
  
  df %>% mutate(Panel = panel)
})

# ----------------------------
# Clean / order
# ----------------------------
combined <- combined %>%
  mutate(
    Distance = as.numeric(Distance),
    Moran_I  = as.numeric(Moran_I),
    Species  = factor(Species, levels = c("Oak", "Beech", "Spruce", "Pine")),
    Panel    = factor(Panel, levels = c("NDVI quantiles (rank)", "soil water potential (kPa)", "transpiration deficit (mm)")),
    Significance = if_else(is.na(Significance), "", as.character(Significance))
  )

# ----------------------------
# Plot: 3 panels + species colors + "*" below points
# ----------------------------
p <- ggplot(combined, aes(x = Distance, y = Moran_I, color = Species, group = Species)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_text(aes(label = Significance), vjust = 2.2, show.legend = FALSE) +
  facet_wrap(~Panel, nrow = 1, scales = "fixed") +
  scale_color_manual(values = sp_cols, drop = FALSE) +
  labs(
    title = "",
    x = "distance (m)",
    y = "Moran's I",
    color = ""
  ) +
  base_theme

p
ggsave("results_rootzone/Figures_till2022/SI_PSI/SI_moran_distance.png", plot = p, width = 12, height = 6, dpi = 300)
