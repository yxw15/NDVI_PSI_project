# ─── 0. Setup ────────────────────────────────────────────────────────────────
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0001/yixuan/NDVI_PSI_project")

# Libraries
library(terra)
library(tidyverse)

# Parameters
years   <- 2003:2022
months  <- c("July", "August")
species <- c("Oak", "Beech", "Spruce", "Pine")
base_dir <- "results_monthly_rootzone"
bin_width <- 3  # PSI bin width in kPa

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

# ─── 2. Reshape to wide format ───────────────────────────────────────────────
tdiff_wide <- tdiff_all %>%
  pivot_wider(
    names_from  = species,
    values_from = transpiration_deficit
  )

head(tdiff_wide)

# ─── 3. Build pairwise tdiff comparisons (x-species vs y-species) ─────────────
tdiff_long <- tdiff_wide %>%
  pivot_longer(
    cols = all_of(species),
    names_to = "species",
    values_to = "tdiff_value"
  )

head(tdiff_long)

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

head(tdiff_pairs)

# ─── 4. Bin x-species tdiff by 50 kPa using cut() ─────────────────────────────
tdiff_summary <- tdiff_pairs %>%
  group_by(species_x) %>%              # do this separately for each panel
  mutate(
    tdiff_bin = cut(
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
  group_by(species_x, species_y, tdiff_bin) %>%
  summarise(
    bin_mean = mean(tdiff_value_x, na.rm = TRUE),
    mean_tdiff_y = mean(tdiff_value_y, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  ) %>%
  filter(count >= 1000)  # optional: remove bins with too few points

# Check result
head(tdiff_summary)


# ─── 6. Colour palette and shapes ───────────────────────────────────────────
cb_palette <- c(
  Oak    = "#E69F00",
  Beech  = "#0072B2",
  Spruce = "#009E73",
  Pine   = "#F0E442"
)

cb_shapes <- c(
  Oak    = 16,
  Beech  = 17,
  Spruce = 15,
  Pine   = 18
)

# ─── 7. Plot: four-panel tdiff–tdiff comparison ──────────────────────────────────
p <- ggplot(tdiff_summary, aes(
  x     = bin_mean,
  y     = mean_tdiff_y,
  color = species_y,
  shape = species_y
)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(aes(group = species_y), linewidth = 1) +  # ensure line connects points by species_y
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ species_x, ncol = 2, scales = "free") +
  scale_color_manual(name = "species", values = cb_palette) +
  scale_shape_manual(name = "species", values = cb_shapes) +
  labs(
    x = expression(paste("transpiration deficit of panel species (mm)")),
    y = expression(paste("mean transpiration deficit of each species (mm)")),
    title = ""
  ) +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title        = element_text(face = "bold", size = 14),
    axis.text         = element_text(color = "black", size = 12),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12),
    legend.position   = "bottom",
    legend.title = element_text(size = 12),
    legend.text       = element_text(size = 12)
  )

print(p)

# ─── Optional save ──────────────────────────────────────────────────────────
ggsave(
  filename = "results_rootzone/Figures_till2022/supplementary/TDiff_TDiff_pairwise_4panel.png",
  plot     = p,
  width    = 10,
  height   = 8,
  dpi      = 300
)
