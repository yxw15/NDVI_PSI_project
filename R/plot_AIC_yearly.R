setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)
library(purrr)

# ----------------------------
# 1) Data prep
# ----------------------------
data <- combined %>%
  filter(month %in% c("July", "August"),
         year < 2023,
         Quantiles > 0) %>%
  drop_na()

# ----------------------------
# 2) Themes + colors
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
# 3) Annual summaries
# ----------------------------
summary_df <- data %>%
  group_by(species, year) %>%
  summarise(
    mean_ndvi = mean(Quantiles, na.rm = TRUE),
    mean_soil_water_potential = mean(soil_water_potential, na.rm = TRUE),
    mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE),
    .groups = "drop"
  )

summary_long <- summary_df %>%
  pivot_longer(
    cols = c(mean_ndvi, mean_soil_water_potential, mean_transpiration_deficit),
    names_to = "variable",
    values_to = "value"
  )

# ----------------------------
# 4) ACF per species × variable (lags 1..5)
# ----------------------------
acf_df <- summary_long %>%
  arrange(species, variable, year) %>%
  group_by(species, variable) %>%
  summarise(
    acf_obj = list(acf(value, lag.max = 7, plot = FALSE)),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    acf_crit = 1.96 / sqrt(n),
    acf_tbl = map(acf_obj, ~ tibble(
      lag = as.numeric(.x$lag[-1]),
      autocorrelation = as.numeric(.x$acf[-1])
    ))
  ) %>%
  select(-acf_obj) %>%
  unnest(acf_tbl)

# optional: best lag (max |ACF|)
best_lag_df <- acf_df %>%
  group_by(species, variable) %>%
  slice_max(order_by = abs(autocorrelation), n = 1, with_ties = FALSE) %>%
  ungroup()


# ----------------------------
# 5) Plot: ACF with species colors + significance as red dashed lines
# ----------------------------
# var_map_expr <- c(
#   mean_ndvi = "NDVI[quantiles]",
#   mean_soil_water_potential = "Psi[~~soil]",
#   mean_transpiration_deficit = "T[d]"
# )

var_map_expr <- c(
  mean_ndvi = "NDVI quantiles (rank)",
  mean_soil_water_potential = "soil water potential (kPa)",
  mean_transpiration_deficit = "transpiration deficit (mm)"
)

acf_df$species <- factor(
  acf_df$species,
  levels = c("Oak", "Beech", "Spruce", "Pine")
)

p_acf <- ggplot(acf_df, aes(lag, autocorrelation, fill = species)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(aes(yintercept =  acf_crit), linetype = "dashed", color = "blue") +
  geom_hline(aes(yintercept = -acf_crit), linetype = "dashed", color = "blue") +
  geom_col(width = 0.8) +
  facet_grid(
    variable ~ species,
    scales = "free_y",
    labeller = labeller(variable = as_labeller(var_map_expr))
    # labeller = labeller(variable = as_labeller(var_map_expr, label_parsed))
  ) +
  scale_fill_manual(values = sp_cols, guide = "none") +
  scale_x_continuous(
    breaks = 1:7,
    # limits = c(0.5, 7.5)   # keeps bars centered and avoids clipping
  ) +
  labs(x = "lag (years)", y = "temporal autocorrelation") +
  base_theme +
  theme(
    strip.text.y = element_text(size=14, angle = -90),
  )

p_acf

ggsave(
  filename = "results_rootzone/Figures_till2022/SI_PSI/SI_ACF_lag_7.png",
  plot = p_acf,
  width = 12,
  height = 10,
  dpi = 300
)

