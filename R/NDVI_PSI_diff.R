setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

library(dplyr)

pine_oak_locations <- all_results_df %>%
  filter(species %in% c("Beech", "Spruce")) %>%
  group_by(x, y) %>%
  filter(n_distinct(species) == 2) %>%
  ungroup()

head(pine_oak_locations)
