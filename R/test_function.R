library(tidyverse)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

NDVI_PSIbin_df <- NDVI_PSIbin(all_results_df)
TDiff_PSIbin_df <- TDiff_PSIbin(all_results_df)

