setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

source("R/functions_month.R")

load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

plot_time_series_Quantiles_PSI_TDiff_year_month_species(final_df,
                                                        "results/Figures/All_month_year_species_Quantiles_PSI_TDiff_species.png")

plot_Quantiles_PSIbin_month_species(final_df,
                                    figure_output = "results/Figures/Quantiles_PSIbin_month_species.png")

plot_Quantiles_TDiffbin_month_species(final_df,
                                      figure_output = "results/Figures/Quantiles_TDiffbin_month_species.png")

plot_TDiff_PSIbin_month_species(final_df,
                                figure_output = "results/Figures/TDiff_PSIbin_month_species.png")

