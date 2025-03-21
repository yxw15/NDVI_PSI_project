setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")


plot_NDVI_TDiff_poly_2_slope(all_results_df, 
                             # "results/key_displays/NDVI_Q_TDiffbin_poly_2_coeff.png",
                             "results/key_displays/NDVI_Q_TDiffbin_poly_2_slope.png")

plot_NDVI_TDiff_poly_3_slope(all_results_df, 
                             #"results/key_displays/NDVI_Q_TDiffbin_coeff.png",
                             "results/key_displays/NDVI_Q_TDiffbin_poly_3_slope.png")

plot_Quantiles_TDiff_exp_slope_coeff(all_results_df, 
                                     "results/key_displays/NDVI_Q_TDiffbin_exp_coeff.png",
                                     "results/key_displays/NDVI_Q_TDiffbin_exp_slope.png")

plot_Quantiles_TDiff_poly_2_slope_coeff(all_results_df, 
                                        "results/key_displays/NDVI_Q_TDiffbin_poly_2_coeff.png",
                                        "results/key_displays/NDVI_Q_TDiffbin_poly_2_slope.png")

plot_TDiff_PSIbin_slope(all_results_df, 
                        # "results/key_displays/TDiff_PSIbin_coeff.png",
                        "results/key_displays/TDiff_PSIbin_slope.png")

plot_TDiff_PSIbin_slope_each(all_results_df, 
                            "results/key_displays/TDiff_PSIbin_coeff_each.png",
                             "results/key_displays/TDiff_PSIbin_slope_each.png")

plot_NDVI_PSIbin_slope(all_results_df, 
                       "results/key_displays/NDVI_Q_PSIbin_coeff.png",
                       "results/key_displays/NDVI_Q_PSIbin_slope.png")


load("results/Data/All_Species_Proportions_PSI_TDiff.RData")

plot_NDVI_PDM_PSIbin_slope(all_results_df, 
                           "results/key_displays/NDVI_PDM_PSIbin_coeff.png",
                           "results/key_displays/NDVI_PDM_PSIbin_slope.png")

plot_Proportions_TDiff_poly_2_slope_coeff(all_results_df, 
                                          "results/key_displays/NDVI_PDM_TDiffbin_poly_2_coeff.png",
                                          "results/key_displays/NDVI_PDM_TDiffbin_poly_2_slope.png")

plot_TDiff_PSIbin_poly_3_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_3_ceoff.png",
                               "results/key_displays/TDiff_PSIbin_poly_3_slope.png")


plot_TDiff_PSIbin_poly_2_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_2_slope.png",
                               "results/key_displays/TDiff_PSIbin_poly_2_coeff.png")

source("R/functions_month.R")

load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

plot_time_series_Quantiles_PSI_TDiff_year_month_species(final_df,
                                                        "results/key_displays/Quantiles_PSI_TDiff_species_time_series.png")

plot_Quantiles_PSIbin_month_species(final_df,
                                    figure_output = "results/key_displays/Quantiles_PSIbin_month.png")

plot_Quantiles_PSIbin_month_species_linear_orig(final_df,
                                                figure_output = "results/key_displays/Quantiles_PSIbin_month_linear_orig.png",
                                                figure_output2 = "results/key_displays/Quantiles_PSIbin_month_parameters_orig.png")

plot_Quantiles_TDiffbin_month_species(final_df,
                                      figure_output = "results/key_displays/Quantiles_TDiffbin_month.png")

plot_TDiff_PSIbin_month_species(final_df,
                                figure_output = "results/key_displays/TDiff_PSIbin_month.png")


plot_Quantiles_TDiffbin_month_species_linear_agg(final_df,
                                                 figure_output = "results/key_displays/Quantiles_TDiffbin_month_linear_agg.png",
                                                 figure_output2 = "results/key_displays/Quantiles_TDiffbin_month_parameters_agg.png")


plot_Quantiles_TDiffbin_month_species_linear_orig(final_df,
                                                  figure_output = "results/key_displays/Quantiles_TDiffbin_month_linear_orig.png",
                                                  figure_output2 = "results/key_displays/Quantiles_TDiffbin_month_parameters_orig.png")

plot_TDiff_PSIbin_month_species_original(final_df,
                                         figure_output = "results/key_displays/TDiff_PSIbin_month_orig.png")


