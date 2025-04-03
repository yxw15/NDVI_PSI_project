setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

#### Time series ####
plot_time_series_NDVI_PSI_TDiff_avg(all_results_df,
                                    "results/key_displays/time_series_Quantiles_PSI_TDiff.png")

plot_time_series_NDVI_PSI_TDiff_species_avg(all_results_df,
                                            "results/key_displays/time_series_Quantiles_PSI_TDiff_species.png")


#### NDVI Q - PSIbin ####

plot_NDVI_Q_PSIbin_exp_slope(all_results_df, 
                             "results/key_displays/NDVI_Q_PSIbin_exp_coeff.png",
                             "results/key_displays/NDVI_Q_PSIbin_exp_slope.png")

plot_NDVI_Q_PSIbin_exp_slope_negPSI(all_results_df, 
                                    "results/key_displays/NDVI_Q_PSIbin_exp_coeff_negPSI.png",
                                    "results/key_displays/NDVI_Q_PSIbin_exp_slope_negPSI.png")

#### NDVI Q - TDiffbin ####

plot_Quantiles_TDiff_linear_slope_coeff(all_results_df, 
                                        "results/key_displays/NDVI_Q_TDiffbin_linear_coeff.png",
                                        "results/key_displays/NDVI_Q_TDiffbin_linear_slope.png")

plot_Quantiles_TDiff_exp_slope_coeff(all_results_df, 
                                     "results/key_displays/NDVI_Q_TDiffbin_exp_coeff.png",
                                     "results/key_displays/NDVI_Q_TDiffbin_exp_slope.png")

plot_Quantiles_TDiff_exp_linear_slope_coeff(all_results_df,
                                            "results/key_displays/NDVI_Q_TDiffbin_linear_coeff.png",
                                            "results/key_displays/NDVI_Q_TDiffbin_exp_coeff.png",
                                            "results/key_displays/NDVI_Q_TDiffbin_exp_linear_slope.png")

plot_Quantiles_TDiff_poly_2_slope_coeff(all_results_df, 
                                        "results/key_displays/NDVI_Q_TDiffbin_poly_2_coeff.png",
                                        "results/key_displays/NDVI_Q_TDiffbin_poly_2_slope.png")

plot_Quantiles_TDiff_poly_3_slope_coeff(all_results_df, 
                                        "results/key_displays/NDVI_Q_TDiffbin_poly_3_coeff.png",
                                        "results/key_displays/NDVI_Q_TDiffbin_poly_3_slope.png")

#### TDiff - PSIbin ####

plot_TDiff_PSIbin_poly_2_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_2_coeff_each.png",
                               "results/key_displays/TDiff_PSIbin_poly_2_slope_each.png")

plot_TDiff_PSIbin_poly_3_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_3_coeff_each.png",
                               "results/key_displays/TDiff_PSIbin_poly_3_slope_each.png")


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


#### AIC ####
plot_NDVI_Q_PSIbin_AIC(all_results_df,  "results/key_displays/NDVI_Q_PSIbin_AIC.png")
plot_NDVI_Q_TDiffbin_AIC(all_results_df,  "results/key_displays/NDVI_Q_TDiffbin_AIC.png")
plot_TDiff_PSIbin_AIC(all_results_df,  "results/key_displays/TDiff_PSIbin_AIC.png")

plot_NDVI_Q_PSIbin_AIC_exp_linear(all_results_df,  "results/key_displays/NDVI_Q_PSIbin_AIC_exp_linear.png")
plot_NDVI_Q_TDiffbin_AIC_exp_linear(all_results_df,  "results/key_displays/NDVI_Q_TDiffbin_AIC_exp_linear.png")
