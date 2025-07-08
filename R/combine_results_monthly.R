setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

species_list <- c("Oak", "Beech", "Spruce", "Pine")
months <- c("April", "May", "June", "July", "August")
month_codes <- c("0423", "0523", "0623", "0723", "0823")

combined_df <- data.frame()

for (sp in species_list) {
  for (i in seq_along(months)) {
    month <- months[i]
    code <- month_codes[i]
    filepath <- paste0("results_monthly/", month, "/", sp, "/NDVI_PSI_TDiff_df_", code, ".RData")
    
    if (file.exists(filepath)) {
      load(filepath)  # loads species_df
      species_df$month <- month  # add month column
      combined_df <- rbind(combined_df, species_df)
    } else {
      warning(paste("File does not exist:", filepath))
    }
  }
}

head(combined_df)
save(combined_df, file = "results/Data/combined_depth50.RData")
