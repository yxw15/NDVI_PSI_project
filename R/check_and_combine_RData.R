# 1) Set working directory and configuration
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Define month-specific settings
months_config <- list(
  April  = list(day = "04-23", NDVI = "../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(day = "05-25", NDVI = "../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(day = "06-26", NDVI = "../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)
months  <- names(months_config)
species <- c("Beech", "Oak", "Spruce", "Pine")
depths  <- c(100, 150)

# 2) Prepare combinations of month, species, depth
combos <- expand.grid(
  month   = months,
  species = species,
  depth   = depths,
  stringsAsFactors = FALSE
)
# Attach day-of-year string for file naming
combos$day <- vapply(combos$month,
                     function(m) months_config[[m]]$day,
                     character(1))

# Helper to build file paths for each combo
make_path <- function(month, species, depth, day) {
  fname <- sprintf("NDVI_PSI_TDiff_%s_depth%d.RData",
                   gsub("-", "", day), depth)
  file.path(sprintf("results_monthly_%d", depth),
            month, species, fname)
}
combos$path <- with(combos,
                    mapply(make_path, month, species, depth, day, USE.NAMES = FALSE)
)

# 3) Ensure output directory exists
out_dir <- "results/Data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 4) Loop over depths, load each file, add month, combine and save
for (d in depths) {
  message("Processing depth = ", d, " cm")
  
  # Filter for current depth
  subset_combos <- subset(combos, depth == d)
  dfs <- vector("list", nrow(subset_combos))
  
  # Load each .RData and add 'month' column
  for (i in seq_len(nrow(subset_combos))) {
    p <- subset_combos$path[i]
    if (!file.exists(p)) {
      warning("Missing file, skipping: ", p)
      next
    }
    tmp_env <- new.env()
    load(p, envir = tmp_env)
    df <- tmp_env$species_df
    # Add a 'month' column from the combos table
    df$month <- subset_combos$month[i]
    dfs[[i]] <- df
  }
  
  # Combine all data.frames and save
  combined_df <- do.call(rbind, dfs)
  out_file <- file.path(out_dir,
                        sprintf("all_species_months_depth%d.RData", d))
  save(combined_df, file = out_file)
  message("  ➜ Saved combined data to ", out_file,
          " (", nrow(combined_df), " rows )")
}
