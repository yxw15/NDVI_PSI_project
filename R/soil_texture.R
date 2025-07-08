setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Create a data frame of KA5 (LUFA) soil texture classes with clay, silt, and sand percentage ranges
soil_data <- data.frame(
  code = c(
    "Ss",  "Su2",
    "St2", "Sl2", "Sl3", "Su3", "Su4",
    "Uu",  "Us",
    "Sl4", "Slu",
    "Uls", "Ut2", "Ut3",
    "Ut4", "Lu",
    "Ls2", "Ls3", "Ls4", "St3", "Lt2",
    "Ts3", "Ts4", "Lts",
    "Tu3", "Tu4",
    "Lt3",
    "Tu2", "Tl", "Ts2", "Tt"
  ),
  clay_lower  = c(
    0, 0,
    5, 5, 8, 0, 0,
    0, 0,
    12, 8,
    8, 8, 12,
    17, 17,
    17, 17, 17, 17, 25,
    35, 25, 25,
    30, 25,
    35,
    45, 45, 45, 65
  ),
  clay_upper  = c(
    5, 5,
    17, 8, 12, 8, 8,
    8, 8,
    17, 17,
    17, 12, 17,
    25, 30,
    25, 25, 25, 25, 35,
    45, 35, 45,
    45, 35,
    45,
    65, 65, 65, 100
    
  ),
  silt_lower  = c(
    0,  10,
    0, 10, 10, 25, 40,
    80, 50,
    10, 40,
    50, 65, 65,
    65, 50,
    40, 30, 15, 0, 30,
    0, 0, 15,
    50, 65,
    30,
    30, 15, 0, 0
  ),
  silt_upper  = c(
    10, 25,
    10, 25, 40, 40, 50,
    100, 80,
    40, 50,
    65, 88, 83,
    75, 65,
    50, 40, 30, 15, 50,
    15, 15, 30,
    65, 75,
    50,
    40, 30, 15, 35
  ),
  sand_lower  = c(
    85, 70,
    72, 67, 48, 52, 42,
    0, 12,
    43, 33,
    18, 0, 0,
    0, 0,
    25, 35, 45, 60, 15,
    40, 50, 25,
    0, 0,
    5,
    0, 5, 10, 0
  ),
  sand_upper  = c(
    100, 90,
    95, 85, 82, 75, 60,
    20, 50,
    78, 52,
    42, 27, 23,
    18, 33,
    43, 53, 68, 83, 35,
    65, 75, 60,
    20, 10,
    35,
    25, 40, 55, 35
  ),
  stringsAsFactors = FALSE
)

# Save the data frame to CSV
write.csv(soil_data, "soil_map/soil_texture_classes.csv", row.names = FALSE)
