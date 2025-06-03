setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Create a data frame of KA5 (LUFA) soil texture classes with clay, silt, and sand percentage ranges
soil_data <- data.frame(
  code = c(
    "Ss",  "Sl2", "Sl3", "Sl4", "Slu", "St2", "St3",
    "Uu",  "Us",  "Ut2", "Ut3", "Ut4", "Uls",
    "Lu",  "Ls2", "Ls3", "Ls4", "Lt2", "Lt3", "Lts",
    "Ts2", "Ts3", "Ts4", "Tu2", "Tu3", "Tu4",
    "Tl",  "Tt"
  ),
  clay_lower  = c(
    0,   5,   8,  12,   8,   5,  17,
    0,   0,   8,  12,  17,   8,
    12,   5,   8,  12,  17,  25,  15,
    8,  12,  17,  25,  40,  50,
    40,  65
  ),
  clay_upper  = c(
    5,   8,  12,  17,  17,  17,  25,
    8,   8,  12,  17,  25,  17,
    17,   8,  12,  17,  25,  40,  30,
    12,  17,  25,  40,  50,  65,
    65, 100
  ),
  silt_lower  = c(
    0,  10,  10,  10,  40,   0,   0,
    80,  50,  65,  65,  65,  50,
    30,  30,  25,  18,  18,  15,  40,
    40,  40,  40,  30,  25,  25,
    15,   0
  ),
  silt_upper  = c(
    10,  25,  40,  40,  50,  10,  15,
    100,  80,  92,  88,  83,  65,
    55,  45,  40,  35,  35,  30,  60,
    60,  60,  60,  50,  50,  50,
    35,  35
  ),
  sand_lower  = c(
    85,  67,  48,  43,  33,  73,  60,
    0,  12,   0,   0,   0,  18,
    13,  42,  40,  42,  30,  10,   5,
    5,   3,   3,   0,   0,   0,
    0,   0
  ),
  sand_upper  = c(
    100,  85,  82,  78,  52,  95,  83,
    20,  50,  27,  23,  18,  42,
    58,  67,  60,  58,  58,  55,  25,
    25,  25,  15,   5,  10,  10,
    5,  35
  ),
  stringsAsFactors = FALSE
)

# Save the data frame to CSV
write.csv(soil_data, "soil_map/soil_texture_classes.csv", row.names = FALSE)
