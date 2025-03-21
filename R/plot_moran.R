setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project/results_spatial")

# Load required libraries
library(ggplot2)
library(dplyr)

# Read CSV files for each parameter at two distances
quantiles_0.2 <- read.csv("Quantiles_spatial_0.2.csv")
quantiles_0.4 <- read.csv("Quantiles_spatial_0.4.csv")
psi_0.2       <- read.csv("PSI_spatial_0.2.csv")
psi_0.4       <- read.csv("PSI_spatial_0.4.csv")
tdiff_0.2     <- read.csv("TDiff_spatial_0.2.csv")
tdiff_0.4     <- read.csv("TDiff_spatial_0.4.csv")

# Add a new column indicating the parameter for each dataset
quantiles_0.2$Parameter <- "Quantiles"
quantiles_0.4$Parameter <- "Quantiles"
psi_0.2$Parameter       <- "PSI"
psi_0.4$Parameter       <- "PSI"
tdiff_0.2$Parameter     <- "TDiff"
tdiff_0.4$Parameter     <- "TDiff"

# Combine all data frames into one
combined_data <- rbind(quantiles_0.2, quantiles_0.4,
                       psi_0.2, psi_0.4,
                       tdiff_0.2, tdiff_0.4)

# Ensure 'distance' is numeric
combined_data$distance <- as.numeric(as.character(combined_data$distance))

# Create a line plot using ggplot2
ggplot(combined_data, aes(x = distance, y = Moran_I, color = Parameter, group = Parameter)) +
  geom_line() +
  geom_point() +
  labs(title = "Moran's I vs Distance for Three Parameters",
       x = "Distance (degrees)",
       y = "Moran's I") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
ggsave("moran_distance.png")
