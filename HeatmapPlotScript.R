# Heatmaps Plot Script

## Load Required Packages
library(tidyr)
library(ggplot2)
library(scales)

## Create Data Frame
# Example data frame - prostate cancer pathways
prostate_data <- data.frame(
  Canonical_Pathways = c(
    "Signaling by VEGF",
    "ESR-mediated signaling",
    "Semaphorin interactions",
    "SUMOylation of transcription cofactors",
    "MicroRNA Biogenesis Signaling Pathway",
    "Thrombin Signaling",
    "Signaling by Rho Family GTPases",
    "Prostate Cancer Signaling"
  ),
  Z_score = c(2.828, 2.449, 2.449, 2.236, 1.667, 2.449, 2.646, 2.236),
  Sample = "Prostate"
)

# Apply line breaks to long pathway names
wrap_pathways <- function(pathways, width = 25) {
  str_wrap(pathways, width = width)
}

prostate_data$Canonical_Pathways <- wrap_pathways(prostate_data$Canonical_Pathways)

# Create heatmap
ggplot(prostate_data, aes(x = Sample, y = Canonical_Pathways, fill = Z_score)) +
  geom_tile(color = "grey80") +
  scale_fill_gradientn(
    colors = c(
      "midnightblue", "#00008B", "#0000CD", "royalblue", "deepskyblue", "lightcyan",
      "white",
      "lightsalmon", "salmon1", "salmon", "firebrick1", "red2", "darkred"
    ),
    values = rescale(c(-3.5, -3, -2.8, -2.5, -2, -1.5, 0, 0.5, 1, 1.5, 2, 2.5, 3.5)),
    limits = c(-3.5, 3.5),
    na.value = "lightgrey"
  ) +
  labs(
    x = "Sample",
    y = "Pathways",
    fill = "Z-score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(face = "bold")
  )

# Save heatmap as a pdf
pdf("Prostate_pathway_heatmap.pdf", width = 10, height = 8)
print(last_plot())  
dev.off()
