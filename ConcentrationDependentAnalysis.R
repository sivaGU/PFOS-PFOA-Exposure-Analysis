# Concentration-Dependent Analysis 

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

## PFOA concentration analysis

# PFOA DEGs by concentration
deg_pfoa <- data.frame(
  Concentration = c(0.02, 0.1, 0.2, 1, 10, 20),
  DEGs = c(107, 101, 159, 115, 240, 276)
)

# Plot PFOA DEGs
ggplot(deg_pfoa, aes(x = Concentration, y = DEGs)) +
  geom_line(color = "blue", size = 1.5) +
  geom_point(color = "red", size = 3) +
  scale_x_log10(breaks = c(0.02, 0.1, 0.2, 1, 10, 20)) +
  labs(x = "PFOA Concentration (µM)", y = "Number of DEGs") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "plain"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

## PFOA pathway z-scores
data_pfoa <- data.frame(
  Pathway = c(
    "Mitochondrial Protein Degradation", 
    "Pathogen Induced Cytokine Storm Signaling Pathway", 
    "Molecular Mechanisms of Cancer", 
    "Serotonin Receptor Signaling", 
    "Pulmonary Fibrosis Idiopathic Signaling Pathway", 
    "Hepatic Fibrosis Signaling Pathway"
  ),
  Concentration_20 = c(3.464, -3.317, -2.985, -2.5, -2.309, -3),
  Concentration_10 = c(3.207, -1.508, -1.604, -1.155, -1.732, -0.775),
  Concentration_1  = c(2, -0.816, -0.333, -0.447, -0.447, -1.342),
  Concentration_0.2 = c(1, 0, 0.632, 1.134, 1.342, 1.134),
  Concentration_0.1 = c(2, -2.236, -0.816, -1.633, -2.646, -1.342),
  Concentration_0.02 = c(NA, NA, -0.707, 0, 0.447, NA)
)

# Convert to long format
data_pfoa_long <- data_pfoa %>%
  pivot_longer(cols = starts_with("Concentration"),
               names_to = "Concentration",
               values_to = "Z_score") %>%
  mutate(Concentration = as.numeric(gsub("Concentration_", "", Concentration)))

# Plot PFOA concentration–response curves
ggplot(data_pfoa_long, aes(x = Concentration, y = Z_score, color = Pathway, group = Pathway)) +
  geom_point(size = 1) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  scale_x_log10(
    breaks = c(0.02, 0.1, 0.2, 1, 10, 20),
    labels = c("0.02", "0.1", "0.2", "1", "10", "20")
  ) +
  labs(x = "PFOA Concentration (µM)", y = "Z-Score", color = "Pathway") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")

# Heatmap for PFOA pathways
ggplot(data_pfoa_long, aes(
  x = factor(Concentration, levels = sort(unique(Concentration))),
  y = str_wrap(Pathway, width = 30),
  fill = Z_score
)) +
  geom_tile(color = "lightgrey") +
  scale_fill_gradientn(
    colors = c('midnightblue', "blue3", 'mediumblue', 'royalblue3', 'deepskyblue',
               'white', 'orangered', 'firebrick1', 'red2', "red3", '#660000'),
    na.value = "lightgrey"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "PFOA Concentration (µM)", y = "Pathways", fill = "Z-Score") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.title.x = element_text(face = "plain", margin = margin(b = 10))
  )

## PFOS concentration analysis

# PFOS DEGs by concentration
deg_pfos <- data.frame(
  Concentration = c(0.02, 0.1, 0.2, 1, 10, 20),
  DEGs = c(79, 108, 220, 80, 543, 700)
)

# Plot PFOS DEGs
ggplot(deg_pfos, aes(x = Concentration, y = DEGs)) +
  geom_line(color = "blue", size = 1.5) +
  geom_point(color = "red", size = 3) +
  scale_x_log10(breaks = c(0.02, 0.1, 0.2, 1, 10, 20)) +
  labs(x = "PFOS Concentration (µM)", y = "Number of DEGs") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "plain"),
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank()
  )

# PFOS pathway z-scores
data_pfos <- data.frame(
  Pathway = c(
    "Mitochondrial Protein Degradation", 
    "Pathogen Induced Cytokine Storm Signaling Pathway", 
    "Molecular Mechanisms of Cancer", 
    "Serotonin Receptor Signaling", 
    "Pulmonary Fibrosis Idiopathic Signaling Pathway", 
    "Hepatic Fibrosis Signaling Pathway"
  ),
  Concentration_20 = c(-0.5, -1.461, -2.534, -2, -0.18, -1.512),
  Concentration_10 = c(1.604, -1.789, -1.567, -2.041, -0.218, -0.943),
  Concentration_1  = c(2.646, NA, NA, 1, -1, NA),
  Concentration_0.2 = c(1.342, -1.89, -2.357, -2.183, -0.277, -0.775),
  Concentration_0.1 = c(-2, NA, -1.342, -1, NA, 0),
  Concentration_0.02 = c(NA, -1, -1.667, NA, 0.447, -1.633)
)

# Convert to long format
data_pfos_long <- data_pfos %>%
  pivot_longer(cols = starts_with("Concentration"),
               names_to = "Concentration",
               values_to = "Z_score") %>%
  mutate(Concentration = as.numeric(gsub("Concentration_", "", Concentration)))

# Plot PFOS concentration–response curves
ggplot(data_pfos_long, aes(x = Concentration, y = Z_score, color = Pathway, group = Pathway)) +
  geom_point(size = 1) +
  geom_smooth(se = FALSE, method = "loess", span = 0.7) +
  scale_x_log10(
    breaks = c(0.02, 0.1, 0.2, 1, 10, 20),
    labels = c("0.02", "0.1", "0.2", "1", "10", "20")
  ) +
  labs(x = "PFOS Concentration (µM)", y = "Z-Score", color = "Pathway") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")

# Heatmap for PFOS pathways
ggplot(data_pfos_long, aes(
  x = factor(Concentration, levels = sort(unique(Concentration))),
  y = str_wrap(Pathway, width = 30),
  fill = Z_score
)) +
  geom_tile(color = "lightgrey") +
  scale_fill_gradientn(
    colors = c('midnightblue', "blue3", 'mediumblue', 'royalblue3', 'deepskyblue',
               'white', 'orangered', 'firebrick1', 'red2', "red3", '#660000'),
    na.value = "lightgrey"
  ) +
  scale_x_discrete(position = "top") +
  labs(x = "PFOS Concentration (µM)", y = "Pathways", fill = "Z-Score") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.title.x = element_text(face = "plain", margin = margin(b = 10))
  )
