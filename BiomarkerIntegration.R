# PFOS/PFOA Biomarker Analysis: Stouffer integration + data visualization

## Load required packages
library(ggplot2)
library(ggrepel)
library(dplyr)

## Functions defined:
# Stouffer integration method
stouffer <- function(x) {
  sum(x) / sqrt(length(x))
}

# Function to perform Stouffer integration and SD calculation
integrate_stouffer <- function(matrix_input, min_non_NA) {
  matrix_input[matrix_input == 'NaN'] <- NA
  matrix_filtered <- matrix_input[rowSums(is.na(matrix_input)) <= min_non_NA, ]
  
  stouffer_coeff <- apply(matrix_filtered, 1, function(x) stouffer(x[!is.na(x)]))
  sd_values <- apply(matrix_filtered, 1, function(x) sd(x[!is.na(x)]))
  
  result <- data.frame(
    stouffer_coefficient = stouffer_coeff,
    sd = sd_values,
    genes = rownames(matrix_filtered),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Load and prepare dataset

# Load dataset
df_raw <- read.csv("combined_data_FC.csv", row.names = 1)

# Reorder columns based on sample type, 
df_ordered <- df_raw[, c(
  "PFOA_monkey_stem_cells", "PFOS_testicular_tumor", 
  "PFOS_human_stem_cells", "PFOA_human_stem_cells", 
  "PFOS_human_liver", "PFOA_human_liver",
  "PFOA_mouse_liver_.0.3.", "PFOA_mouse_liver_.0.5.", 
  "PFOA_mouse_liver", "PFOS_prostate_tumor"
)]

## PFOS and PFOA combined integration + volcano plot

# Prepare subsets
monkeySC <- df_ordered[,1, drop = FALSE]
humanSC <- df_ordered[,3:4]
humanL  <- df_ordered[,5:6]
mouseL  <- df_ordered[,7:9]

# Create integrated matrix
mat_combined <- cbind(
  apply(humanL, 1, function(x) stouffer(x[!is.na(x)])),
  apply(humanSC, 1, function(x) stouffer(x[!is.na(x)])),
  monkeySC,
  apply(mouseL, 1, function(x) stouffer(x[!is.na(x)]))
)
colnames(mat_combined) <- c("human liver", "human stem cells", "monkey stem cells", "mouse liver")

# Perform integration
result_combined <- integrate_stouffer(mat_combined, min_non_NA = 2)

# Categorize expression
result_combined$expression <- "Not Significant"
result_combined$expression[result_combined$sd < 1.5 & result_combined$stouffer_coefficient < -1.5] <- "Downregulated & Low Sd"
result_combined$expression[result_combined$sd < 1.5 & result_combined$stouffer_coefficient > 2.5] <- "Upregulated & Low Sd"
result_combined$expression[result_combined$sd > 2 & result_combined$stouffer_coefficient < -1.5]  <- "Downregulated & High Sd"
result_combined$expression[result_combined$sd > 2 & result_combined$stouffer_coefficient > 2.5]  <- "Upregulated & High Sd"

# Top 10 genes
top_up_combined <- result_combined %>%
  filter(expression %in% c("Upregulated & Low Sd", "Upregulated & High Sd")) %>%
  arrange(desc(stouffer_coefficient)) %>%
  head(10)

top_down_combined <- result_combined %>%
  filter(expression %in% c("Downregulated & Low Sd", "Downregulated & High Sd")) %>%
  arrange(stouffer_coefficient) %>%
  head(10)

# Add labels
result_combined$label <- NA
result_combined$label[result_combined$genes %in% c(top_up_combined$genes, top_down_combined$genes)] <- 
  result_combined$genes[result_combined$genes %in% c(top_up_combined$genes, top_down_combined$genes)]

# Save result of top genes
write.csv(result_combined, "PFOS_PFOA_integration.csv", row.names = FALSE)

# Create volcano plot

png('PFOS_PFOA_volcano_plot.png', width = 8000, height = 5000, res = 600)

ggplot(data = result_combined, aes(x = stouffer_coefficient, y = sd, color = expression, label = label)) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      'Downregulated & Low Sd' = 'darkturquoise',
      'Upregulated & Low Sd' = 'darkorange',
      'Downregulated & High Sd' = 'royalblue',
      'Upregulated & High Sd' = 'red'
    ),
    na.value = "lightgray"
  ) +
  geom_text_repel(
    data = subset(result_combined, !is.na(label)),
    aes(label = paste0("bold('", label, "')")),
    box.padding = 0.5, max.overlaps = 18, colour = 'black', parse = TRUE
  ) +
  labs(color = "Gene Expression",
       x = "Integrated Fold Change (Stouffer Method)",
       y = "Standard Deviation") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  )

dev.off()

# ============================================================
# PFOA-only integration

# Prepare subset
mat_pfoa <- cbind(
  df_ordered[,1],  # monkey stem cells
  df_ordered[,4],  # human stem cells
  df_ordered[,6],  # human liver
  df_ordered[,7:9] # mouse liver doses
)
colnames(mat_pfoa) <- c("monkey stem cells", "human stem cells", "human liver", "mouse liver (0.3)", "mouse liver (0.05)", "mouse liver")

# Perform integration
result_pfoa <- integrate_stouffer(mat_pfoa, min_non_NA = 2)

# Save result
write.csv(result_pfoa, "PFOA_integration.csv", row.names = FALSE)

# ============================================================
# PFOS-only integration

# Prepare subset
mat_pfos <- cbind(
  df_ordered[,3], # human stem cells
  df_ordered[,5]  # human liver
)
colnames(mat_pfos) <- c("human stem cells", "human liver")

# Perform integration
result_pfos <- integrate_stouffer(mat_pfos, min_non_NA = 1)

# Save result
write.csv(result_pfos, "PFOS_integration.csv", row.names = FALSE)

# ============================================================
# Biomarker expression figures (calculate z-scores from stouffer integration and SD)

# Load biomarker dataset
biomarker_data <- read.csv("biomarkers.csv")

# Calculate Z-scores based on stouffer_coefficient
mean_stouffer <- mean(biomarker_data$stouffer_coefficient)
sd_stouffer <- sd(biomarker_data$stouffer_coefficient)

biomarker_data$z_score <- (biomarker_data$stouffer_coefficient - mean_stouffer) / sd_stouffer

# Filter genes within z-score thresholds
threshold <- 2
biomarker_data_filtered <- biomarker_data[abs(biomarker_data$z_score) <= threshold, ]

# Define colors for bar plot
biomarker_data_filtered$color <- ifelse(biomarker_data_filtered$z_score > 0, "red3", "royalblue3")

# Create biomarker bar plot
ggplot(biomarker_data_filtered, aes(x = reorder(genes, z_score), y = z_score, fill = color)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  labs(x = "Biomarkers", y = "Z-Score") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 12)
  )