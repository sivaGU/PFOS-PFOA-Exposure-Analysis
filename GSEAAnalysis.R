# KEGG GSEA Analysis for PFOS/PFOA Exposed Samples

## Load Required Packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(biomaRt)
library(KEGGREST)
library(EnrichmentBrowser)
library(stringr)
library(reshape2)

## Define Organism 
# used: org.Hs.eg.db (human), org.Mm.eg.db (mouse), org.Mmu.eg.db (monkey)
# set organism database
organism_db <- org.Hs.eg.db  

#  Load and Prepare Data
# function to clean DESeq2 result data
clean_deseq_data <- function(filepath, columns_to_remove) {
  df <- read.csv(filepath, header = TRUE)
  df <- df[, !(names(df) %in% columns_to_remove)]
  return(df)
}

# Open  file for analysis
# example for PFOS human liver:
df <- clean_deseq_data("PFOSHumanliverddsres.csv", c("baseMean", "lfcSE", "stat", "padj"))

# Prepare gene list
colnames(df) <- c("gene_symbol", "log2", "pvalue")
original_gene_list <- df$log2
names(original_gene_list) <- df$gene_symbol
gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)

# Map gene symbols to ENTREZ IDs
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism_db)
ids <- ids[!duplicated(ids$SYMBOL), ]

# Merge with original data
df2 <- df[df$gene_symbol %in% ids$SYMBOL, ]
df2$ENTREZID <- ids$ENTREZID

# Create kegg_gene_list
kegg_gene_list <- df2$log2
names(kegg_gene_list) <- df2$ENTREZID
kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)

## Run GSEA: gseKEGG
# Set KEGG organism code: "hsa" (human), "mmu" (mouse), "mcc" (monkey)
kegg_organism <- "hsa"

kk2 <- gseKEGG(
  geneList = kegg_gene_list,
  organism = kegg_organism,
  minGSSize = 2,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  keyType = "ncbi-geneid"
)

# Save GSEA results
gsea_results <- as.data.frame(kk2@result) %>% select(ID, Description, NES, p.adjust)
write.csv(gsea_results, "GSEA_results.csv", row.names = FALSE)

## Visualization
# dotplot
dotplot(kk2, showCategory = 10, title = "Enriched Pathways", split = ".sign") +
  facet_grid(. ~ .sign)

# CNET plot
cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)

# ridge plot
ridgeplot(kk2) + labs(x = "Enrichment Distribution")

## Merge the multiple GSEA result files for each dataset

# Merge CSVs
file_list <- list.files(pattern = "GSEA_results.csv")
merged_data <- Reduce(function(x, y) merge(x, y, by = "Description", all = TRUE),
                      lapply(file_list, read.csv))

# Keep pathways present in 2+ datasets
pathway_counts <- rowSums(!is.na(merged_data[, -1]))
filtered_data <- merged_data[pathway_counts >= 2, ]

write.csv(filtered_data, "Merged_GSEA_Results.csv", row.names = FALSE)

## Heatmap of Selected Pathways

# example: select pathways by category
selected_pathways <- c(
  "Chemical carcinogenesis - DNA adducts", "Wnt signaling pathway", "Hippo signaling pathway",
  "Arachidonic acid metabolism", "Fatty acid metabolism", "Peroxisome",
  "Graft-versus-host disease", "Systemic lupus erythematosus", "Autoimmune thyroid disease",
  "PPAR signaling pathway", "Linoleic acid metabolism", "alpha-Linolenic acid metabolism",
  "Protein processing in endoplasmic reticulum", "Carbon metabolism"
)

# Filter data
selected_data <- filtered_data[filtered_data$Description %in% selected_pathways, ]

# Melt data for heatmap
heatmap_data <- melt(selected_data, id.vars = "Description", variable.name = "Sample", value.name = "NES")

# Format long names
heatmap_data$Description <- str_wrap(heatmap_data$Description, width = 30)

# Heatmap plot
ggplot(heatmap_data, aes(x = Sample, y = Description, fill = NES)) +
  geom_tile(color = "grey80") +
  scale_fill_gradientn(
    colors = c('midnightblue', "blue3", 'mediumblue', 'royalblue3', 'deepskyblue', 
               'white', 'orangered', 'firebrick1', 'red2', "red3", '#660000'),
    limits = c(-2, 2),
    na.value = "lightgrey"
  ) +
  labs(
    x = "Samples",
    y = "Pathways",
    fill = "NES"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(face = "bold")
  ) +
  scale_x_discrete(
    labels = c(
      "NES_PFOA.0.05.ML" = "PFOA_ML(0.05)", 
      "NES_PFOA.0.3.ML" = "PFOA_ML(0.3)", 
      "NES_PFOAHeSC" = "PFOA_HeSC",
      "NES_PFOAHL" = "PFOA_HL",
      "NES_MMeSC" = "Monkey_SC",
      "NES_PFOAML" = "Monkey_Liver",
      "NES_PFOSHeSC" = "PFOS_HeSC",
      "NES_PFOSHL" = "PFOS_HL"
    )
  )

# Save heatmap
pdf("GSEA_heatmap.pdf", width = 10, height = 8)
print(last_plot())
dev.off()
