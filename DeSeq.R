# RNA-Seq Data Processing and DESeq2 Differential Expression Analysis

# Load required packages
library(DESeq2)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(pathview)
library(enrichplot)
library(limma)
library(ggplot2)
library(tidyverse)
library(airway)
library(pheatmap)
library(msigdbr)

### Prepare Counts and Sample Info for Each Dataset

## PFOS Human Embryonic Stem Cells (HeSC) Data
counts_data <- read.csv("HeSCRawGeneCounts.csv", row.names = 1)
counts_data <- counts_data[, c("read_count_E6_D16_DMSO_1", "read_count_E6_D16_DMSO_2",
                               "read_count_E6_D16_PFOS_1", "read_count_E6_D16_PFOS_2")]
colData <- read.csv("GSE236956_sample_info.csv", row.names = 1)
## --------------------------

## PFOS Testes Tumor Data
counts_data <- read.csv("PFOSTestesTumor_raw_counts.csv", row.names = 1)

# Convert Entrez IDs to Gene Symbols
entrez_ids <- rownames(counts_data)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = entrez_ids, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

# Add gene symbols
counts_data$gene_symbol <- gene_symbols

# Optionally, filter out genes without symbols
counts_data <- counts_data[!is.na(rownames(counts_data)), ]

# Reload counts with gene symbols
counts_data <- read.csv("PFOSTestesTumor.csv", row.names = 1)
counts_data <- counts_data[, c("GSM8158094", "GSM8158096", "GSM8158097", "GSM8158098",
                               "GSM8158099", "GSM8158100")]
colData <- read.csv("TestesTumorSampleInfo.csv", row.names = 1)
## --------------------------

## PFOA Mouse Liver (GSE119441)
counts_data <- read.csv("9441PfoaMouseLiver_raw_counts.csv", row.names = 1)
colData <- read.csv("9441pfoamouseliversampleinfo.csv", row.names = 1)
## --------------------------

## PFOA Mouse Liver (GSE212294, multiple doses)
counts_data <- read.csv("GSE212294data.csv", row.names = 1)

# Convert Ensembl IDs to gene symbols (mouse)
ensembl_ids <- rownames(counts_data)
gene_symbols <- mapIds(org.Mm.eg.db, 
                       keys = ensembl_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

counts_data$gene_symbol <- gene_symbols

# Load counts mapped to gene symbols
counts_data <- read.csv("PFOAmouseliverrawcountsGSE212294.csv", row.names = 1)

# Subset counts for specific conditions (choose based on sample set)
# Example: PFOA 0.3 mg/kg WT
counts_data <- counts_data[, c("PFASdil4", "PFASdil63", "PFASdil84", "PFASdil86",
                               "PFASdil16", "PFASdil80", "PFASdil95", "PFASdil96")]

# Load corresponding sample info
colData <- read.csv("0.3pfoaWTmlsampleinfo.csv", row.names = 1)
## --------------------------

## PFOS/PFOA Human Liver Data
counts_data <- read.csv("PFOS_PFOAHumanLiverRawCounts.csv", row.names = 1)

# Subset counts: PFOA samples
counts_data <- counts_data[, c("DMSO_0_657", "DMSO_0_658", "DMSO_0_659", "DMSO_0_660", 
                               "DMSO_0_651", "DMSO_0_652", "DMSO_0_655", "DMSO_0_656", 
                               "DMSO_0_747", "DMSO_0_751", "DMSO_0_649", "DMSO_0_650", 
                               "DMSO_0_653", "DMSO_0_654", "DMSO_0_745", "DMSO_0_746", 
                               "DMSO_0_749", "DMSO_0_750", "PFOA_20_595", "PFOA_20_599",
                               "PFOA_20_691")]

# Subset counts: PFOS samples
counts_data <- counts_data[, c("DMSO_0_657", "DMSO_0_658", "DMSO_0_659", "DMSO_0_660",
                               "DMSO_0_651", "DMSO_0_652", "DMSO_0_655", "DMSO_0_656",
                               "DMSO_0_747", "DMSO_0_751", "DMSO_0_649", "DMSO_0_650",
                               "DMSO_0_653", "DMSO_0_654", "DMSO_0_745", "DMSO_0_746",
                               "DMSO_0_749", "DMSO_0_750", "PFOS_20_594", "PFOS_20_598",
                               "PFOS_20_690", "PFOS_20_694")]

# Load sample information
colData <- read.csv("PFOAHumanLiverSampleInfo.csv", row.names = 1)
# or
colData <- read.csv("PFOSHumanLiverSampleInfo.csv", row.names = 1)
## --------------------------

## DESeq2 Analysis

# Check for NA values
sum(is.na(counts_data))

# Remove rows with any NA
counts_data <- na.omit(counts_data)

# Check for non-integer counts
sum(floor(counts_data) != counts_data)

# Round counts if necessary
counts_data <- round(counts_data)

# Check sample matching
if (!all(colnames(counts_data) == rownames(colData))) {
  stop("Mismatch between counts data columns and sample info rows.")
}

if (ncol(counts_data) != nrow(colData)) {
  stop("Number of samples mismatch between counts data and sample info.")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ PFAS)

# Filter genes with low counts (less than 10 across samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set reference level for PFAS variable
dds$PFAS <- relevel(dds$PFAS, ref = "Control")

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# View first few results
head(res)

# Save results
write.csv(res, "datasetname_dds.csv", row.names = TRUE)
