# PFOS-PFOA-Exposure-Analysis

All the R codes used in the Uncovering the Tumorigenic Blueprint of PFOS and PFOA Through Multi-Organ Transcriptomic Analysis of Biomarkers, Pathways, and Therapeutic Targets study

Files:

1. DeSeq: RNA sequencing workflow (input files = raw count data from GEO) 
2. GSEA Analysis: Gene set enrichment analysis using KEGG database for each dataset + integration of pathways (input file = excel table S1)
3. Biomarker Integration: Stouffer integration method used to determine top DEGs + biomarker figures (input file = excel table S1)
4. Concentration-Dependent Analysis: Analysis and figures of concentration dependent response of PFOS and PFOA exposed human liver samples (input files = excel table S1) 
5. Heatmap Plot Script: Script used for creating all heatmaps in this study
6. Meta-Anlaysis_PvalComb: P-value combination (Stouffer + Fisher) using the metaRNAseq package. (input files = excel table S1)
