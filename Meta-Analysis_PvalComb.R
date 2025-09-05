## OpenPackages
library(metaRNASeq)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(clusterProfiler)
library(msigdbr)

# for reproducible label placement in ggrepel
set.seed(42)  

## define function of signed z-score from pvalue, will be used later:
signed_z_from_p <- function(p_two_sided, sign_dir){
  # convert two-sided p + direction to a signed Z
  z_mag <- qnorm(p = pmin(0.999999, p_two_sided/2), lower.tail = FALSE)
  z_mag * sign(sign_dir)
}

## Load each study differential expression results (gene, fold change, pval)
human_liver_PFOS <- read.csv("PFOS Human Liver.csv")
human_liver_PFOA <- read.csv("PFOA Human Liver.csv")
human_ESC_PFOS   <- read.csv("PFOS Human Stem Cells.csv")
human_ESC_PFOA   <- read.csv("PFOA Human Stem Cells.csv")
mouse_liver_PFOA <- read.csv("PFOA Mouse Liver.csv")
mouse_liver_PFOA_005 <- read.csv("PFOA Mouse Liver 0.05.csv")
mouse_liver_PFOA_03  <- read.csv("PFOA Mouse Liver 0.3.csv")

## Group datasets for meta-analysis arms
# PFOS (human)
studies_PFOS_human <- list(
  human_liver_PFOS = human_liver_PFOS,
  human_ESC_PFOS   = human_ESC_PFOS)

# PFOA (human)
studies_PFOA_human <- list(
  human_liver_PFOA = human_liver_PFOA,
  human_ESC_PFOA   = human_ESC_PFOA)

# PFOA (mouse)
studies_PFOA_mouse <- list(
  mouse_liver_PFOA     = mouse_liver_PFOA,
  mouse_liver_PFOA_005 = mouse_liver_PFOA_005,
  mouse_liver_PFOA_03  = mouse_liver_PFOA_03)

## Weights by number of biological replicates (per study)
# if PFAS and control sizes differ, calculate harmonic mean per study
nrep_PFOS_human  <- c(human_liver_PFOS=7, human_ESC_PFOS=4)
nrep_PFOA_human  <- c(human_liver_PFOA=5, human_ESC_PFOA=4)
nrep_PFOA_mouse  <- c(mouse_liver_PFOA=10, mouse_liver_PFOA_005=8, mouse_liver_PFOA_03=8)

## Convert two-sided p to one-sided using fold-change sign
one_sided <- function(p_two_sided, log2FC){
  ifelse(log2FC >= 0, p_two_sided/2, 1 - p_two_sided/2)}

## Harmonize per-study DE tables to (gene, log2FC, pvalue)
harmonize_de <- function(df, study_name){
  nm  <- names(df); nml <- tolower(nm)
  
  # locate gene column (many common variants supported)
  gene_idx <- which(nml %in% c("gene","genes","geneid","gene_id","symbol","gene.symbol","genename","gene.name"))[1]
  if (is.na(gene_idx)) gene_idx <- if (is.character(df[[1]])) 1 else NA_integer_
  if (is.na(gene_idx)) stop(sprintf("[%s] Can't find a gene column; saw: %s", study_name, paste(nm, collapse=", ")))
  gene_col <- nm[gene_idx]
  
  # locate fold-change/log2FC column; convert linear FC to log2 if needed
  lfc_idx <- grep("(log2|log_fc|logfc|log2fold)", nml, perl = TRUE)
  if (length(lfc_idx) == 0) lfc_idx <- grep("(fold).*?(change)|(^|[._-])fc($|[._-])", nml, perl = TRUE)
  if (length(lfc_idx) == 0) stop(sprintf("[%s] Can't find fold-change column; saw: %s", study_name, paste(nm, collapse=", ")))
  lfc_col <- nm[lfc_idx[1]]
  
  # p-value or stat (Wald/t) column
  p_idx  <- grep("(p\\.?(value)?|pval)", nml, perl = TRUE)
  stat_i <- grep("(^|[._-])(stat|wald|z|t)($|[._-])", nml, perl = TRUE)
  have_p <- length(p_idx) > 0; have_s <- length(stat_i) > 0
  if (!have_p && !have_s) stop(sprintf("[%s] Need pvalue or stat; columns: %s", study_name, paste(nm, collapse=", ")))
  
  gene_vec <- df[[gene_col]]
  fc_raw   <- suppressWarnings(as.numeric(df[[lfc_col]]))
  
  # heuristically convert linear FC to log2FC if values look like ratios
  log2FC <- fc_raw
  if (!any(is.na(fc_raw))) {
    if (all(fc_raw > 0) && stats::median(fc_raw, na.rm = TRUE) > 0.25 && stats::median(fc_raw, na.rm = TRUE) < 4) {
      eps <- 1e-12
      log2FC <- log2(pmax(fc_raw, eps))}}
  
  # two-sided p from pvalue or from test statistic
  if (have_p) {
    p_two <- suppressWarnings(as.numeric(df[[nm[p_idx[1]]]]))} else {
    statv <- suppressWarnings(as.numeric(df[[nm[stat_i[1]]]]))
    p_two <- 2*pnorm(-abs(statv))}
  
  # keep best (lowest p) row per gene
  out <- tibble::tibble(
    gene   = as.character(gene_vec),
    log2FC = as.numeric(log2FC),
    pvalue = as.numeric(p_two)) |>
    dplyr::filter(!is.na(.data$gene), !is.na(.data$log2FC), !is.na(.data$pvalue), is.finite(.data$pvalue)) |>
    dplyr::group_by(.data$gene) |>
    dplyr::slice_min(order_by = .data$pvalue, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  if (nrow(out) == 0) stop(sprintf("[%s] No valid rows after harmonization.", study_name))
  message(sprintf("[%s] Harmonized: %d genes", study_name, nrow(out)))
  out}

## Meta-analysis runner (Stouffer + Fisher) with replicate weights
meta_run <- function(studies_list, nrep_vec, min_k = 2){
  # harmonize each study
  studies_list <- purrr::imap(studies_list, ~ harmonize_de(.x, .y))
  
  # long table with one-sided p
  long <- purrr::imap_dfr(studies_list, ~ .x |>
                            dplyr::mutate(study = .y,
                                          p1 = one_sided(.data$pvalue, .data$log2FC)))
  
  # wide p-value matrix; keep genes present in >= min_k studies
  pwide <- long |>
    dplyr::select(.data$gene, .data$study, .data$p1) |>
    tidyr::pivot_wider(names_from = .data$study, values_from = .data$p1) |>
    dplyr::mutate(n_non_na = rowSums(!is.na(dplyr::across(-.data$gene)))) |>
    dplyr::filter(.data$n_non_na >= min_k) |>
    dplyr::select(-.data$n_non_na)
  
  genes <- pwide$gene
  if (length(genes) == 0) stop("No genes present in >= min_k studies. Lower min_k or check inputs.")
  
  # collect per-study vectors with matching weights
  study_names <- names(studies_list)
  indpval <- lapply(study_names, function(nm) pwide[[nm]])
  names(indpval) <- study_names
  keep <- which(sapply(indpval, function(v) any(!is.na(v))))
  indpval <- indpval[keep]; study_names <- names(indpval)
  
  if (!all(study_names %in% names(nrep_vec))) {
    stop(sprintf("Weights missing for: %s", paste(setdiff(study_names, names(nrep_vec)), collapse=", ")))
  }
  nrep_use <- unname(nrep_vec[study_names])
  
  # metaRNASeq: Stouffer (invnorm) + Fisher, BH FDR
  st <- metaRNASeq::invnorm(indpval = indpval, nrep = nrep_use, BHth = 0.999)
  fi <- metaRNASeq::fishercomb(indpval = indpval, BHth = 0.999)
  
  # directional votes across studies for each gene
  dir_wide <- long |>
    dplyr::mutate(dir = sign(.data$log2FC)) |>
    dplyr::select(.data$gene, .data$study, .data$dir) |>
    tidyr::pivot_wider(names_from = .data$study, values_from = .data$dir)
  
  vote <- dir_wide |>
    dplyr::mutate(
      n_up   = rowSums(dplyr::across(-.data$gene, ~ .x ==  1), na.rm = TRUE),
      n_down = rowSums(dplyr::across(-.data$gene, ~ .x == -1), na.rm = TRUE)
    ) |>
    dplyr::select(.data$gene, .data$n_up, .data$n_down)
  
  # final meta table per arm
  tibble::tibble(
    gene               = genes,
    p_meta_stouffer    = st$rawpval,
    FDR_meta_stouffer  = p.adjust(st$rawpval, method = "BH"),
    p_meta_fisher      = fi$rawpval,
    FDR_meta_fisher    = p.adjust(fi$rawpval, method = "BH")) |>
    dplyr::left_join(vote, by = "gene") |>
    dplyr::arrange(.data$p_meta_stouffer)}

## Run meta for each arm
meta_PFOS_human  <- meta_run(studies_PFOS_human,  nrep_PFOS_human)
meta_PFOA_human  <- meta_run(studies_PFOA_human,  nrep_PFOA_human)
meta_PFOA_mouse  <- meta_run(studies_PFOA_mouse,  nrep_PFOA_mouse)

## Save per-arm meta tables (in supplementary)
write.csv(meta_PFOS_human,  "Supp_meta_PFOS_human.csv",  row.names = FALSE)
write.csv(meta_PFOA_human,  "Supp_meta_PFOA_human.csv",  row.names = FALSE)
write.csv(meta_PFOA_mouse,  "Supp_meta_PFOA_mouse.csv",  row.names = FALSE)

##############################################
## Monkey single-study processing (included alongside meta arms)
monkey_DE <- read.csv("PFOA Monkey Stem Cells.csv", check.names = FALSE)

# harmonize to (gene, log2FC, pvalue)
monkey_std <- harmonize_de(monkey_DE, "monkey_ESC_PFOA")

# add FDR and signed Z; keep minimal columns
monkey_DE2 <- monkey_std %>%
  dplyr::mutate(
    FDR      = p.adjust(pvalue, method = "BH"),
    dir_score = sign(log2FC),
    p         = pvalue,
    signed_z  = signed_z_from_p(pvalue, dir_score),
    group     = "PFOA – Monkey (single study)"
  ) %>%
  dplyr::select(gene, p, FDR, signed_z, group)

write.csv(monkey_DE2, "mokeysinglestudyprocessing.csv")

##############################################
## Robust DESeq2 consensus table (count across individual datasets)
df <- readr::read_csv("Deseq2data.csv")
colnames(df)[1] <- "gene"  # ensure first column is 'gene'

## create a summary for each gene in the study 
# calculate avg FC,median FC,  assign overall direction of FC, then rank the genes
gene_summary <- df %>%
  rowwise() %>%
  mutate(
    n_datasets = sum(!is.na(c_across(-gene))),  
    mean_FC    = mean(c_across(-gene), na.rm = TRUE),
    median_FC  = median(c_across(-gene), na.rm = TRUE),
    direction  = ifelse(mean_FC > 0, "Up", "Down")) %>%
  ungroup() %>%
  arrange(desc(n_datasets))

#filter to retain genes only available in >4 datasets
robust_genes <- gene_summary %>% filter(n_datasets >= 4)
write.csv(robust_genes, "robust_significant_genes4.csv", row.names = FALSE)

##############################################
## Combine meta results across arms + monkey
meta_files <- list(
  PFOS_human = "Supp_meta_PFOS_human.csv",
  PFOA_human = "Supp_meta_PFOA_human.csv",
  PFOA_mouse = "Supp_meta_PFOA_mouse.csv")

meta_list <- purrr::imap(meta_files, ~ read.csv(.x) %>% mutate(dataset = .y))
meta_long <- bind_rows(meta_list)

# keep FDR<0.05 per arm
meta_long_sig <- meta_long %>% filter(FDR_meta_stouffer < 0.05)

# add monkey as an arm (consistent columns)
monkey_meta <- monkey_DE2 %>%
  transmute(
    gene,
    FDR_meta_stouffer = FDR,
    dataset = "PFOA_monkey")

# append monkey and re-filter significant rows
meta_long <- dplyr::bind_rows(meta_long, monkey_meta)
meta_long_sig <- meta_long %>% filter(FDR_meta_stouffer < 0.05)

# per-gene summary across arms
meta_summary <- meta_long_sig %>%
  group_by(gene) %>%
  summarise(
    n_meta_sig = n(),                      
    min_FDR    = min(FDR_meta_stouffer)) %>%
  arrange(min_FDR)

# load robust table for overlap (already written above)
robust <- read.csv("robust_significant_genes4.csv")

# intersection of robust DESeq2 with meta summary
overlap <- robust %>% inner_join(meta_summary, by = "gene")
head(overlap)  # quick sanity check

###########################
## Venn diagram: meta analysis vs indivdual analysis
meta_genes    <- unique(meta_long_sig$gene)
deseq2_genes <- unique(robust$gene)

n_meta    <- length(meta_genes)
n_deseq2  <- length(deseq2_genes)
n_overlap <- length(intersect(meta_genes, deseq2_genes))
jaccard   <- n_overlap / length(union(meta_genes, deseq2_genes))
message(sprintf("Meta: %d | DESeq2 robust: %d | Overlap: %d | Jaccard: %.3f",
                n_meta, n_deseq2, n_overlap, jaccard))

venn_fit <- eulerr::euler(list(
  `Meta (union)`   = meta_genes,
  `DESeq2 robust`  = deseq2_genes
))
plot(venn_fit,
     quantities = TRUE,
     labels     = list(col = "grey20", font = 2),
     fills      = list(fill = c("#4C97D8", "#D96B6B"), alpha = 0.6),
     edges      = list(lwd = 1.2, col = "grey25"))

###########################
## Meta consensus direction across arms (+ monkey)
meta_dir <- meta_long_sig %>%
  mutate(dir_arm = dplyr::case_when(
    !is.na(n_up) & !is.na(n_down) ~ sign(n_up - n_down),
    TRUE ~ NA_real_
  )) %>%
  dplyr::select(gene, dataset, dir_arm)

monkey_dir <- monkey_DE2 %>%
  transmute(gene, dataset = "PFOA_monkey", dir_arm = sign(signed_z))

meta_dir_all <- meta_dir %>%
  full_join(monkey_dir, by = c("gene","dataset"), suffix = c("", ".monk")) %>%
  transmute(gene, dataset, dir_arm = coalesce(dir_arm.monk, dir_arm))

dir_meta <- meta_dir_all %>%
  group_by(gene) %>%
  summarise(dir_meta = sign(sum(dir_arm, na.rm = TRUE)), .groups = "drop") %>%
  mutate(dir_plot = ifelse(dir_meta >= 0, "Up", "Down"))

## Top 50 genes (by min meta FDR) for barplot
top50 <- overlap %>%
  mutate(
    min_FDR      = as.numeric(min_FDR),
    min_FDR_safe = pmax(min_FDR, 1e-300),
    neglog10FDR  = -log10(min_FDR_safe)
  ) %>%
  arrange(min_FDR) %>%
  slice_head(n = 50) %>%
  left_join(dir_meta %>% dplyr::select(gene, dir_plot), by = "gene") %>%
  mutate(dir_plot = ifelse(is.na(dir_plot),
                           ifelse(mean_FC >= 0, "Up", "Down"),
                           dir_plot),
         gene = factor(gene, levels = rev(unique(gene))))

## Bar plot of top 50 (color by consensus direction; label # DESeq2 datasets)
max_x <- max(top50$neglog10FDR, na.rm = TRUE)
ggplot(top50, aes(x = neglog10FDR, y = gene, fill = dir_plot)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_datasets), hjust = -0.15, size = 3.1) +
  scale_fill_manual(values = c(Down = "royalblue", Up = "red3")) +
  coord_cartesian(xlim = c(0, max_x * 1.08)) +
  labs(x = expression(-log[10]("meta FDR (min across arms)")),
       y = NULL, fill = "Direction") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right")

############################
## significant pathways of meta analysis
# open singifciant meta genes
meta_genes <- unique(meta_long_sig$gene)

## Hallmark gene set
m_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

## Reactome gene set
m_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

## run ORA (no universe provided -> uses all symbols in TERM2GENE as background)
ek_h <- enricher(gene = meta_genes, TERM2GENE = m_h)
ek_rx <- enricher(gene = meta_genes, TERM2GENE = m_reactome)

## Create dot plots
# hallmark plot
dotplot(ek_h, showCategory = 20) +
  theme(axis.text.y = element_text(size = 8))
#reactome plot       
dotplot(ek_rx, showCategory = 20)+
  theme(axis.text.y = element_text(size = 8))
