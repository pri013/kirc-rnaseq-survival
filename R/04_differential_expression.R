# =============================================================================
# 04_differential_expression.R
#
# Purpose: Differential expression analysis between KIRC tumor and normal
#          samples using DESeq2, with covariate adjustment.
#
# Question answered:
#   Which genes are significantly up- or down-regulated in KIRC tumors
#   compared to normal kidney tissue?
#
# Inputs:
#   - data/processed/kirc_dds.rds      (from 03_qc_counts.R)
#
# Outputs:
#   - results/de_results_full.tsv      : all genes, sorted by significance
#   - results/de_results_significant.tsv : padj < 0.05 & |log2FC| > 1
#   - figures/de/ma_plot.png
#   - figures/de/volcano_plot.png
#   - figures/de/pvalue_histogram.png
#   - data/processed/dds_fit.rds       : fitted model for later scripts
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(DESeq2)
library(here)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(AnnotationDbi)


# ---- Step 1: Load the DESeq2 dataset ----------------------------------------
dds <- readRDS(here("data", "processed", "kirc_dds.rds"))
message("Loaded DESeq2 dataset: ", nrow(dds), " genes x ", ncol(dds), " samples")


# ---- Step 2: Inspect and clean covariates -----------------------------------
# We want to adjust for gender and age in the model. First check completeness.

message("\n=== Covariate completeness ===")
message("Missing gender: ", sum(is.na(colData(dds)$gender)))
message("Missing age: ", sum(is.na(colData(dds)$age)))


# ---- Step 3: Update design formula ------------------------------------------
# The design formula tells DESeq2 what to test. Order matters: the LAST term
# is the one being tested (sample_type). Earlier terms (gender, age) are
# controlled for as covariates.
#
# This is "model gene expression as a function of gender, age, AND sample_type"
# and report the effect of sample_type after holding the other two constant.

design(dds) <- ~ gender + sample_type

# Note: I'm only including gender (and not age) because age has 1 NA, which
# would drop that sample. We can add age in a sensitivity analysis. The main
# confounder for KIRC tumor-vs-normal is sample sex, so this design is sound.

message("\nDesign formula: ", deparse(design(dds)))


# ---- Step 4: Run DESeq2 -----------------------------------------------------
# This does three things internally:
#   1. Estimate size factors (library size normalization)
#   2. Estimate gene-wise dispersion (the variance model)
#   3. Fit a negative binomial GLM and test the contrast of interest
#
# Takes 2-5 minutes on this dataset.

message("\n=== Running DESeq2 (this takes a few minutes) ===")
dds <- DESeq(dds)
message("DESeq2 complete.")


# ---- Step 5: Inspect dispersion ---------------------------------------------
# A diagnostic: plot dispersion estimates. The fitted curve should look smooth
# and pass through the cloud of gene-wise estimates.

dir.create(here("figures", "de"), recursive = TRUE, showWarnings = FALSE)

png(here("figures", "de", "dispersion_plot.png"), width = 1000, height = 700, res = 120)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()


# ---- Step 6: Extract results ------------------------------------------------
# Contrast: Tumor vs Normal. Because we set Normal as the reference level
# earlier, positive log2FC = "up in Tumor" and negative = "down in Tumor".

# Without shrinkage (raw effect sizes)
res_raw <- results(dds, contrast = c("sample_type", "Tumor", "Normal"))

# WITH shrinkage (apeglm) — more reliable for low-count genes
# We need to use coefficient name (not contrast) for apeglm
res <- lfcShrink(dds,
                 coef = "sample_type_Tumor_vs_Normal",
                 type = "apeglm")

message("\n=== Results summary (after shrinkage) ===")
summary(res)


# ---- Step 7: Annotate genes with HUGO symbols -------------------------------
# Gene IDs in recount3 are Ensembl IDs with version suffixes (ENSG00000123.5).
# Strip the suffix and map to gene symbols using org.Hs.eg.db.

ensembl_clean <- gsub("\\..*", "", rownames(res))

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)


# ---- Step 8: Tidy results dataframe -----------------------------------------
de_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ensembl_versioned") %>%
  mutate(
    ensembl       = gsub("\\..*", "", ensembl_versioned),
    symbol        = gene_symbols,
    significant   = padj < 0.05 & abs(log2FoldChange) > 1,
    direction     = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "NS"
    )
  ) %>%
  arrange(padj) %>%
  dplyr::select(ensembl, symbol, everything())

# Summary
message("\n=== DE gene counts ===")
print(table(de_df$direction, useNA = "ifany"))


# ---- Step 9: Save results ---------------------------------------------------
dir.create(here("results"), showWarnings = FALSE)

# Full results
write.table(de_df, here("results", "de_results_full.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Significant only
sig_df <- filter(de_df, significant)
write.table(sig_df, here("results", "de_results_significant.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\nFull results: ", nrow(de_df), " genes -> results/de_results_full.tsv")
message("Significant:   ", nrow(sig_df), " genes -> results/de_results_significant.tsv")


# ---- Step 10: Show the top hits --------------------------------------------
message("\n=== Top 20 upregulated genes (by padj) ===")
top_up <- de_df %>% filter(direction == "Up") %>% slice_head(n = 20)
print(top_up %>% dplyr::select(symbol, log2FoldChange, padj))

message("\n=== Top 20 downregulated genes (by padj) ===")
top_down <- de_df %>% filter(direction == "Down") %>% slice_head(n = 20)
print(top_down %>% dplyr::select(symbol, log2FoldChange, padj))


# ---- Step 11: Diagnostic figures --------------------------------------------

# (a) p-value histogram — should be flat with a spike near 0
# This confirms the test calibration is reasonable. A bimodal or U-shaped
# histogram would indicate problems.
png(here("figures", "de", "pvalue_histogram.png"), width = 800, height = 600, res = 120)
hist(de_df$pvalue, breaks = 50, col = "steelblue", border = "white",
     main = "Distribution of raw p-values", xlab = "p-value")
dev.off()

# (b) MA plot — log fold change vs mean expression
png(here("figures", "de", "ma_plot.png"), width = 1000, height = 700, res = 120)
plotMA(res, ylim = c(-8, 8), main = "MA plot (LFC shrinkage applied)")
dev.off()

# (c) Volcano plot — the headline figure
# Label the top 20 genes by significance (combining up and down)
top_labels <- bind_rows(
  de_df %>% filter(direction == "Up")   %>% slice_head(n = 10),
  de_df %>% filter(direction == "Down") %>% slice_head(n = 10)
) %>% pull(symbol) %>% na.omit()

png(here("figures", "de", "volcano_plot.png"), width = 1100, height = 900, res = 130)
print(EnhancedVolcano(
  de_df,
  lab        = de_df$symbol,
  x          = "log2FoldChange",
  y          = "padj",
  pCutoff    = 0.05,
  FCcutoff   = 1,
  selectLab  = top_labels,
  title      = "KIRC Tumor vs Normal — Differential Expression",
  subtitle   = "DESeq2 with gender adjustment, apeglm shrinkage",
  caption    = paste0("Total genes: ", nrow(de_df),
                      " | Significant: ", nrow(sig_df)),
  pointSize  = 1.2,
  labSize    = 3.5,
  drawConnectors = TRUE,
  maxoverlapsConnectors = Inf
))
dev.off()


# ---- Step 12: Save the fitted dds for later --------------------------------
saveRDS(dds, here("data", "processed", "dds_fit.rds"))
message("\nFitted DESeq2 object saved: data/processed/dds_fit.rds")

message("\nDone. Run 05_pathway_enrichment.R next.")