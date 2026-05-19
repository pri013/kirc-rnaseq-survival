# =============================================================================
# 03_qc_counts.R
#
# Purpose: Quality-control the RNA-seq counts before differential expression.
#          Filter low-count genes, transform counts, run PCA, identify any
#          outlier samples.
#
# Inputs:
#   - data/interim/kirc_rse.rds          (from 01_load_data.R)
#   - data/processed/kirc_clinical.tsv   (from 02_explore_data.R)
#
# Outputs:
#   - data/processed/kirc_dds.rds        : DESeq2 dataset (filtered, ready)
#   - data/processed/kirc_vst.rds        : VST-transformed counts
#   - figures/qc/*.png                   : QC figures
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(SummarizedExperiment)
library(DESeq2)
library(recount3)       # for transform_counts()
library(here)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# ---- Step 1: Load inputs ----------------------------------------------------
rse_kirc <- readRDS(here("data", "interim", "kirc_rse.rds"))
clinical <- read.table(here("data", "processed", "kirc_clinical.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

message("Loaded RSE: ", nrow(rse_kirc), " genes x ", ncol(rse_kirc), " samples")
message("Loaded clinical: ", nrow(clinical), " rows")

# ---- Step 2: Scale raw counts -----------------------------------------------
# recount3 stores counts as "raw_counts" but they're actually base-pair coverage
# summed across the gene. We need to convert to read-equivalent counts using
# transform_counts(). This is recount3's standard procedure.
# DESeq2 requires the 'counts' assay to be first in the assays list.

message("\n=== Scaling raw counts to read-equivalents ===")
scaled_counts <- transform_counts(rse_kirc)
assays(rse_kirc) <- list(counts = scaled_counts, raw_counts = assay(rse_kirc, "raw_counts"))
message("Assays now: ", paste(assayNames(rse_kirc), collapse = ", "))

# ---- Step 3: Attach clinical metadata to RSE colData ------------------------
# We want sample_type, vital_status, age, gender, stage, grade easily accessible
# alongside the counts. Match by sample_id (which is colnames(rse_kirc)).

stopifnot(all(colnames(rse_kirc) == clinical$sample_id))

# Add a clean sample_type factor with Normal as the reference level
# (so log2FC in DESeq2 is Tumor vs Normal, the direction we want)
colData(rse_kirc)$sample_type <- factor(
  ifelse(clinical$sample_type == "Primary Tumor", "Tumor",
         ifelse(clinical$sample_type == "Solid Tissue Normal", "Normal", NA)),
  levels = c("Normal", "Tumor")  # Normal first = reference
)

colData(rse_kirc)$gender <- factor(clinical$gender)
colData(rse_kirc)$age    <- clinical$age

message("\n=== Sample type distribution in RSE ===")
print(table(colData(rse_kirc)$sample_type, useNA = "ifany"))


# ---- Step 4: Drop ambiguous samples -----------------------------------------
# Keep only Tumor and Normal samples for the DE analysis.

keep_samples <- !is.na(colData(rse_kirc)$sample_type)
message("\nDropping ", sum(!keep_samples), " ambiguous samples")
rse_kirc <- rse_kirc[, keep_samples]
message("Kept: ", ncol(rse_kirc), " samples")


# ---- Step 5: Library size check ---------------------------------------------
# Library size = total reads per sample. Wildly different sizes can indicate
# technical issues. We expect ~30-80 million reads for TCGA RNA-seq.

lib_sizes <- colSums(assay(rse_kirc, "counts"))
message("\n=== Library size summary (millions) ===")
print(summary(lib_sizes / 1e6))

# Save a library-size figure
dir.create(here("figures", "qc"), recursive = TRUE, showWarnings = FALSE)

lib_df <- data.frame(
  sample = colnames(rse_kirc),
  lib_size_M = lib_sizes / 1e6,
  sample_type = colData(rse_kirc)$sample_type
)

p_lib <- ggplot(lib_df, aes(x = sample_type, y = lib_size_M, fill = sample_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  labs(title = "Library size by sample type",
       x = NULL, y = "Library size (millions of reads)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(here("figures", "qc", "library_size_boxplot.png"),
       p_lib, width = 6, height = 4, dpi = 150)


# ---- Step 6: Build DESeq2 dataset -------------------------------------------
# DESeqDataSet is the core data structure for DESeq2. The design formula tells
# it which variable we're testing — here, sample_type (Tumor vs Normal).
# We'll add covariates (age, gender) in the DE script.

message("\n=== Building DESeqDataSet ===")

# DESeq2 needs integer counts
mode(assay(rse_kirc, "counts")) <- "integer"

dds <- DESeqDataSet(rse_kirc, design = ~ sample_type)


# ---- Step 7: Filter low-count genes -----------------------------------------
# A gene is "informative" if it has reasonable expression in a reasonable number
# of samples. We require at least 10 counts in at least 10 samples — gentle
# enough to keep tissue-specific genes but strict enough to drop pure noise.

message("\n=== Filtering low-count genes ===")
message("Before filter: ", nrow(dds), " genes")

keep_genes <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep_genes, ]

message("After filter: ", nrow(dds), " genes")
message("Removed: ", sum(!keep_genes), " low-count genes")


# ---- Step 8: Variance-stabilizing transformation ----------------------------
# Raw counts are skewed (a few highly expressed genes dominate). For PCA and
# heatmaps we need a transformation that puts genes on comparable scales and
# stabilizes variance.
# vst() is fast and works well for samples with relatively similar library
# sizes (which is true for TCGA).
# blind = FALSE means we let the transformation know about the experimental
# design — recommended for downstream visualization of the contrast of
# interest.

message("\n=== Computing VST (variance-stabilizing transformation) ===")
message("This takes 1-3 minutes...")
vsd <- vst(dds, blind = FALSE)
message("Done. VST matrix: ", nrow(vsd), " genes x ", ncol(vsd), " samples")


# ---- Step 9: PCA -------------------------------------------------------------
# PCA on the top variable genes — we expect tumor and normal to separate clearly
# along PC1. If they don't, something is wrong.

message("\n=== Computing PCA ===")

# Use DESeq2's built-in plotPCA; returnData=TRUE gives us the underlying values
pca_data <- plotPCA(vsd, intgroup = "sample_type", returnData = TRUE, ntop = 500)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Normal" = "#4C72B0", "Tumor" = "#C44E52")) +
  labs(
    title = "PCA — top 500 most variable genes",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = NULL
  ) +
  theme_minimal()

ggsave(here("figures", "qc", "pca_sample_type.png"),
       p_pca, width = 7, height = 5, dpi = 150)

message("PC1 explains ", percent_var[1], "% of variance")
message("PC2 explains ", percent_var[2], "% of variance")


# ---- Step 10: Sample-to-sample distance heatmap -----------------------------
# An alternative view: how similar is each sample to every other sample?
# Tumors should cluster with tumors, normals with normals.

message("\n=== Computing sample distance heatmap ===")

# Use a random subset of 100 samples — full 600x600 heatmap is too cluttered
set.seed(42)
sub_idx <- sample(ncol(vsd), 100)
sample_dists <- dist(t(assay(vsd)[, sub_idx]))
sample_dist_matrix <- as.matrix(sample_dists)

# Annotation: color rows by sample type
annot <- data.frame(
  SampleType = colData(vsd)$sample_type[sub_idx],
  row.names = colnames(vsd)[sub_idx]
)

annot_colors <- list(SampleType = c(Normal = "#4C72B0", Tumor = "#C44E52"))

png(here("figures", "qc", "sample_distance_heatmap.png"),
    width = 1200, height = 1000, res = 120)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_row = annot,
  annotation_col = annot,
  annotation_colors = annot_colors,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Sample-to-sample distance (random subset of 100 samples)"
)
dev.off()

message("Heatmap saved")


# ---- Step 11: Save the DESeq2 dataset and VST ------------------------------
saveRDS(dds, here("data", "processed", "kirc_dds.rds"))
saveRDS(vsd, here("data", "processed", "kirc_vst.rds"))

message("\nDESeq2 dataset saved: data/processed/kirc_dds.rds")
message("VST saved: data/processed/kirc_vst.rds")

message("\nDone. Run 04_differential_expression.R next.")