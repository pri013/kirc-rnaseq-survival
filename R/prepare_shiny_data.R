# =============================================================================
# prepare_shiny_data.R
#
# Prepares a single .rds file containing all data the Shiny app needs:
#   - DE results (with symbols, log2FC, padj)
#   - VST expression matrix (genes x tumor samples)
#   - Clinical data for tumor samples
#   - Univariate Cox results
#   - Multivariate Cox results
#   - GSEA Hallmark results
#
# Output: data/processed/shiny_data.rds
# =============================================================================

library(here)
library(SummarizedExperiment)
library(dplyr)

# Load saved analysis outputs
vsd      <- readRDS(here("data", "processed", "kirc_vst.rds"))
clinical <- read.table(here("data", "processed", "kirc_clinical.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
de_df    <- read.table(here("results", "de_results_full.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
uni_cox  <- read.table(here("results", "survival_univariate_cox.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mv_cox   <- read.table(here("results", "survival_multivariate_cox.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gsea_h   <- read.table(here("results", "gsea_hallmarks.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Tumor samples with valid survival data
tumor_clinical <- clinical %>%
  filter(sample_type == "Primary Tumor",
         !is.na(os_time_days),
         !is.na(os_event),
         os_time_days > 0)

# Subset expression matrix to tumor samples only
expr_tumor <- assay(vsd)[, tumor_clinical$sample_id]

# Build a gene lookup table: symbol -> ensembl with version, log2FC, padj
gene_lookup <- de_df %>%
  filter(!is.na(symbol), symbol != "") %>%
  dplyr::select(symbol, ensembl, log2FoldChange, padj, direction)

# Map ensembl in expression matrix
expr_ensembl_versioned <- rownames(expr_tumor)
expr_ensembl_clean     <- gsub("\\..*", "", expr_ensembl_versioned)

# Save it all into one list
shiny_data <- list(
  de_df          = de_df,
  gene_lookup    = gene_lookup,
  expr_tumor     = expr_tumor,
  expr_ensembl_clean = expr_ensembl_clean,
  tumor_clinical = tumor_clinical,
  uni_cox        = uni_cox,
  mv_cox         = mv_cox,
  gsea_h         = gsea_h
)

# Save
output_path <- here("data", "processed", "shiny_data.rds")
saveRDS(shiny_data, output_path)

message("Saved shiny_data to: ", output_path)
message("File size: ",
        round(file.info(output_path)$size / 1024^2, 1), " MB")