# =============================================================================
# 06_survival_analysis.R
#
# Purpose: Identify genes whose expression independently predicts overall
#          survival in KIRC patients.
#
# Question answered:
#   Of the differentially expressed genes, which predict survival after
#   adjusting for clinical covariates (stage, grade, age, sex)?
#
# Method:
#   1. Subset to tumor samples with complete survival data
#   2. Univariate Cox regression on top DE genes
#   3. Kaplan-Meier curves for top hits
#   4. Multivariate Cox adjusting for stage, grade, age, sex
#   5. Build and evaluate a composite risk score
#
# Inputs:
#   - data/processed/kirc_vst.rds      (VST-transformed expression)
#   - data/processed/kirc_clinical.tsv (clinical/survival data)
#   - results/de_results_full.tsv      (DE results)
#
# Outputs:
#   - results/survival_univariate_cox.tsv
#   - results/survival_multivariate_cox.tsv
#   - figures/survival/ : KM curves, forest plot, risk score figures
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(survival)
library(survminer)
library(broom)
library(DESeq2)
library(SummarizedExperiment)


# ---- Step 1: Load inputs ----------------------------------------------------
vsd      <- readRDS(here("data", "processed", "kirc_vst.rds"))
clinical <- read.table(here("data", "processed", "kirc_clinical.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
de_df    <- read.table(here("results", "de_results_full.tsv"),
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

message("Loaded VST: ", nrow(vsd), " genes x ", ncol(vsd), " samples")
message("Loaded clinical: ", nrow(clinical), " rows")
message("Loaded DE results: ", nrow(de_df), " genes")


# ---- Step 2: Build the survival analysis cohort -----------------------------
# Keep tumor samples with:
#   - non-missing vital_status, os_time, os_event
#   - positive survival time (zeros are unusable)

surv_clin <- clinical %>%
  filter(sample_type == "Primary Tumor",
         !is.na(os_time_days),
         !is.na(os_event),
         os_time_days > 0) %>%
  mutate(
    # Recode stage to a simpler factor (collapse subgroups like "Stage IA" -> "Stage I")
    stage_clean = case_when(
      grepl("^Stage I[^IV]", stage)  ~ "Stage I",
      grepl("^Stage II[^I]", stage)  ~ "Stage II",
      grepl("^Stage III", stage)     ~ "Stage III",
      grepl("^Stage IV", stage)      ~ "Stage IV",
      TRUE                           ~ NA_character_
    ),
    stage_clean = factor(stage_clean,
                         levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    grade_clean = factor(grade,
                         levels = c("G1", "G2", "G3", "G4")),
    gender      = factor(gender)
  )

message("\n=== Cohort summary ===")
message("Tumor samples with survival: ", nrow(surv_clin))
message("Deaths (events): ", sum(surv_clin$os_event))
message("Median follow-up (months): ",
        round(median(surv_clin$os_time_months), 1))


# ---- Step 3: Align expression matrix with clinical cohort -------------------
# Subset VST to the cohort. Order matters: keep them aligned.

expr_mat <- assay(vsd)[, surv_clin$sample_id]
stopifnot(all(colnames(expr_mat) == surv_clin$sample_id))

message("\nExpression matrix aligned: ",
        nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")


# ---- Step 4: Pick candidate genes for survival testing ---------------------
# Testing all 31,000 genes is excessive; multiple testing burden destroys power.
# Better: test only the strongly DE genes (the "biologically interesting" set).
# We pick the top ~500 by combined effect + significance.

# Match DE results to genes in expression matrix (need version-stripped Ensembl)
de_df <- de_df %>%
  mutate(ensembl_versioned = gsub("\\..*", "", rownames(vsd)[
    match(ensembl, gsub("\\..*", "", rownames(vsd)))
  ]))

# Get top DE genes by absolute log2FC and significance
top_de <- de_df %>%
  filter(padj < 1e-10,            # very significant
         abs(log2FoldChange) > 2, # strong effect
         !is.na(symbol),
         symbol != "") %>%
  arrange(padj) %>%
  slice_head(n = 500)

message("\nCandidate genes for survival testing: ", nrow(top_de))


# ---- Step 5: Run univariate Cox regression for each gene -------------------
# For each gene: median-split tumor samples into "high" vs "low" expression,
# fit Cox model, extract HR / CI / p-value.
#
# Median-split is the convention in oncology biomarker work and is
# easier to interpret than continuous models. It loses statistical
# information but matches how clinical thresholds are used in practice.

message("\n=== Running univariate Cox regression on ", nrow(top_de),
        " candidate genes ===")

# Match the candidate genes to the expression matrix
match_idx <- match(top_de$ensembl,
                   gsub("\\..*", "", rownames(expr_mat)))

# Drop genes not in expression matrix (rare)
keep <- !is.na(match_idx)
top_de <- top_de[keep, ]
match_idx <- match_idx[keep]

# Function to fit Cox for one gene and return tidy result
fit_cox_one <- function(gene_idx, gene_info) {
  expr <- as.numeric(expr_mat[gene_idx, ])
  group <- factor(ifelse(expr > median(expr), "high", "low"),
                  levels = c("low", "high"))
  
  fit <- coxph(Surv(os_time_months, os_event) ~ group, data = surv_clin)
  s <- summary(fit)
  
  data.frame(
    ensembl   = gene_info$ensembl,
    symbol    = gene_info$symbol,
    log2FC_DE = gene_info$log2FoldChange,
    padj_DE   = gene_info$padj,
    HR        = s$conf.int[1, "exp(coef)"],
    HR_lower  = s$conf.int[1, "lower .95"],
    HR_upper  = s$conf.int[1, "upper .95"],
    cox_p     = s$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

uni_cox <- purrr::map_dfr(seq_along(match_idx), function(i) {
  tryCatch(
    fit_cox_one(match_idx[i], top_de[i, ]),
    error = function(e) NULL
  )
})

# Adjust p-values for multiple testing
uni_cox <- uni_cox %>%
  mutate(cox_padj = p.adjust(cox_p, method = "BH")) %>%
  arrange(cox_p)

n_sig <- sum(uni_cox$cox_padj < 0.05, na.rm = TRUE)
message("Genes with significant survival association (padj < 0.05): ", n_sig)

message("\n=== Top 15 prognostic genes (univariate Cox) ===")
print(uni_cox %>%
        slice_head(n = 15) %>%
        dplyr::select(symbol, log2FC_DE, HR, HR_lower, HR_upper, cox_p, cox_padj))

# Save
write.table(uni_cox, here("results", "survival_univariate_cox.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ---- Step 6: Pick the top genes for visualization --------------------------
# Top genes that are both DE and prognostic

top_prognostic <- uni_cox %>%
  filter(cox_padj < 0.05) %>%
  arrange(cox_p) %>%
  slice_head(n = 8)

message("\nTop 8 prognostic genes for visualization:")
print(top_prognostic$symbol)


# ---- Step 7: Kaplan-Meier curves -------------------------------------------
dir.create(here("figures", "survival"), recursive = TRUE, showWarnings = FALSE)

plot_km_one <- function(ens, sym) {
  idx <- which(gsub("\\..*", "", rownames(expr_mat)) == ens)[1]
  if (is.na(idx)) return(NULL)
  
  expr <- as.numeric(expr_mat[idx, ])
  group <- factor(ifelse(expr > median(expr), "High", "Low"),
                  levels = c("Low", "High"))
  df <- data.frame(time = surv_clin$os_time_months,
                   event = surv_clin$os_event,
                   group = group)
  
  fit <- survfit(Surv(time, event) ~ group, data = df)
  
  ggsurvplot(
    fit, data = df,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = c("#4C72B0", "#C44E52"),
    title = paste0("KIRC overall survival by ", sym, " expression"),
    xlab = "Months",
    legend.title = paste(sym, "expression"),
    legend.labs = c("Low (≤ median)", "High (> median)"),
    ggtheme = theme_minimal()
  )
}

# Save KM curve for each top prognostic gene
for (i in seq_len(nrow(top_prognostic))) {
  ens <- top_prognostic$ensembl[i]
  sym <- top_prognostic$symbol[i]
  km <- plot_km_one(ens, sym)
  if (!is.null(km)) {
    png(here("figures", "survival", paste0("km_", sym, ".png")),
        width = 900, height = 800, res = 130)
    print(km)
    dev.off()
  }
}

message("\nKaplan-Meier curves saved for ", nrow(top_prognostic), " genes")


# ---- Step 8: Multivariate Cox — does the gene predict independently? -------
# Test whether expression predicts survival AFTER adjusting for clinical vars.
# Only test the genes that were univariate-significant.

message("\n=== Running multivariate Cox (adjusting for stage, grade, age, sex) ===")

# Subset cohort to samples with complete clinical info
mv_cohort <- surv_clin %>%
  filter(!is.na(stage_clean), !is.na(grade_clean),
         !is.na(age), !is.na(gender))

message("Samples with complete clinical data for MV Cox: ", nrow(mv_cohort))

mv_expr <- expr_mat[, mv_cohort$sample_id]
sig_uni <- uni_cox %>% filter(cox_padj < 0.05) %>% arrange(cox_p)

fit_mv_one <- function(gene_idx, gene_info) {
  expr <- as.numeric(mv_expr[gene_idx, ])
  group <- factor(ifelse(expr > median(expr), "high", "low"),
                  levels = c("low", "high"))
  d <- cbind(mv_cohort, group = group)
  
  fit <- coxph(
    Surv(os_time_months, os_event) ~ group + stage_clean + grade_clean +
      age + gender,
    data = d
  )
  s <- summary(fit)
  gene_row <- which(rownames(s$coefficients) == "grouphigh")
  if (length(gene_row) == 0) return(NULL)
  
  data.frame(
    ensembl  = gene_info$ensembl,
    symbol   = gene_info$symbol,
    HR_mv    = s$conf.int[gene_row, "exp(coef)"],
    HR_lower = s$conf.int[gene_row, "lower .95"],
    HR_upper = s$conf.int[gene_row, "upper .95"],
    p_mv     = s$coefficients[gene_row, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

mv_indices <- match(sig_uni$ensembl,
                    gsub("\\..*", "", rownames(mv_expr)))

mv_cox <- purrr::map_dfr(seq_along(mv_indices), function(i) {
  if (is.na(mv_indices[i])) return(NULL)
  tryCatch(
    fit_mv_one(mv_indices[i], sig_uni[i, ]),
    error = function(e) NULL
  )
})

mv_cox <- mv_cox %>%
  mutate(padj_mv = p.adjust(p_mv, method = "BH")) %>%
  arrange(p_mv)

message("\nGenes independently prognostic (MV padj < 0.05): ",
        sum(mv_cox$padj_mv < 0.05, na.rm = TRUE))

message("\n=== Top 15 independently prognostic genes ===")
print(mv_cox %>%
        slice_head(n = 15) %>%
        dplyr::select(symbol, HR_mv, HR_lower, HR_upper, p_mv, padj_mv))

write.table(mv_cox, here("results", "survival_multivariate_cox.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ---- Step 9: Forest plot of top prognostic genes ---------------------------
forest_data <- mv_cox %>%
  filter(padj_mv < 0.05) %>%
  arrange(HR_mv) %>%
  slice_head(n = 20)

if (nrow(forest_data) > 0) {
  p_forest <- ggplot(forest_data,
                     aes(x = HR_mv, y = reorder(symbol, HR_mv))) +
    geom_point(size = 3, color = "#2E5BBA") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                   height = 0.25, color = "#2E5BBA") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    scale_x_log10() +
    labs(
      title = "Multivariate Cox — independent prognostic genes",
      subtitle = "Adjusted for stage, grade, age, sex",
      x = "Hazard Ratio (log scale)",
      y = NULL,
      caption = "HR > 1: high expression → worse outcome"
    ) +
    theme_minimal(base_size = 11)
  
  ggsave(here("figures", "survival", "forest_multivariate.png"),
         p_forest, width = 8, height = 6, dpi = 150)
}


# ---- Step 10: Build a composite risk score ---------------------------------
# Sum z-scored expression of the top independent prognostic genes, weighted
# by sign of effect (positive HR -> add, negative HR -> subtract).
# Then split patients into high/low risk and look at survival difference.

risk_genes <- mv_cox %>%
  filter(padj_mv < 0.05) %>%
  arrange(p_mv) %>%
  slice_head(n = 10)

if (nrow(risk_genes) >= 3) {
  ridx <- match(risk_genes$ensembl,
                gsub("\\..*", "", rownames(expr_mat)))
  ridx <- ridx[!is.na(ridx)]
  
  # Z-score each gene across samples
  rexpr <- t(scale(t(expr_mat[ridx, , drop = FALSE])))
  
  # Sign by direction of HR (HR > 1 means high expr = bad -> add as-is)
  signs <- sign(log(risk_genes$HR_mv))
  signs <- signs[seq_len(nrow(rexpr))]
  
  risk_score <- as.numeric(t(signs) %*% rexpr)
  
  risk_df <- data.frame(
    sample_id = colnames(expr_mat),
    risk_score = risk_score
  ) %>%
    inner_join(surv_clin %>% dplyr::select(sample_id, os_time_months, os_event),
               by = "sample_id") %>%
    mutate(risk_group = factor(ifelse(risk_score > median(risk_score),
                                      "High risk", "Low risk"),
                               levels = c("Low risk", "High risk")))
  
  # KM curve for risk score
  fit <- survfit(Surv(os_time_months, os_event) ~ risk_group, data = risk_df)
  
  png(here("figures", "survival", "risk_score_km.png"),
      width = 900, height = 800, res = 130)
  print(ggsurvplot(
    fit, data = risk_df,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = c("#4C72B0", "#C44E52"),
    title = "KIRC overall survival by composite transcriptomic risk score",
    subtitle = paste0("Built from top ", nrow(risk_genes),
                      " independently prognostic genes"),
    xlab = "Months",
    ggtheme = theme_minimal()
  ))
  dev.off()
  
  # Cox HR for the risk score itself
  rfit <- coxph(Surv(os_time_months, os_event) ~ risk_group, data = risk_df)
  message("\n=== Risk score HR ===")
  print(summary(rfit)$conf.int)
}

message("\nDone. Run 07_visualizations.R next.")