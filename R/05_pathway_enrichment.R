# =============================================================================
# 05_pathway_enrichment.R
#
# Purpose: Identify biological pathways and processes enriched among the
#          differentially expressed genes from 04_differential_expression.R.
#
# Question answered:
#   What biological pathways are dysregulated in KIRC tumors?
#
# Methods:
#   - GSEA (Gene Set Enrichment Analysis) using fgsea, against:
#       * MSigDB Hallmarks (50 curated hallmark pathways)
#       * KEGG pathways
#   - ORA (Over-Representation Analysis) of significant DE genes via
#     clusterProfiler against GO Biological Process
#
# Inputs:
#   - results/de_results_full.tsv   (from 04_differential_expression.R)
#
# Outputs:
#   - results/gsea_hallmarks.tsv
#   - results/gsea_kegg.tsv
#   - results/ora_go_bp.tsv
#   - figures/pathway/ : enrichment plots
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(here)
library(dplyr)
library(tibble)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


# ---- Step 1: Load DE results ------------------------------------------------
de_df <- read.table(
  here("results", "de_results_full.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

message("Loaded DE results: ", nrow(de_df), " genes")


# ---- Step 2: Build the ranked gene list for GSEA ----------------------------
# GSEA needs a single ranked vector: numeric scores named by gene.
# Best practice: rank by the DESeq2 'stat' column (signed test statistic) which
# combines fold change and significance — gives sharper rankings than log2FC.
# But after lfcShrink with apeglm we lost 'stat' (it isn't computed for
# shrunken estimates). So we use a proxy: log2FoldChange * -log10(pvalue).
#
# This is a well-accepted alternative and preserves direction + significance.

ranking <- de_df %>%
  filter(!is.na(symbol),
         !is.na(log2FoldChange),
         !is.na(pvalue),
         pvalue > 0,                 # avoid log10(0)
         symbol != "") %>%
  mutate(rank_score = log2FoldChange * -log10(pvalue)) %>%
  # If duplicate symbols exist, keep the most extreme score
  group_by(symbol) %>%
  slice_max(abs(rank_score), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(rank_score))

ranks <- setNames(ranking$rank_score, ranking$symbol)

message("\nRanked gene list: ", length(ranks), " genes")
message("Top of ranking (most up): ", paste(head(names(ranks), 5), collapse = ", "))
message("Bottom of ranking (most down): ", paste(tail(names(ranks), 5), collapse = ", "))


# ---- Step 3: Get Hallmark gene sets from MSigDB -----------------------------
# msigdbr makes MSigDB accessible from R. Hallmarks = category "H".

message("\n=== Loading MSigDB Hallmark gene sets ===")
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")

# Convert to a named list: name=pathway, value=character vector of gene symbols
hallmarks <- split(msig_h$gene_symbol, msig_h$gs_name)

message("Hallmark sets loaded: ", length(hallmarks))


# ---- Step 4: Run GSEA against Hallmarks -------------------------------------
message("\n=== Running GSEA against Hallmarks ===")
set.seed(42)
gsea_h <- fgsea(
  pathways = hallmarks,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500,
  nproc    = 1   # 1 to keep reproducible
)

# Sort by significance
gsea_h <- gsea_h %>%
  as_tibble() %>%
  arrange(padj)

message("Hallmark pathways tested: ", nrow(gsea_h))
message("Significant (padj < 0.05): ", sum(gsea_h$padj < 0.05))

message("\n=== Top 10 enriched Hallmark pathways ===")
print(gsea_h %>%
        slice_head(n = 10) %>%
        dplyr::select(pathway, NES, padj, size))


# ---- Step 5: Run GSEA against KEGG ------------------------------------------
message("\n=== Loading KEGG gene sets ===")
msig_kegg <- msigdbr(species = "Homo sapiens", category = "C2",
                     subcollection = "CP:KEGG_LEGACY")

kegg_list <- split(msig_kegg$gene_symbol, msig_kegg$gs_name)
message("KEGG pathways loaded: ", length(kegg_list))

message("\n=== Running GSEA against KEGG ===")
set.seed(42)
gsea_kegg <- fgsea(
  pathways = kegg_list,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500,
  nproc    = 1
) %>% as_tibble() %>% arrange(padj)

message("Top 10 enriched KEGG pathways:")
print(gsea_kegg %>%
        slice_head(n = 10) %>%
        dplyr::select(pathway, NES, padj, size))


# ---- Step 6: ORA against GO Biological Process ------------------------------
# This complements GSEA by asking "are any GO BP terms over-represented among
# the significant DE genes?"

message("\n=== Running ORA (GO Biological Process) ===")

# Map significant gene symbols to Entrez IDs (clusterProfiler prefers these)
sig_symbols <- de_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1, !is.na(symbol)) %>%
  pull(symbol)

sig_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = sig_symbols,
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID)

# Background universe (Entrez IDs for all tested genes)
universe_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de_df$symbol[!is.na(de_df$symbol)]),
  columns = "ENTREZID",
  keytype = "SYMBOL"
) %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID)

ora_go <- enrichGO(
  gene          = sig_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

message("Significant GO BP terms: ", nrow(ora_go))

ora_go_df <- as.data.frame(ora_go) %>%
  arrange(p.adjust) %>%
  as_tibble()

message("\n=== Top 15 enriched GO BP terms ===")
print(ora_go_df %>%
        slice_head(n = 15) %>%
        dplyr::select(Description, Count, p.adjust))


# ---- Step 7: Save results ---------------------------------------------------
# fgsea results have a leading-edge gene column that is a list — collapse to
# strings for writing to TSV.

write_gsea_tsv <- function(df, path) {
  df %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    write.table(path, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_gsea_tsv(gsea_h, here("results", "gsea_hallmarks.tsv"))
write_gsea_tsv(gsea_kegg, here("results", "gsea_kegg.tsv"))

write.table(
  ora_go_df,
  here("results", "ora_go_bp.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("\nResults tables saved to results/")


# ---- Step 8: Figures --------------------------------------------------------
dir.create(here("figures", "pathway"), recursive = TRUE, showWarnings = FALSE)

# (a) Hallmark dot plot — the headline figure
hallmark_top <- gsea_h %>%
  filter(padj < 0.05) %>%
  slice_head(n = 25) %>%
  arrange(NES) %>%
  mutate(
    pathway_clean = gsub("HALLMARK_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    direction = ifelse(NES > 0, "Up in Tumor", "Down in Tumor")
  )

p_hallmark <- ggplot(hallmark_top,
                     aes(x = NES,
                         y = reorder(pathway_clean, NES),
                         size = -log10(padj),
                         color = direction)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c("Up in Tumor" = "#C44E52",
                                "Down in Tumor" = "#4C72B0")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "MSigDB Hallmarks — GSEA enrichment in KIRC Tumor vs Normal",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    size = expression(-log[10] ~ padj),
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave(here("figures", "pathway", "hallmark_dotplot.png"),
       p_hallmark, width = 10, height = 7, dpi = 150)


# (b) Top enrichment plot for the SINGLE top pathway (usually HYPOXIA)
top_pw <- gsea_h$pathway[1]
png(here("figures", "pathway", paste0("gsea_top_", top_pw, ".png")),
    width = 1000, height = 600, res = 130)
print(plotEnrichment(hallmarks[[top_pw]], ranks) +
        labs(title = top_pw,
             subtitle = paste0("NES = ", round(gsea_h$NES[1], 2),
                               ", padj = ", signif(gsea_h$padj[1], 3))))
dev.off()


# (c) Top enrichment plot for the most depleted pathway
gsea_h_down <- gsea_h %>% filter(NES < 0) %>% arrange(padj)
if (nrow(gsea_h_down) > 0) {
  top_down_pw <- gsea_h_down$pathway[1]
  png(here("figures", "pathway", paste0("gsea_top_down_", top_down_pw, ".png")),
      width = 1000, height = 600, res = 130)
  print(plotEnrichment(hallmarks[[top_down_pw]], ranks) +
          labs(title = top_down_pw,
               subtitle = paste0("NES = ", round(gsea_h_down$NES[1], 2),
                                 ", padj = ", signif(gsea_h_down$padj[1], 3))))
  dev.off()
}


# (d) KEGG dot plot
kegg_top <- gsea_kegg %>%
  filter(padj < 0.05) %>%
  slice_head(n = 20) %>%
  arrange(NES) %>%
  mutate(
    pathway_clean = gsub("KEGG_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    direction = ifelse(NES > 0, "Up in Tumor", "Down in Tumor")
  )

if (nrow(kegg_top) > 0) {
  p_kegg <- ggplot(kegg_top,
                   aes(x = NES,
                       y = reorder(pathway_clean, NES),
                       size = -log10(padj),
                       color = direction)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Up in Tumor" = "#C44E52",
                                  "Down in Tumor" = "#4C72B0")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "KEGG pathways — GSEA enrichment in KIRC Tumor vs Normal",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = expression(-log[10] ~ padj),
      color = NULL
    ) +
    theme_minimal(base_size = 11)
  
  ggsave(here("figures", "pathway", "kegg_dotplot.png"),
         p_kegg, width = 10, height = 7, dpi = 150)
}

message("\nFigures saved to figures/pathway/")
message("\nDone. Run 06_survival_analysis.R next.")