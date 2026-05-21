# =============================================================================
# KIRC RNA-seq Dashboard
#
# Interactive exploration of differential expression and survival analysis
# results from TCGA-KIRC.
# =============================================================================

library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(survival)
library(survminer)
library(here)


# ---- Load data --------------------------------------------------------------
shiny_data <- readRDS(here("data", "processed", "shiny_data.rds"))

de_df          <- shiny_data$de_df
gene_lookup    <- shiny_data$gene_lookup
expr_tumor     <- shiny_data$expr_tumor
expr_ensembl   <- shiny_data$expr_ensembl_clean
tumor_clinical <- shiny_data$tumor_clinical
uni_cox        <- shiny_data$uni_cox
mv_cox         <- shiny_data$mv_cox
gsea_h         <- shiny_data$gsea_h


# ---- UI ---------------------------------------------------------------------
ui <- page_navbar(
  title = "KIRC RNA-seq Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  # ---- Tab 1: Overview ----
  nav_panel(
    title = "Overview",
    icon = icon("info-circle"),
    div(
      style = "max-width: 900px; margin: auto; padding: 20px;",
      h2("Clear Cell Renal Cell Carcinoma (KIRC) — TCGA Analysis"),
      p("Interactive exploration of differential expression and survival analysis from",
        strong("543 KIRC tumors"), "and", strong("72 matched normal kidney samples"),
        "in The Cancer Genome Atlas (TCGA-KIRC)."),
      hr(),
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(
          title = "Differentially Expressed Genes",
          value = format(sum(de_df$direction != "NS"), big.mark = ","),
          showcase = icon("dna"),
          theme = "primary"
        ),
        value_box(
          title = "Significant Pathways (Hallmarks)",
          value = sum(gsea_h$padj < 0.05),
          showcase = icon("project-diagram"),
          theme = "success"
        ),
        value_box(
          title = "Independent Prognostic Genes",
          value = sum(mv_cox$padj_mv < 0.05, na.rm = TRUE),
          showcase = icon("heart-pulse"),
          theme = "danger"
        )
      ),
      hr(),
      h4("Research Questions"),
      tags$ol(
        tags$li(strong("Differential expression:"),
                "Which genes are dysregulated in tumors vs normal?"),
        tags$li(strong("Pathway dysregulation:"),
                "Which biological pathways are coordinately altered?"),
        tags$li(strong("Clinical relevance:"),
                "Which genes independently predict patient survival?")
      ),
      hr(),
      h4("Key findings"),
      tags$ul(
        tags$li("Canonical KIRC signature confirmed: HIF-driven hypoxia (CA9, NDUFA4L2, EGLN3), Warburg metabolic shift, and immune infiltration"),
        tags$li("Composite 10-gene risk score: ", strong("HR = 2.72 (95% CI: 1.98–3.75, p < 0.0001)"),
                " for overall survival")
      ),
      hr(),
      p(em("Use the tabs above to explore individual genes, browse DE results, or view enriched pathways."))
    )
  ),
  
  # ---- Tab 2: Gene Explorer ----
  nav_panel(
    title = "Gene Explorer",
    icon = icon("search"),
    layout_sidebar(
      sidebar = sidebar(
        selectizeInput(
          "gene_input",
          "Select a gene (HUGO symbol):",
          choices = NULL,
          selected = "CA9",
          options = list(placeholder = "Type a gene symbol...")
        ),
        hr(),
        uiOutput("gene_summary")
      ),
      navset_card_tab(
        nav_panel("Expression",
                  plotlyOutput("expr_plot", height = "500px")),
        nav_panel("Survival (KM)",
                  plotOutput("km_plot", height = "500px")),
        nav_panel("DE statistics",
                  uiOutput("de_stats"))
      )
    )
  ),
  
  # ---- Tab 3: DE Browser ----
  nav_panel(
    title = "DE Browser",
    icon = icon("table"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        sliderInput("padj_thresh", "padj threshold (max):",
                    min = 0, max = 0.1, value = 0.05, step = 0.005),
        sliderInput("lfc_thresh", "|log2 fold change| ≥",
                    min = 0, max = 5, value = 1, step = 0.5),
        radioButtons("direction_filter", "Direction:",
                     choices = c("Both", "Up", "Down"),
                     selected = "Both"),
        hr(),
        textOutput("de_count")
      ),
      card(
        card_header("Differentially Expressed Genes"),
        DTOutput("de_table")
      )
    )
  ),
  
  # ---- Tab 4: Pathways ----
  nav_panel(
    title = "Pathways",
    icon = icon("diagram-project"),
    card(
      card_header("MSigDB Hallmark Pathway Enrichment (GSEA)"),
      plotlyOutput("pathway_plot", height = "600px")
    )
  )
)


# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  
  # Populate gene input dropdown server-side (faster for large lists)
  updateSelectizeInput(
    session, "gene_input",
    choices = sort(unique(gene_lookup$symbol)),
    selected = "CA9",
    server = TRUE
  )
  
  # Reactive: lookup selected gene's ensembl + expression
  selected_gene_data <- reactive({
    req(input$gene_input)
    gene_info <- gene_lookup %>% filter(symbol == input$gene_input) %>% slice(1)
    if (nrow(gene_info) == 0) return(NULL)
    
    idx <- which(expr_ensembl == gene_info$ensembl)[1]
    if (is.na(idx)) return(NULL)
    
    expr_vals <- as.numeric(expr_tumor[idx, ])
    
    list(
      info  = gene_info,
      expr  = expr_vals,
      idx   = idx
    )
  })
  
  # Gene summary in sidebar
  output$gene_summary <- renderUI({
    g <- selected_gene_data()
    req(g)
    tagList(
      h5(g$info$symbol),
      p("Ensembl: ", code(g$info$ensembl)),
      p("log2FoldChange: ", strong(round(g$info$log2FoldChange, 2))),
      p("padj: ", strong(format.pval(g$info$padj, digits = 3))),
      p("Direction: ", strong(g$info$direction))
    )
  })
  
  # Expression boxplot (tumor only — need normals to compare; let's load on demand)
  output$expr_plot <- renderPlotly({
    g <- selected_gene_data()
    req(g)
    df <- data.frame(
      sample = colnames(expr_tumor),
      expression = g$expr
    )
    p <- ggplot(df, aes(y = expression)) +
      geom_boxplot(fill = "#C44E52", alpha = 0.6, width = 0.4) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
      labs(
        title = paste0(g$info$symbol, " expression in KIRC tumors (n=",
                       ncol(expr_tumor), ")"),
        y = "Variance-stabilized expression",
        x = NULL
      ) +
      theme_minimal()
    ggplotly(p)
  })
  
  # KM curve
  output$km_plot <- renderPlot({
    g <- selected_gene_data()
    req(g)
    
    df <- data.frame(
      time = tumor_clinical$os_time_months,
      event = tumor_clinical$os_event,
      group = factor(
        ifelse(g$expr > median(g$expr), "High", "Low"),
        levels = c("Low", "High")
      )
    )
    
    fit <- survfit(Surv(time, event) ~ group, data = df)
    
    ggsurvplot(
      fit, data = df,
      pval = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      palette = c("#4C72B0", "#C44E52"),
      title = paste0("KIRC overall survival by ", g$info$symbol, " expression"),
      xlab = "Months",
      legend.title = paste(g$info$symbol, "expression"),
      legend.labs = c("Low (≤ median)", "High (> median)"),
      ggtheme = theme_minimal()
    )$plot
  })
  
  # DE statistics card
  output$de_stats <- renderUI({
    g <- selected_gene_data()
    req(g)
    
    uni_row <- uni_cox %>% filter(symbol == g$info$symbol)
    mv_row  <- mv_cox  %>% filter(symbol == g$info$symbol)
    
    tagList(
      h4("Differential expression"),
      p("log2FC: ", round(g$info$log2FoldChange, 3),
        " — padj: ", format.pval(g$info$padj, digits = 3)),
      hr(),
      h4("Univariate Cox"),
      if (nrow(uni_row) > 0) {
        tagList(
          p("HR: ", round(uni_row$HR[1], 2),
            " (95% CI: ", round(uni_row$HR_lower[1], 2),
            "–", round(uni_row$HR_upper[1], 2), ")"),
          p("p-value: ", format.pval(uni_row$cox_p[1], digits = 3))
        )
      } else p(em("Gene not in survival analysis candidate set.")),
      hr(),
      h4("Multivariate Cox (adjusted for stage, grade, age, sex)"),
      if (nrow(mv_row) > 0) {
        tagList(
          p("Adjusted HR: ", round(mv_row$HR_mv[1], 2),
            " (95% CI: ", round(mv_row$HR_lower[1], 2),
            "–", round(mv_row$HR_upper[1], 2), ")"),
          p("p-value: ", format.pval(mv_row$p_mv[1], digits = 3))
        )
      } else p(em("Not significant in univariate Cox, so not tested in multivariate."))
    )
  })
  
  # DE Browser
  filtered_de <- reactive({
    df <- de_df %>%
      filter(padj <= input$padj_thresh,
             abs(log2FoldChange) >= input$lfc_thresh,
             !is.na(symbol))
    if (input$direction_filter == "Up")   df <- df %>% filter(log2FoldChange > 0)
    if (input$direction_filter == "Down") df <- df %>% filter(log2FoldChange < 0)
    df
  })
  
  output$de_count <- renderText({
    paste0(nrow(filtered_de()), " genes match filters")
  })
  
  output$de_table <- renderDT({
    filtered_de() %>%
      dplyr::select(symbol, log2FoldChange, padj, direction) %>%
      mutate(
        log2FoldChange = round(log2FoldChange, 2),
        padj = format.pval(padj, digits = 3)
      ) %>%
      datatable(
        filter = "top",
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE
      )
  })
  
  # Pathway plot
  output$pathway_plot <- renderPlotly({
    pw <- gsea_h %>%
      filter(padj < 0.05) %>%
      slice_head(n = 25) %>%
      mutate(
        pathway_clean = gsub("HALLMARK_", "", pathway),
        pathway_clean = gsub("_", " ", pathway_clean),
        direction = ifelse(NES > 0, "Up in Tumor", "Down in Tumor"),
        neg_log10_padj = -log10(padj)
      )
    
    p <- ggplot(pw,
                aes(x = NES,
                    y = reorder(pathway_clean, NES),
                    size = neg_log10_padj,
                    color = direction,
                    text = paste0("Pathway: ", pathway_clean,
                                  "<br>NES: ", round(NES, 2),
                                  "<br>padj: ", signif(padj, 3),
                                  "<br>Genes: ", size))) +
      geom_point(alpha = 0.85) +
      scale_color_manual(values = c("Up in Tumor" = "#C44E52",
                                    "Down in Tumor" = "#4C72B0")) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      labs(
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        color = NULL,
        size = "-log10(padj)"
      ) +
      theme_minimal(base_size = 11)
    
    ggplotly(p, tooltip = "text")
  })
}


# ---- Run --------------------------------------------------------------------
shinyApp(ui, server)