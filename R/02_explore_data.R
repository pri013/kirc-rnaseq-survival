# =============================================================================
# 02_explore_data.R
#
# Purpose: Inspect the TCGA-KIRC RSE to understand sample metadata, sample
#          types, and clinical variables before any analysis.
#
# Inputs:
#   - data/interim/kirc_rse.rds (from 01_load_data.R)
#
# Outputs:
#   - data/processed/kirc_clinical.tsv : cleaned sample-level clinical metadata
#   - figures/eda/                     : exploratory figures
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(SummarizedExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)


# ---- Step 1: Load the RSE ---------------------------------------------------
rse_path <- here("data", "interim", "kirc_rse.rds")
stopifnot(file.exists(rse_path))

message("Loading RSE from: ", rse_path)
rse_kirc <- readRDS(rse_path)

message("Loaded: ", nrow(rse_kirc), " genes x ", ncol(rse_kirc), " samples")


# ---- Step 2: Overall structure ----------------------------------------------
message("\n=== Class: ", class(rse_kirc))
message("\n=== Assays: ", paste(assayNames(rse_kirc), collapse = ", "))
message("\n=== Gene ID format (first 5 genes):")
print(head(rownames(rse_kirc), 5))


# ---- Step 3: Sample type breakdown ------------------------------------------
message("\n=== Sample type (cgc_sample_sample_type) ===")
print(table(rse_kirc$tcga.cgc_sample_sample_type, useNA = "ifany"))

message("\n=== gdc_cases.samples.sample_type ===")
print(table(rse_kirc$tcga.gdc_cases.samples.sample_type, useNA = "ifany"))


# ---- Step 4: Explore colData (sample metadata) ------------------------------
coldata <- as.data.frame(colData(rse_kirc))
message("\n=== colData dimensions: ", nrow(coldata), " samples x ",
        ncol(coldata), " metadata columns")

clinical_cols <- grep(
  "age|gender|sex|stage|grade|race|days_to|vital|tumor|sample_type",
  colnames(coldata),
  ignore.case = TRUE,
  value = TRUE
)

message("\n=== Clinically relevant columns (", length(clinical_cols), " found):")
print(clinical_cols)

# ---- Step 5: Map clinical columns to standard names -------------------------
# We use the GDC-harmonized columns where available — these are the current
# source. Names confirmed from inspection of colData.

cols_needed <- list(
  sample_type           = "tcga.cgc_sample_sample_type",
  vital_status          = "tcga.gdc_cases.diagnoses.vital_status",
  days_to_death         = "tcga.gdc_cases.diagnoses.days_to_death",
  days_to_last_followup = "tcga.gdc_cases.diagnoses.days_to_last_follow_up",
  age                   = "tcga.xml_age_at_initial_pathologic_diagnosis",   # CHANGED
  gender                = "tcga.gdc_cases.demographic.gender",
  stage                 = "tcga.cgc_case_pathologic_stage",
  grade                 = "tcga.xml_neoplasm_histologic_grade",
  race                  = "tcga.gdc_cases.demographic.race"
)

# Sanity check — make sure every column we listed actually exists
missing_cols <- setdiff(unlist(cols_needed), colnames(coldata))
if (length(missing_cols) > 0) {
  stop("Missing columns in coldata: ",
       paste(missing_cols, collapse = ", "))
}

message("\n=== Mapped columns ===")
print(cols_needed)

# ---- Step 6: Build a tidy clinical table ------------------------------------
clinical_df <- data.frame(
  sample_id          = colnames(rse_kirc),
  sample_type        = coldata[[cols_needed$sample_type]],
  vital_status       = coldata[[cols_needed$vital_status]],
  days_to_death      = suppressWarnings(as.numeric(coldata[[cols_needed$days_to_death]])),
  days_to_followup   = suppressWarnings(as.numeric(coldata[[cols_needed$days_to_last_followup]])),
  age                = suppressWarnings(as.numeric(coldata[[cols_needed$age]])),
  gender             = coldata[[cols_needed$gender]],
  stage              = coldata[[cols_needed$stage]],
  grade              = coldata[[cols_needed$grade]],
  race               = coldata[[cols_needed$race]],
  stringsAsFactors   = FALSE
)


# ---- Step 7: Compute survival time + event ----------------------------------
clinical_df <- clinical_df %>%
  mutate(
    os_time_days = ifelse(tolower(vital_status) == "dead",
                          days_to_death,
                          days_to_followup),
    os_time_months = os_time_days / 30.44,
    os_event = case_when(
      tolower(vital_status) == "dead" ~ 1L,
      tolower(vital_status) == "alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Inspect vital_status values
message("\n=== Vital status values ===")
print(table(clinical_df$vital_status, useNA = "ifany"))


# ---- Step 8: Completeness check ---------------------------------------------
message("\n=== Missing values per column ===")
print(sapply(clinical_df, function(x) sum(is.na(x))))

tumor_survival <- clinical_df %>%
  filter(sample_type == "Primary Tumor",
         !is.na(os_time_days),
         !is.na(os_event),
         os_time_days > 0)

message("\n=== Tumor samples with complete survival data: ", nrow(tumor_survival))
message("=== Deaths: ", sum(tumor_survival$os_event == 1))
message("=== Median follow-up (months): ",
        round(median(tumor_survival$os_time_months, na.rm = TRUE), 1))


# ---- Step 9: Visualize sample composition -----------------------------------
dir.create(here("figures", "eda"), recursive = TRUE, showWarnings = FALSE)

p_sample_type <- clinical_df %>%
  count(sample_type) %>%
  ggplot(aes(x = reorder(sample_type, n), y = n, fill = sample_type)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.1) +
  coord_flip() +
  labs(title = "Sample type distribution", x = NULL, y = "Count") +
  theme_minimal()

ggsave(here("figures", "eda", "sample_type_distribution.png"),
       p_sample_type, width = 7, height = 4, dpi = 150)

p_stage <- clinical_df %>%
  filter(sample_type == "Primary Tumor", !is.na(stage)) %>%
  count(stage) %>%
  ggplot(aes(x = stage, y = n, fill = stage)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(title = "Tumor stage (Primary Tumor samples)", x = "Stage", y = "Count") +
  theme_minimal()

ggsave(here("figures", "eda", "stage_distribution.png"),
       p_stage, width = 7, height = 4, dpi = 150)

p_age <- clinical_df %>%
  filter(!is.na(age)) %>%
  ggplot(aes(x = age, fill = sample_type)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  labs(title = "Age at diagnosis by sample type",
       x = "Age (years)", y = "Count", fill = NULL) +
  theme_minimal()

ggsave(here("figures", "eda", "age_distribution.png"),
       p_age, width = 7, height = 4, dpi = 150)

message("\nFigures saved to figures/eda/")


# ---- Step 10: Save clinical table -------------------------------------------
output_path <- here("data", "processed", "kirc_clinical.tsv")
write.table(clinical_df, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

message("\nClinical table saved to: ", output_path)
message("Rows: ", nrow(clinical_df), ", Cols: ", ncol(clinical_df))

message("\nDone.")