# =============================================================================
# 01_load_data.R
#
# Purpose: Load TCGA-KIRC RNA-seq data from recount3 and save locally.
#
# Inputs:
#   - None (downloads from recount3 server)
#
# Outputs:
#   - data/interim/kirc_rse.rds : RangedSummarizedExperiment with counts + metadata
#
# What this script does:
#   1. Connect to recount3, list available projects
#   2. Subset to the TCGA-KIRC project
#   3. Download counts and metadata as a RangedSummarizedExperiment (RSE) object
#   4. Inspect the data (sample count, structure)
#   5. Save the RSE locally so we never re-download
# =============================================================================


# ---- Load libraries ---------------------------------------------------------
# recount3: the package that pulls uniformly-processed RNA-seq data
# SummarizedExperiment: the data structure that holds counts + metadata
# here: makes file paths work regardless of where the script is run from
# dplyr: for any data manipulation

library(recount3)
library(SummarizedExperiment)
library(here)
library(dplyr)


# ---- Step 1: Get the list of available projects -----------------------------
# This queries the recount3 metadata server and returns a tibble of every
# project available, across SRA, GTEx, and TCGA.

message("Fetching list of available projects from recount3...")
human_projects <- available_projects(organism = "human")

# Quick peek
message("Total projects available: ", nrow(human_projects))
message("Project sources: ",
        paste(unique(human_projects$file_source), collapse = ", "))


# ---- Step 2: Find TCGA-KIRC -------------------------------------------------
# Filter to just TCGA-KIRC. There should be exactly one row.

kirc_info <- subset(
  human_projects,
  project == "KIRC" & file_source == "tcga"
)

# Sanity check
stopifnot(nrow(kirc_info) == 1)
message("Found TCGA-KIRC: ", kirc_info$n_samples, " samples expected")


# ---- Step 3: Download counts and metadata -----------------------------------
# create_rse() does the actual download. It creates a RangedSummarizedExperiment
# (RSE) — an R object that holds three things in one place:
#   1. The counts matrix (genes x samples)
#   2. The sample metadata (sample type, etc.) — accessed via colData()
#   3. The gene annotation — accessed via rowRanges()
#
# This downloads several files (~hundreds of MB), so first run takes a few
# minutes. After that it's cached.

message("Downloading TCGA-KIRC counts and metadata (this may take a few minutes)...")
rse_kirc <- create_rse(kirc_info)


# ---- Step 4: Inspect the data -----------------------------------------------
# Confirm we got what we expected.

message("\n=== RSE summary ===")
print(rse_kirc)

message("\n=== Number of genes: ", nrow(rse_kirc))
message("=== Number of samples: ", ncol(rse_kirc))
message("\n=== Available assays: ",
        paste(assayNames(rse_kirc), collapse = ", "))

# Look at sample type breakdown
# (Primary Tumor, Solid Tissue Normal, etc.)
message("\n=== Sample type breakdown:")
print(table(rse_kirc$tcga.cgc_sample_sample_type, useNA = "ifany"))


# ---- Step 5: Save locally ---------------------------------------------------
# Save the RSE to disk so we don't have to re-download next time.
# Using .rds (R's native binary format — small and fast).
# We save to data/interim/ because the data has been "fetched and assembled"
# but not yet cleaned or filtered for analysis.

output_path <- here("data", "interim", "kirc_rse.rds")
message("\nSaving RSE to: ", output_path)
saveRDS(rse_kirc, output_path)

message("\nDone. Run 02_qc_eda.R next.")