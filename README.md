# KIRC RNA-seq: Differential Expression and Survival Analysis

End-to-end RNA-seq analysis of clear cell renal cell carcinoma (KIRC) using TCGA data, built to demonstrate a full computational genomics workflow — from raw sequencing data through clinical interpretation.

## Research Question

Which genes and biological pathways are dysregulated in KIRC tumors compared to matched normal kidney tissue, and which of those dysregulated genes independently predict patient overall survival after adjusting for clinical covariates?

## Approach

Three layers, each answering one part of the question:

1. **Pipeline (Nextflow)** — containerized FASTQ-to-counts pipeline demonstrating reproducible large-scale data processing
2. **Differential expression and pathway analysis (R + Python)** — DESeq2 / PyDESeq2 with GSEA against MSigDB Hallmarks
3. **Survival analysis** — univariate and multivariate Cox regression to test whether top dysregulated genes predict outcomes independently of stage, grade, age, and sex

## Dataset

- **TCGA-KIRC** via the `recount3` resource (uniformly processed)
- ~530 tumor samples, ~70 matched normal samples
- Clinical and survival metadata from cBioPortal
- Subset of raw FASTQ files from GDC for pipeline demonstration

## Tech Stack

- **R:** DESeq2, recount3, clusterProfiler, fgsea, survival, survminer, shiny
- **Python:** PyDESeq2, gseapy, lifelines, streamlit, pandas
- **Pipeline:** Nextflow (DSL2), Docker, Salmon, FastQC, fastp, MultiQC
- **Reporting:** Quarto

## Status

🚧 In active development.

- [x] Project scaffold
- [ ] R environment and data loading
- [ ] QC and exploratory analysis
- [ ] Differential expression
- [ ] Pathway enrichment
- [ ] Survival analysis
- [ ] Nextflow pipeline
- [ ] Python replication
- [ ] Interactive dashboard
- [ ] Final report

## Repository Structure
kirc-rnaseq-survival/
├── R/                  # R analysis scripts (numbered in execution order)
├── data/               # Data (gitignored)
│   ├── raw/            # Original downloaded data
│   ├── interim/        # Partially processed
│   └── processed/      # Analysis-ready
├── notebooks/          # Exploratory analysis
├── reports/            # Quarto reports
├── figures/            # Generated figures
├── results/            # Output tables (DE results, etc.)
├── docs/               # Documentation
└── tests/              # Sanity-check scripts

## Reproducibility

Detailed reproduction instructions will be added as the project develops.

## License

MIT