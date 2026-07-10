# miRNA Processing, Differential Expression, and Network Analysis Pipeline for HR-MDS AZA Response

This repository contains a reproducible, end-to-end bioinformatics pipeline designed to process single-end small RNA-seq data and execute downstream differential expression, pathway enrichment, and network visualization. 

The workflow is tailored for analyzing bone marrow mononuclear cells (BMMCs) from High-Risk Myelodysplastic Syndrome (HR-MDS) patients collected pre- and post-Azacitidine (AZA) treatment to discover collective microRNA regulatory dynamics associated with clinical response.

---

## Pipeline Architecture Overview

The pipeline is split into two logical phases:
1. **Upstream Preprocessing (BASH/Python):** Raw FASTQ processing, quality control, deterministic exact-match alignment against mature human miRNAs, and UMI deduplication.
2. **Downstream Statistical Analysis (R):** Interaction-design differential expression via `DESeq2`, threshold-free functional enrichment via `fgsea`, and custom network/driver expression plotting.

---

## Prerequisites & Dependencies

To execute this pipeline, ensure your system has a Linux/BASH environment with the following tools installed and accessible in your `$PATH`:

### Command Line Tools
* **GNU Parallel** (for multi-sample parallel processing)
* **FastQC** (v0.12.1)
* **UMI-tools** (v1.1.6)
* **Cutadapt** (v5.2)
* **BWA** (v0.7.19)
* **SAMtools** (v1.21)

### Python Environment
* **Python 3.x**
* `pandas`

### R Environment & Packages (v4.4+)
* `DESeq2` (v1.44.0)
* `tidyverse` (v2.0.0)
* `fgsea` (v1.34.2)
* `ggrepel`, `ggraph`, `tidygraph`, `patchwork`, `scales`, `writexl`, `yaml`

---

## Repository Structure

```text
MDS-miRNA-AZA-pipeline/
├── README.md                 # Pipeline documentation
├── config.yaml               # Modular workspace, input/output, and threshold configs
├── format_miRBase_FASTA.sh   # Utility script to format miRBase mature sequence headers
├── BWA_index_creation.sh     # Utility script to build reference genome index
├── miRNA_pipeline.sh         # Core automated multi-sample preprocessing shell pipeline
├── count_combination.py      # Python script compiling individual counts into an expression matrix
└── gsea_on_mirs.r            # Master R script executing DEA, GSEA, and network generation
