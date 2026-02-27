# snRNA-seq Pipeline — DNAJC5/CSP Knockout in Mouse Neurons

[![Language: R](https://img.shields.io/badge/R-%E2%89%A54.5-276DC3.svg?logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5-2E9FDF.svg)](https://satijalab.org/seurat/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%E2%89%A53.18-85BB65.svg)](https://bioconductor.org/)
[![Reproducible: renv](https://img.shields.io/badge/Reproducible-renv-forestgreen.svg)](https://rstudio.github.io/renv/)
[![License: Academic](https://img.shields.io/badge/License-Academic%20use-blue.svg)](#license)

> Modular R pipeline for single-nucleus RNA-seq analysis of **DNAJC5/CSP knockout** effects in two distinct mouse neuronal populations, from 10X CellRanger output through clustering, differential expression, and multi-layered functional enrichment.

## Overview

End-to-end snRNA-seq pipeline investigating the transcriptomic consequences of **DNAJC5** (Cysteine String Protein / CSP) knockout in mouse neurons. Two transgenic lines are analyzed in parallel:

- **UBC-Cre × CSP-Flox** — pan-neuronal knockout (~30,000 nuclei, 4 samples)
- **PV-Cre × CSP-Flox** — GABAergic interneuron-specific knockout

DNAJC5 is a synaptic co-chaperone whose loss causes adult neuronal ceroid lipofuscinosis (CLN4), a fatal neurodegenerative disease. This pipeline identifies cell-type-specific transcriptional responses and dysregulated pathways downstream of CSP loss.

---

## Pipeline Architecture

```mermaid
flowchart TD
    A["📥 10X CellRanger output · filtered_feature_bc_matrix"] --> B

    subgraph QC ["1 · Quality Control & Doublet Detection"]
        B["Load with Seurat::Read10X · 4 samples per line"]
        B --> C["QC metrics · nCount_RNA · nFeature_RNA · %MT · Complexity"]
        C --> D["Doublet detection · scDblFinder · ~8% removal per sample"]
        D --> E["Cell filtering · UMI · gene count · MT thresholds"]
    end

    subgraph NORM ["2 · Normalisation & Integration"]
        E --> F["SCTransform · regularised negative binomial regression"]
        F --> G["Multi-sample integration · batch correction · HVG selection"]
    end

    subgraph DIM ["3 · Dimensionality Reduction & Clustering"]
        G --> H["PCA · 30 components · elbow selection"]
        H --> I["UMAP · 2D projection · cluster exploration"]
        I --> J["Louvain clustering · FindNeighbors k=30 · res=0.8"]
    end

    subgraph DE ["4 · Differential Expression"]
        J --> K["FindAllMarkers · Wilcoxon · min.pct=0.1 · log2FC > 0.25"]
        K --> L["Volcano plots per cluster · UP / DOWN gene lists"]
    end

    subgraph ENRICH ["5 · Functional Enrichment"]
        L --> M["GO · enrichGO · gseGO · BP · CC · MF"]
        L --> N["KEGG · enrichKEGG · gseKEGG · pathview"]
        L --> O["STRING PPI networks · PANTHER · EnrichR"]
    end

    subgraph OUT ["6 · Output"]
        M & N & O --> P["Publication figures · TIFF + PDF · Excel tables"]
    end

    style QC fill:#1e3a5f,color:#fff,stroke:#3b82f6
    style NORM fill:#1e3a1e,color:#fff,stroke:#22c55e
    style DIM fill:#2a1e3a,color:#fff,stroke:#8b5cf6
    style DE fill:#3a1e1e,color:#fff,stroke:#ef4444
    style ENRICH fill:#3a2a1e,color:#fff,stroke:#f59e0b
    style OUT fill:#1e2a2a,color:#fff,stroke:#06b6d4
```

---

## Key Numbers — UBC-Cre Dataset

| Step | Value |
|---|---|
| Input cells (4 samples) | ~32,000 |
| Doublet rate (scDblFinder) | ~8% per sample |
| Cells retained after QC | ~30,000 |
| Clusters identified | 36 |
| Differential expression tests | Wilcoxon, per cluster vs rest |
| Enrichment methods | GO (ORA + GSEA) · KEGG · STRING · PANTHER · EnrichR |
| Output figures | 190+ (TIFF + PDF) |

---

## Example Output

### Quality Control

| | |
|---|---|
| ![DNAJC5 expression](docs/images/002_Dnajc5_Expression.png) | ![Genes vs cells](docs/images/006_Genes_vs_Cell.png) |
| *DNAJC5 expression across conditions* | *QC — genes detected vs UMI counts* |
| ![Complexity](docs/images/008_Complexity.png) | ![Elbow plot](docs/images/034_Elblowplot.png) |
| *Cell complexity metric* | *Elbow plot — PCA dimension selection* |

### Clustering & Differential Expression

| | |
|---|---|
| ![UMAP clusters](docs/images/050_UMAP_Dims_30_Res_0.8.png) | ![UMAP WT vs KO](docs/images/057_UMAP_ident.png) |
| *UMAP — 36 clusters at 30 PCs / res=0.8* | *UMAP — WT vs KO cell distribution* |
| ![Integrated UMAP](docs/images/056_UMAP.png) | ![Heatmap cluster markers](docs/images/062_HeatMap_cluster_plot.png) |
| *Integrated UMAP across 4 samples* | *Heatmap — top markers per cluster* |

---

## Repository Structure

```
.
├── README.md
├── renv.lock                          # Dependency lock file (R 4.5.2)
├── single_cell.Rproj                  # RStudio project
│
├── code/                              # Modular R scripts
│   ├── 00_packages.R                  # Dependency management (pak)
│   ├── 01_sc_functions.R              # Core QC utilities & plot export
│   ├── 02_vulcano_plots.R             # Volcano plot generation
│   ├── 03_GO.R                        # GO over-representation analysis
│   ├── 04_strings.R                   # STRING PPI network analysis
│   ├── 05_gse.R                       # GSEA (GO + KEGG ranked lists)
│   ├── 06_Heatmap.R                   # ComplexHeatmap visualisation
│   ├── 07_HeatMap_GO_types.R          # GO-category heatmaps
│   ├── 08_EnrichR.R                   # EnrichR enrichment
│   ├── 09_gseKEGG.R                   # KEGG pathway GSEA
│   ├── ABA_sc_ref.R                   # Allen Brain Atlas reference
│   ├── Clusters_splitted_libraries.R  # Per-library independent clustering
│   ├── Doublets_Finders.R             # scDblFinder doublet detection
│   └── global_variables.R             # Thresholds & organism parameters
│
├── rmds/                              # R Markdown analysis notebooks
│   ├── Single_Cell_10X_Integrated_functions_SCT - UBC_Cre.Rmd   # UBC-Cre (pan-neuronal)
│   ├── Single_Cell_10X_Integrated_functions_SCT - PV_Cre.Rmd    # PV-Cre (GABAergic)
│   ├── Clustering Association_FindAllMarkers.Rmd                  # Cluster annotation
│   └── ...                            # Additional integration/merge notebooks
│
└── docs/
    └── images/                        # Representative output figures (PNG)
```

---

## R Scripts Reference

### Core pipeline (`code/`)

| Script | Purpose |
|---|---|
| `00_packages.R` | Install/load all dependencies via `pak` |
| `01_sc_functions.R` | `library_summary()`, `generate_qc_plots()`, `save_plot()` — QC metrics and dual-format export |
| `02_vulcano_plots.R` | `perform_vulcano()` — ggplot2 volcano plots with ggrepel labels |
| `03_GO.R` | `perform_enrichGO()` — clusterProfiler GO over-representation, BP/CC/MF |
| `04_strings.R` | STRINGdb PPI network retrieval and visualisation |
| `05_gse.R` | `process_gene_list()` + GSEA ranked-list pipeline for GO and KEGG |
| `06_Heatmap.R` | ComplexHeatmap of top DEGs per cluster; tracks DNAJC5 expression |
| `07_HeatMap_GO_types.R` | GO-category-specific expression heatmaps |
| `08_EnrichR.R` | Multi-library enrichment via enrichR (GO, KEGG, Reactome, WikiPathways) |
| `09_gseKEGG.R` | `gseKEGG()` with pathview pathway diagrams |
| `Doublets_Finders.R` | scDblFinder per-sample doublet scoring and removal |
| `Clusters_splitted_libraries.R` | Independent per-sample UMAP + Louvain clustering |
| `ABA_sc_ref.R` | Allen Brain Atlas integration for cell type reference annotation |
| `global_variables.R` | Centralised thresholds: `p_val`, `FC`, `kegg_organism`, `species` |

### Configuration (`global_variables.R`)

```r
p_val          <- 0.05          # Adjusted p-value threshold
FC             <- 0.25          # log2FC threshold
kegg_organism  <- "mmu"         # Mus musculus KEGG code
species        <- 10090         # NCBI taxonomy ID
organism       <- "org.Mm.eg.db"
keyType        <- "UNIPROT"
```

---

## Reproducing the Analysis

### 1. Restore the R environment

```r
# Install renv if needed
install.packages("renv")

# Restore all packages from the lockfile (R 4.5.2)
renv::restore()
```

### 2. Configure paths and run

Open the appropriate RMarkdown notebook and set the `data_path` variable to point to your CellRanger output directory:

```r
# UBC-Cre (pan-neuronal) analysis
# Expects rawdata/SC/{sample}_results/filtered_feature_bc_matrix/
rmarkdown::render("rmds/Single_Cell_10X_Integrated_functions_SCT - UBC_Cre.Rmd")

# PV-Cre (GABAergic) analysis
# Expects rawdata/Pv_nuclei/{sample}/filtered/
rmarkdown::render("rmds/Single_Cell_10X_Integrated_functions_SCT - PV_Cre.Rmd")
```

### 3. Run enrichment scripts

Enrichment scripts are sourced from within the RMarkdown notebooks. Individual modules can be run independently:

```r
source("code/03_GO.R")       # GO enrichment
source("code/05_gse.R")      # GSEA
source("code/09_gseKEGG.R")  # KEGG pathway analysis
source("code/04_strings.R")  # STRING PPI networks
```

### Data note

Raw 10X CellRanger outputs (~70 GB) and processed Seurat objects are not versioned. The `renv.lock` file fully specifies the computational environment.

---

## Tech Stack

| Layer | Tools |
|---|---|
| Single-cell framework | Seurat v5, sctransform |
| Doublet detection | scDblFinder |
| Cell type annotation | SingleR, Azimuth, Allen Brain Atlas |
| Differential expression | FindAllMarkers (Wilcoxon), MAST |
| GO enrichment | clusterProfiler (enrichGO, gseGO) |
| KEGG analysis | clusterProfiler (enrichKEGG, gseKEGG), pathview |
| PPI networks | STRINGdb |
| Multi-database enrichment | enrichR (GO, KEGG, Reactome, WikiPathways) |
| Pathway classification | rbioapi (PANTHER) |
| Visualisation | ggplot2, ComplexHeatmap, patchwork |
| Annotation | org.Mm.eg.db, biomaRt, AnnotationHub |
| Reproducibility | renv (lockfile) |

---

## Author

**Santiago López Begines, PhD**
Neuroscientist → Data Scientist
[Portfolio](https://slopezbegines.github.io/projects/single-cell/) · [GitHub](https://github.com/SLopezBegines) · [LinkedIn](https://linkedin.com/in/santibegines) · [ORCID](https://orcid.org/0000-0001-8809-8919)

---

## License

Code available for educational and research purposes with attribution. Raw data not included. Biological datasets are the property of the originating research group.
