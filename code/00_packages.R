# Function to load libraries from CRAN using pak
install_and_load_library <- function(lib_names) {
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
  }

  missing_libs <- lib_names[!sapply(lib_names, requireNamespace, quietly = TRUE)]

  if (length(missing_libs) > 0) {
    pak::pkg_install(missing_libs, dependencies = TRUE)
  }

  # Load all libraries
  for (lib in lib_names) {
    library(lib, character.only = TRUE)
  }
}

# Function to load libraries from Bioconductor using pak
install_and_load_library_bioconductor <- function(lib_names) {
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
  }

  missing_libs <- lib_names[!sapply(lib_names, requireNamespace, quietly = TRUE)]

  if (length(missing_libs) > 0) {
    # pak handles Bioconductor packages automatically with bioc:: prefix
    bioc_libs <- paste0("bioc::", missing_libs)
    pak::pkg_install(bioc_libs, dependencies = TRUE, ask = FALSE)
  }

  # Load all libraries
  for (lib in lib_names) {
    library(lib, character.only = TRUE)
  }
}

# Usage example
libraries <- c(
  "BiocManager", "R.utils", "Seurat", "plan", "tidyverse", "data.table",
  "kableExtra", "ggrepel", "promises", "furrr", "purrr",
  "readxl", "future", "future.apply", "sctransform",
  "presto", "BPCells", "SeuratData", "SeuratData", "Azimuth",
  "SeuratWrappers", "glmGamPoi", "Signac", "writexl", "patchwork",
  "remotes", "tictoc", "openxlsx", "clustree"
)

libraries_bioconductor <- c(
  "DEP", "SummarizedExperiment", "SingleR", "MAST",
  "genekitr", "clusterProfiler", "enrichplot",
  "AnnotationDbi", "topGO", "STRINGdb", "AnnotationHub",
  "rrvgo", "europepmc", "Rgraphviz", "pathview", "DOSE",
  "viridis", "gson", "tidytree", "celldex", "glmGamPoi",
  "biomaRt", "org.Mm.eg.db", "limma", "ComplexHeatmap",
  "edgeR", "STRINGdb", "scDblFinder", "KEGGREST"
)

install_and_load_library(libraries)
install_and_load_library_bioconductor(libraries_bioconductor)

# Install GitHub packages with remotes (more stable for complex builds)
github_packages <- c(
  "satijalab/seurat-data",
  "satijalab/azimuth",
  "satijalab/seurat-wrappers",
  "mojaveazure/seurat-disk",
  "bnprks/BPCells/r"
)

# Other GitHub packages
pak::pkg_install("satijalab/seurat-data")
pak::pkg_install("satijalab/azimuth")
pak::pkg_install("satijalab/seurat-wrappers")
pak::pkg_install("mojaveazure/seurat-disk")
for (pkg in github_packages) {
  pkg_name <- basename(pkg)
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    remotes::install_github(pkg, upgrade = "never", quiet = TRUE)
  }
}

# Remove variables from GlobalEnvironment
rm(libraries, libraries_bioconductor, install_and_load_library, install_and_load_library_bioconductor)
