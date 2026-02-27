# There are two main tools to find Doublets:

# scDblFinder ####
# https://github.com/plger/scDblFinder

# ---------------------- #
# scRNAseq doublet detection tutorial
# ---------------------- #

#library(scDblFinder)
tic()
doublet_folder <- file.path(output_path, 'Doublet_detection') # subfolder for Doublet detection results
dir.create(doublet_folder, recursive = TRUE, showWarnings = FALSE) 

# Doublet detection =================================================================
## scDblFinder -------------------------------
# Run scDbltFinder
# Función para procesar cada objeto Seurat
# Lista de objetos Seurat

# Función para procesar cada objeto Seurat
process_seurat <- function(seu, seu_name) {
  # Convertir a SingleCellExperiment
  sce <- as.SingleCellExperiment(seu)
  
  # Ejecutar scDblFinder
  sce <- scDblFinder(sce)
  
  # Extraer resultados y agregarlos al objeto Seurat
  meta_scdblfinder <- sce@colData@listData %>% 
    as.data.frame() %>%
    dplyr::select(starts_with('scDblFinder'))
  
  rownames(meta_scdblfinder) <- sce@colData@rownames
  seu <- AddMetaData(object = seu, metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))
  
  # Generar resumen de dobletes
  doublets_summary <- seu@meta.data %>% 
    group_by(scDblFinder.class) %>% 
    summarise(total_count = n(), .groups = 'drop') %>% 
    ungroup() %>%
    mutate(percent = paste0(round(100 * total_count / sum(total_count), 2), '%'),
           Sample = seu_name)
  
  # Guardar resumen en archivo individual
  write.table(doublets_summary, file = file.path(doublet_folder, paste0(seu_name, '_scDblFinder_doublets_summary.txt')), 
              quote = FALSE, row.names = FALSE, sep = '\t')
  
  # Remover doublets y dejar solo singlets
  seu_clean <- subset(seu, subset = scDblFinder.class == 'singlet')
  
  # Devolver lista con objeto filtrado, resumen y el objeto sin filtrar
  return(list(seurat = seu_clean, doublets_summary = doublets_summary, raw_seurat = seu))
}

# Aplicar la función a cada objeto en la lista con lapply
results <- lapply(names(filter_seurat_objects), function(name) process_seurat(filter_seurat_objects[[name]], name))

# Extraer objetos Seurat filtrados
doublets_seurat_objects <- lapply(results, function(x) x$seurat)
names(doublets_seurat_objects) <- names(filter_seurat_objects)

# Extraer objetos Seurat sin filtrar (para gráficos)
raw_seurat_list <- lapply(results, function(x) x$raw_seurat)
names(raw_seurat_list) <- names(filter_seurat_objects)

# Extraer resúmenes de dobletes y combinarlos en un solo dataframe
doublets_summary_list <- lapply(results, function(x) x$doublets_summary)
combined_doublets_summary <- do.call(rbind, doublets_summary_list)

# Guardar el resumen combinado en un archivo
write.table(combined_doublets_summary, file = file.path(doublet_folder, "combined_scDblFinder_summary.txt"), 
            quote = FALSE, row.names = FALSE, sep = '\t')

# Mostrar el dataframe resumen en pantalla
print(combined_doublets_summary)


# Export results ####
save(doublets_seurat_objects, file = paste0(output_path, "RData/", "doublets_seurat_objects", ".RData"))  # Save using the name
load(paste0(output_path,"RData/doublets_seurat_objects.RData"))


# Función para generar y guardar violin plots
generate_violin_plots <- function(filter_seurat_objects, output_folder) {
  for (name in names(filter_seurat_objects)) {
    seu <- filter_seurat_objects[[name]]
    
    # Generar gráfico de violín
    plot <- VlnPlot(seu, split.by = "scDblFinder.class",
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal"), 
                    ncol = 3, pt.size = 0) + theme(legend.position = 'right')
    
    # Imprimir el gráfico en pantalla
    print(plot)
    
    # Guardar el gráfico
    save_plot(paste0(name,  "_violin_plot_Doublets"), plot)
  }
}

# Ejecutar la función para generar gráficos usando los objetos sin filtrar
generate_violin_plots(raw_seurat_list, doublet_folder)

# Ejecutar la función para generar gráficos
generate_violin_plots(doublets_seurat_objects, doublet_folder)
#----------------------------#

toc()

