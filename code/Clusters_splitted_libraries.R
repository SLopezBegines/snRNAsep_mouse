output_path <- paste0(output_path, "splitted_objects/")
dir.create(output_path, recursive = TRUE,showWarnings = FALSE)

DefaultAssay(SCT_seurat_objects$c1) <- "SCT"
DefaultAssay(SCT_seurat_objects$c2) <- "SCT"
DefaultAssay(SCT_seurat_objects$ko1) <- "SCT"
DefaultAssay(SCT_seurat_objects$ko1) <- "SCT"

# Supongamos que ya tienes 4 objetos Seurat (tras SCTransform) en una lista:
# Por ejemplo: lista_seurat <- list(lib1, lib2, lib3, lib4)

# 1. Realizar PCA, UMAP, encontrar vecinos y clusters en cada objeto
SCT_seurat_objects <- lapply(SCT_seurat_objects, function(obj){
  # Ejecutamos PCA (sin imprimir mensajes)
  #obj <- RunPCA(obj, verbose = TRUE)
  # Ejecutamos UMAP usando los primeros 30 PCs (puedes ajustar este parámetro)
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", verbose = TRUE)
  # Encontramos vecinos y definimos clusters
  obj <- FindNeighbors(obj, dims = 1:30, verbose = TRUE)
  obj <- FindClusters(obj, resolution = 0.8, verbose = TRUE)
    return(obj)
})

# 2. Visualizar el UMAP de cada librería
# Usaremos purrr::imap para tener acceso al índice o nombre de cada objeto
plots <- purrr::imap(SCT_seurat_objects, ~{
  DimPlot(.x, reduction = "umap", label = TRUE,pt.size = 0.1) +
    ggtitle(paste("UMAP de la librería", .y))
})
save_plot("UMAP_c1", plots$c1)
save_plot("UMAP_c2", plots$c2)
save_plot("UMAP_ko1", plots$ko1)
save_plot("UMAP_ko2", plots$ko2)


# 3. Mostrar los gráficos
# Puedes imprimirlos uno a uno o combinarlos en un panel, por ejemplo, usando patchwork
# Imprimir de forma individual:
print(plots[[1]])
print(plots[[2]])
print(plots[[3]])
print(plots[[4]])

# Opcional: combinar en un solo panel (requiere el paquete patchwork)
# install.packages("patchwork")  # Si no lo tienes instalado
library(patchwork)
plot_wrapped <- wrap_plots(plots)
save_plot("UMAP_all", plot_wrapped)
print(plot_wrapped)
  
## FindAllMarkers

plan("sequential")  # Usa 6 núcleos (ajusta según tu PC)
options(future.globals.maxSize = 14000 * 1024^2) # 16 GB limit
gc()
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
#DefaultAssay(data_select) <- "integrated" #SCT
#data_select <- PrepSCTFindMarkers(data_select)
Idents(SCT_seurat_objects$c1) <- "seurat_clusters"  # Set cluster identities. Look for differences between cluster. Dont look for differences between genotype.
Idents(SCT_seurat_objects$c2) <- "seurat_clusters"
Idents(SCT_seurat_objects$ko1) <- "seurat_clusters"
Idents(SCT_seurat_objects$ko2) <- "seurat_clusters"




#Idents(data_select) <- "condition"  # Set condition identities. Look for differences between genotype and dont look for differenes between clusters


# 1. Encontrar los marcadores de todos los clusters para cada objeto en la lista
markers_list <- lapply(SCT_seurat_objects, function(obj){
  FindAllMarkers(obj, 
                 test.use = "wilcox",      # Test estadístico: wilcox, wilcox_limma, LR, MAST, DESeq2
                 only.pos = FALSE, 
                 min.pct = 0.1, 
                 logfc.threshold = 0.1)
})

# Si tus objetos tienen nombres (por ejemplo, "lib1", "lib2", ...), puedes asignarlos a markers_list
names(markers_list) <- names(SCT_seurat_objects)

# 2. Escribir cada uno de los resultados en un archivo Excel separado
# Se utilizará purrr::imap para iterar sobre markers_list y usar el nombre o índice (.y) en el nombre del archivo
imap(markers_list, ~{
  output_file <- paste0(output_path, "tables/cluster.all.markers_", .y, ".xlsx")
  write_xlsx(.x, path = output_file)
})



# Cerebellum markers ####


# Leer el archivo Excel (asumiendo que la primera hoja contiene los datos)
file_path <- "rawdata/cb_cell_types.xlsx"
data <- read_excel(file_path)

# Transformar los datos: agrupar subtipos similares y pivotar
cerebellum_data <- data %>%
  dplyr::select(gene, subtype, avgExpr) %>%  # Seleccionar solo las columnas necesarias
  dplyr::mutate(subtype = str_replace(subtype, "_\\d+$", "")) %>%  # Eliminar sufijos _1, _2, etc.
  dplyr::group_by(gene, subtype) %>%  # Agrupar por gen y subtipo
  dplyr::summarise(avgExpr = sum(avgExpr, na.rm = TRUE), .groups = "drop") %>%  # Sumar valores de expresión
  pivot_wider(
    names_from = subtype,  # Usar la columna "subtype" como nombres de columnas
    values_from = avgExpr,  # Usar los valores de la columna "avgExpr" como contenido de las celdas
    values_fill = list(avgExpr = 0)  # Asignar 0 a los valores NA
  )%>%
  dplyr::mutate(Purkinje = coalesce(Purkinje_Aldoc, 0) + coalesce(Purkinje_Anti_Aldoc, 0)) %>%  # Combinar columnas Purkinje_Aldoc y Purkinje_Anti_Aldoc
  dplyr::select(-Purkinje_Aldoc, -Purkinje_Anti_Aldoc)  # Eliminar las columnas originales
# Verificar las primeras filas
head(cerebellum_data)
# Revisar dimensiones
dim(cerebellum_data)

# Convertir a matriz
cerebellum_matrix <- as.matrix(cerebellum_data[ , -1, with = FALSE])  # Asumiendo que la primera columna son los nombres de genes
rownames(cerebellum_matrix) <- cerebellum_data[[1]]  # Asignar nombres de genes a las filas

# Crear referencia de SingleR
# Crear un nuevo objeto SummarizedExperiment con logcounts
reference_cerebellum <- SummarizedExperiment(
  assays = list(logcounts = log1p(cerebellum_matrix)), # Aplicar log-transformation
  rowData = DataFrame(gene = rownames(cerebellum_matrix)),
  colData = DataFrame(label.main = colnames(cerebellum_matrix))
)
# Revisar las etiquetas
colData(reference_cerebellum)

# Confirmar que `label.main` contiene los nombres de los tipos celulares
reference_cerebellum$label.main <- colnames(cerebellum_matrix)


# Create markers list from FindAllMarkers ####

clusters_markers <- purrr::imap(markers_list, function(marker_df, name) {
  
  # Extraer los genes significativos encontrados por FindAllMarkers
  significant_genes <- marker_df %>% 
    dplyr::filter(p_val_adj < 0.05) %>% 
    dplyr::select(gene) %>% 
    dplyr::distinct() %>% 
    dplyr::pull(gene)
  
  # Obtener el objeto Seurat correspondiente usando el nombre
  seurat_obj <- SCT_seurat_objects[[name]]
  
  # Extraer la matriz de expresión para los genes significativos.
  # Se utiliza 'slot = "data"' para obtener los datos normalizados/log-transformados.
  expression_matrix <- GetAssayData(seurat_obj, slot = "data")[significant_genes, ]
  
  # Extraer las identidades (clusters) de las células
  cluster_ids <- Idents(seurat_obj)
  
  # Crear un objeto SummarizedExperiment con la matriz de expresión
  cluster_markers_se <- SummarizedExperiment(
    assays = list(logcounts = as.matrix(expression_matrix)),
    rowData = DataFrame(gene = rownames(expression_matrix)),
    colData = DataFrame(cluster = cluster_ids)
  )
  
  return(cluster_markers_se)
})


## Comparar con base de datos de referencia ####


# Usando lapply para iterar sobre la lista clusters_markers
annotations_cerebellum_list <- lapply(clusters_markers, function(cluster_marker) {
  
  # Extraer los identificadores de clúster del objeto SummarizedExperiment
  cluster_ids <- colData(cluster_marker)$cluster
  
  # Aplicar SingleR con el objeto SummarizedExperiment, el objeto de referencia y los labels correspondientes
  annotation <- SingleR(test = cluster_marker, 
                        ref = reference_cerebellum, 
                        labels = reference_cerebellum$label.main,
                        clusters = cluster_ids)
  
  return(annotation)
})


# Cargar librerías necesarias
library(Seurat)
library(tidyverse)  # Para dplyr, purrr, etc.
library(ggplot2)

# Se asume que:
# SCT_seurat_objects: lista nombrada de objetos Seurat procesados.
# annotations_cerebellum_list: lista nombrada de resultados SingleR, con nombres coincidentes.
# (Opcional) cluster_markers_df_list: lista nombrada de data.frames con información de marcadores, 
#   donde cada data.frame contiene al menos una columna "cluster" con IDs numéricos.
# output_path, tiff_extension, pdf_extension, image_number están definidos previamente.

# Procesar cada objeto de la lista SCT_seurat_objects
plots_list <- purrr::imap(SCT_seurat_objects, function(data_select, name) {
  
  # 1. Obtener la anotación SingleR correspondiente
  annotation <- annotations_cerebellum_list[[name]]
  
  # 2. Renombrar los clústeres en el objeto Seurat según las etiquetas asignadas por SingleR
  new.cluster.cerebellum <- annotation$labels
  names(new.cluster.cerebellum) <- levels(data_select)
  data_select_cerebellum <- RenameIdents(data_select, new.cluster.cerebellum)
  
  # 3. Convertir la anotación a data.frame para luego unirla con el data.frame de marcadores (si lo tienes)
  annotations_cerebellum_df <- as.data.frame(annotation)
  annotations_cerebellum_df$cluster <- rownames(annotations_cerebellum_df)  # Añade los IDs de clúster
  
  # 4. (Opcional) Actualizar el data.frame de marcadores con la etiqueta asignada por SingleR
  # Se asume que cada objeto tiene su correspondiente data.frame en cluster_markers_df_list
  if (exists("cluster_markers_df_list") && !is.null(cluster_markers_df_list[[name]])) {
    cluster_markers_df <- cluster_markers_df_list[[name]]
    
    # Suponemos que la columna "cluster" en cluster_markers_df es numérica y comienza en 0.
    # Se extraen las etiquetas (cell types) de la anotación y se asignan a cada clúster.
    cell_types_cerebellum <- annotations_cerebellum_df$labels
    cluster_markers_df$labels_cerebellum <- cell_types_cerebellum[cluster_markers_df$cluster + 1]
    
    # Actualizar la lista (opcional)
    cluster_markers_df_list[[name]] <<- cluster_markers_df
  }
  
  # 5. Visualizar el UMAP con las nuevas identidades
  plot <- DimPlot(data_select_cerebellum, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle(paste("UMAP -", name))
  print(plot)
  
  # 6. Guardar la figura en TIFF y PDF
  filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), "_UMAP_cerebellum_cluster_plot_", name)
  ggsave(filename = paste0(filename, tiff_extension), plot = plot, width = 12, height = 6, units = "in")
  ggsave(filename = paste0(filename, pdf_extension), plot = plot, width = 12, height = 6, units = "in")
  
  # 7. Incrementar el contador global de imagenes
  image_number <<- image_number + 1
  
  # Devolver el objeto Seurat modificado (con las identidades renombradas)
  return(data_select_cerebellum)
})



# Al finalizar, 'plots_list' contendrá los objetos Seurat modificados de cada librería.



features = c("Ppp1r17", "Pcp2", "Aldoc","L7","Gabra6", "Eomes", "Lypd6", "Prkcd", "Klhl1", "Lgi2", "Gdf10", "Aqp4", "Mobp", "Ppfibp1", "Dcn", "Kcnj8", "Ttr", "Mrc1", "C1qa", "Flt1", "Foxj1", "Calb1", "Pvalb", "Sst", "Calb2", "Tph2",  "Gfap", "Plp1","Mbp","Aif1", "Cdk11b", "Cx3cr1", "Pecam1", "Cldn5", "Slc17a7", "Lamp5", "Cux2", "Ndst4", "Rorb", "Tox", "Hs3st2", "Tle4", "Fexf2", "Gad1",  "Vip", "Cacng4", "Siglech")


markers_list <- lapply(SCT_seurat_objects, function(obj) {
  
  features_presents <- features[features %in% rownames(obj)]
  print(features_presents)
  
  for (feature in features_presents) {
    expr_values <- FetchData(obj, vars = feature)[,1]
    
    if (length(unique(expr_values)) == 1) {
      warning(paste("Skipping", feature, "because all cells have the same value:", unique(expr_values)))
      next
    }
    
    # VlnPlot
    plot <- VlnPlot(obj, features = feature, split.by = "library")
    print(plot)
    save_plot(paste0("Vln_plot_", feature), plot)
    
    # FeaturePlot
    feature_plot <- FeaturePlot(obj, features = feature) 
    print(feature_plot)
    save_plot(paste0("Feature_plot_Markers_", feature), feature_plot)
  }
})
