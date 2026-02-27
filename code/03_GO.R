# Librerías ####
source("./code/00_packages.R")

# Variables Globales ####
tiff_extension <- ".tiff"
pdf_extension <- ".pdf" # Formato vectorial
output_path <- "./output/" # Ruta base para guardar resultados

# Umbrales para el filtrado ####
p_val <- 0.05
FC <- 0.25


# Aplicar filtrado a todas las tablas
filtered_list <- purrr::map(cluster_list, filter_dataframe)
# Función para dividir en UP y DOWN
split_dataframe <- function(filtered_df) {
  list(
    UP = filtered_df %>% dplyr::filter(direction == "UP"),
    DOWN = filtered_df %>% dplyr::filter(direction == "DOWN")
  )
}

# Dividir cada elemento de filtered_list y asignar nombres significativos
split_filtered_list <- purrr::imap(filtered_list, function(df, cluster_name) {
  split_result <- split_dataframe(df)
  # Nombrar cada resultado como Cluster_X_UP o Cluster_X_DOWN
  purrr::imap(split_result, function(split_df, direction) {
    list_name <- paste0(cluster_name, "_", direction)
    split_df
  })
})

# Aplanar la lista de listas en una sola lista
split_filtered_list <- purrr::flatten(split_filtered_list)

# Asignar nombres significativos a los elementos de la lista final
names(split_filtered_list) <- unlist(purrr::imap(filtered_list, function(x, cluster_name) {
  c(paste0(cluster_name, "_UP"), paste0(cluster_name, "_DOWN"))
}))


# Enriquecimiento GO ####
perform_enrichGO <- function(gene_set) {
  if (nrow(gene_set) == 0) return(NULL)
  clusterProfiler::enrichGO(
    gene = gene_set$gene,
    OrgDb = "org.Mm.eg.db",
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "fdr"
  )
}

# Aplicar enriquecimiento GO
results_list_enrichGO_ALL <- purrr::map(split_filtered_list, perform_enrichGO)

# Filtrar elementos NULL
results_list_enrichGO_ALL <- purrr::compact(results_list_enrichGO_ALL)

# Extraer resultados válidos
df_list <- purrr::map(results_list_enrichGO_ALL, function(x) {
  if (!inherits(x, "enrichResult")) return(NULL) # Verifica que sea un enrichResult
  x@result # Extrae el componente result
})

# Eliminar resultados vacíos
df_list <- purrr::compact(df_list)

# Filtrar elementos NULL
results_list_enrichGO_ALL <- purrr::compact(results_list_enrichGO_ALL)

# Extraer resultados válidos
df_list <- purrr::map(results_list_enrichGO_ALL, function(x) {
  if (!inherits(x, "enrichResult")) return(NULL) # Verifica que sea un enrichResult
  x@result # Extrae el componente result
})

# Eliminar resultados vacíos
df_list <- purrr::compact(df_list)

# Funciones de visualización ####
plot_results <- function(results, plot_type, output_prefix) {
  purrr::walk2(results, names(results), function(res, name) {
    if (is.null(res) || nrow(res@result) == 0) {
      message("No data available for ", name)
      return(NULL)
    }
    
    plot <- switch(
      plot_type,
      "barplot" = barplot(res, showCategory = 30) + ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free"),
      "dotplot" = enrichplot::dotplot(res, showCategory = 30) + ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free"),
      #"cnetplot" = enrichplot::cnetplot(simplify(res), circular = TRUE, colorEdge = TRUE),
      "heatmap" = enrichplot::heatplot(res),
      "upsetplot" = enrichplot::upsetplot(res),
      stop("Invalid plot type")
    )
    
    filename <- paste0(output_prefix, "_", plot_type, "_", gsub(" ", "_", name), tiff_extension)
    ggsave(filename, plot, width = 8, height = 6, units = "in", dpi = 300)
    filename <- paste0(output_prefix, "_", plot_type, "_", gsub(" ", "_", name), pdf_extension)
    ggsave(filename, plot, width = 8, height = 6, units = "in", dpi = 300)
    print(plot)
  })
}

# Aplicar y guardar cada tipo de visualización
output_prefix <- paste0(output_path, "figures/enrichGO/")
plot_results(results_list_enrichGO_ALL, "barplot", output_prefix)
plot_results(results_list_enrichGO_ALL, "dotplot", output_prefix)
#plot_results(results_list_enrichGO_ALL, "cnetplot", output_prefix)
plot_results(results_list_enrichGO_ALL, "heatmap", output_prefix)
plot_results(results_list_enrichGO_ALL, "upsetplot", output_prefix)



# Lolliplots ####

# Make dataframe list from GO results #
results_df_enrichGO_ALL <- map(results_list_enrichGO_ALL, ~ .x@result)

results_df_enrichGO_ALL <- lapply(results_df_enrichGO_ALL, function(df) {
  df <- df %>%
    mutate(GeneRatio = strsplit(GeneRatio, "/") %>%
             map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2])))
  return(df)
})

lolliplot <- function(data_name, df_list, file_prefix = NULL) {
  non_empty_indices <- which(sapply(df_list, function(x) nrow(x) > 0))
  
  if (length(non_empty_indices) == 0) {
    message("All data frames are empty. No plots generated.")
    return(NULL)
  }
  
  if (length(non_empty_indices) < data_name) {
    message("Selected index is out of range. No plot generated.")
    return(NULL)
  }
  
  df <- df_list[[non_empty_indices[data_name]]]
  
  if (nrow(df) == 0) {
    message("Selected data frame is empty. No plot generated.")
    return(NULL)
  }
  
  plot_name <- gsub("^names_", "", names(df_list)[non_empty_indices[data_name]])  # Remove "names_"
  plot_name <- gsub("_vs_", " vs ", plot_name)  # Replace "_vs_" with " vs "
  plot_name <- gsub("_", " ", plot_name)  # Replace "_" with " "
  plot_name <- gsub("UP$", "UP", plot_name)  # Remove "UP" from the end
  
  plot <- df %>%
    dplyr::mutate(Description = reorder(Description, GeneRatio)) %>%
    top_n(15, GeneRatio) %>%
    ggplot(aes(x = Description,
               y = GeneRatio,
               colour = p.adjust)) +
    geom_segment(aes(x = Description,
                     xend = Description,
                     y = 0,
                     yend = GeneRatio)) +
    geom_point(aes(size = Count), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    facet_wrap(ONTOLOGY ~ ., scale = "free", ncol=1)+
    coord_flip() +
    theme_gray() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "Gene Ratio") + 
    ggtitle(paste("Plot for", plot_name))
  
  if (!is.null(file_prefix)) {
    tiff_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".tiff")
    pdf_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename = tiff_filename, plot = plot, device = "tiff",width = 12, height = 12, units = "in")
    ggsave(filename = pdf_filename, plot = plot, device = "pdf", width = 12, height = 12, units = "in")
  }
  print(plot)
  return(plot) 
}

# Apply lolliplot function
for (i in seq_along(results_df_enrichGO_ALL)) {
  lolliplot_result <- lolliplot(i, results_df_enrichGO_ALL, file_prefix = paste0(output_path,"figures/","enrichGO/Lolliplot_0",i))
}

