# Preprocessing and Quality control

# 2. Data preprocessing for each library
## QC and selecting cells for further analysis
'Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include
## Identify mitochondrial genes
The number of unique genes detected in each cell.
Low-quality cells or empty droplets will often have very few genes
Cell doublets or multiplets may exhibit an aberrantly high gene count
Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
The percentage of reads that map to the mitochondrial genome
Low-quality / dying cells often exhibit extensive mitochondrial contamination
We calculate mitochondrial QC metrics with the "PercentageFeatureSet" function, which calculates the percentage of counts originating from a set of features
We use the set of all genes starting with "MT-" as a set of mitochondrial genes'


# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)  # For combining plots

#Save Plots

# Función para guardar gráficos en TIFF y PDF
save_plot <- function(plotname, plot, width = 8, height = 6) {
  # Verificar si plot es un objeto válido de ggplot
  if (!inherits(plot, "ggplot")) {
    warning(paste("El objeto", plotname, "no es un gráfico ggplot válido. Se omite la guardado."))
    return(NULL)
  }
  
  # Construcción del nombre de archivo con contador
  filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), "_", plotname)
  
  # Guardar en TIFF y PDF
  tryCatch({
    ggsave(paste0(filename, tiff_extension), plot, width = width, height = height, units = "in")
    ggsave(paste0(filename, pdf_extension), plot, width = width, height = height, units = "in")
    
    # Incrementar el contador global de imágenes
    image_number <<- image_number + 1
  }, error = function(e) {
    warning(paste("No se pudo guardar la imagen:", filename, "Error:", e$message))
  })
}



# QC summary
# Función para obtener métricas de calidad de cada objeto Seurat
library_summary <- function(seurat_list) {
  # Crear una lista para almacenar los resultados
  results <- lapply(names(seurat_list), function(name) {
    seurat_obj <- seurat_list[[name]]
    # Calcular las métricas
    number_of_reads <- sum(seurat_obj$nCount_RNA)
    estimated_nuclei <- ncol(seurat_obj)
    mean_reads_per_nucleus <- mean(seurat_obj$nCount_RNA)
    median_genes_per_nucleus <- median(seurat_obj$nFeature_RNA)
    
    # Retornar un dataframe para cada objeto
    data.frame(
      Library = name,
      Number_of_Reads = number_of_reads,
      Estimated_Number_of_Nuclei = estimated_nuclei,
      Mean_Reads_Per_Nucleus = mean_reads_per_nucleus,
      Median_Genes_Per_Nucleus = median_genes_per_nucleus
    )
  })
  
  # Combinar los resultados en un solo dataframe
  do.call(rbind, results)
}




# QC plots ####
generate_qc_plots <- function(seurat_object, seurat_object_name, image_number, output_path, tiff_extension = ".tiff", pdf_extension = ".pdf", suffix = "") {
  # Generate QC violin plot
  QC_VlnPlot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribosomal"), ncol = 2)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  # Generate feature scatter plots
  plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  plot3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.ribosomal")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  plot4 <- FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.mt")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  plot5 <- FeatureScatter(seurat_object, feature1 = "nFeature_RNA", feature2 = "percent.ribosomal")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  plot6 <- FeatureScatter(seurat_object, feature1 = "percent.mt", feature2 = "percent.ribosomal")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  combined_plot <- plot1 + plot2 + plot3 + plot4 + plot5 + plot6
  
  # Print the plots
  print(QC_VlnPlot)
  print(combined_plot)
  
  # Set filenames with suffix
  vln_filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), "_QC_VlnPlot_", suffix, "_", seurat_object_name)
  scatter_filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), "_Features_Plot_", suffix, "_", seurat_object_name)
  
  # Save plots
  ggsave(paste0(vln_filename, tiff_extension), QC_VlnPlot, width = 8, height = 8, units = "in")
  ggsave(paste0(vln_filename, pdf_extension), QC_VlnPlot, width = 8, height = 8, units = "in")
  ggsave(paste0(scatter_filename, tiff_extension), combined_plot, width = 8, height = 8, units = "in")
  ggsave(paste0(scatter_filename, pdf_extension), combined_plot, width = 8, height = 8, units = "in")
  
  # Return incremented figure counter
  return(image_number + 1)
}

# Normalization ####
# Define the normalization and scaling function
normalize_seurat_objects <- function(seurat_objects, log_file = "process_log.txt") {
  
  # Apply processing to each Seurat object in parallel
  norm_seurat_objects <- future_lapply(names(seurat_objects), function(name) {
    start_time <- Sys.time()
    obj <- seurat_objects[[name]]  # Access the Seurat object
    
    tryCatch({
      # Log start
      write(paste(Sys.time(), "- Processing:", name), file = log_file, append = TRUE)
      
      # Retrieve all genes
      all.genes <- rownames(obj)
      
      # Normalize, find variable features, and scale the data
      norm <- NormalizeData(obj, normalization.method = "LogNormalize", 
                            scale.factor = 10000) %>%
        FindVariableFeatures(selection.method = "vst", 
                             nfeatures = 2000)
      
      # Scale data
      norm <- ScaleData(norm, features = all.genes)
      
      # Save the normalized object
      save(norm, file = paste0(output_path, "RData/", "norm_", name, ".RData"))
      
      # Capture metadata
      processing_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      metadata_row <- data.frame(
        object_name = name,
        num_genes = nrow(norm),
        num_cells = ncol(norm),
        processing_time = processing_time,
        stringsAsFactors = FALSE
      )
      
      # Return both normalized object and metadata
      list(norm = norm, metadata = metadata_row)
      
    }, error = function(e) {
      # Log error
      write(paste(Sys.time(), "- Error processing:", name, ":", e$message), file = log_file, append = TRUE)
      NULL  # Return NULL in case of error
    })
  })
  
  # Separate normalized objects and metadata
  valid_results <- Filter(Negate(is.null), norm_seurat_objects)  # Remove NULL results
  norm_objects <- setNames(
    lapply(valid_results, `[[`, "norm"),
    names(seurat_objects)  # Use original names
  )
  
  metadata_list <- lapply(valid_results, `[[`, "metadata")
  
  # Combine metadata
  metadata <- do.call(rbind, metadata_list)
  
  # Print metadata to R Markdown
  print(metadata)
  
  # Save metadata as CSV
  write.csv(metadata, file = paste0(output_path, "tables/metadata.csv"), row.names = FALSE)
  
  return(norm_objects)
}




# HV Plots ####
# Function to generate highly variable (HV) gene plots
generate_HV_plots <- function(seurat_object, seurat_object_name, image_number, output_path, 
                              tiff_extension = ".tiff", pdf_extension = ".pdf") {
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(seurat_object), 10)
  
  # Plot variable features
  plot <- VariableFeaturePlot(seurat_object) %>% LabelPoints(points = top10, repel = TRUE) +
    ggtitle(paste("Highly Variable Genes: ", seurat_object_name))  # Add title with object name
  
  # Print plot to R Markdown
  print(plot)
  
  # Save plots
  filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), 
                     "_Highly_Variable_Genes_", seurat_object_name)
  ggsave(paste0(filename, tiff_extension), plot, width = 8, height = 6, units = "in")
  ggsave(paste0(filename, pdf_extension), plot, width = 8, height = 6, units = "in")
  
  # Return incremented figure counter
  return(image_number + 1)
}

# Add PCA ####
# Function to perform PCA and save visualizations
add_pca <- function(norm_seurat_objects, output_path, image_number, 
                    tiff_extension = ".tiff", pdf_extension = ".pdf") {
  # Create a list to store the updated Seurat objects
  updated_seurat_objects <- list()
  
  # Apply PCA and save plots for each Seurat object
  for (obj_name in names(norm_seurat_objects)) {
    obj <- norm_seurat_objects[[obj_name]]  # Retrieve Seurat object by name
    
    tryCatch({
      # Run PCA
      obj <- RunPCA(obj, features = VariableFeatures(object = obj))
      
      # Print PCA summary (first 5 dimensions and 5 features per dimension)
      print(obj[["pca"]], dims = 1:5, nfeatures = 5)
      
      # Visualize PCA loadings
      plot <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")
      print(plot)
      
      # Save plot to file
      filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), 
                         "_PCA_Loadings_", obj_name)
      ggsave(paste0(filename, tiff_extension), plot, width = 8, height = 6, units = "in")
      ggsave(paste0(filename, pdf_extension), plot, width = 8, height = 6, units = "in")
      
      # Increment the figure counter
      image_number <<- image_number + 1  # Update global figure counter
      
      # Add the updated object to the list
      updated_seurat_objects[[obj_name]] <- obj
      
    }, error = function(e) {
      message(paste("Error processing PCA for", obj_name, ":", e$message))
    })
  }
  
  # Ensure the returned list retains the names of the input objects
  if (length(updated_seurat_objects) < length(norm_seurat_objects)) {
    warning("Some objects failed during PCA processing and are not included in the output.")
  }
  
  return(updated_seurat_objects)
}

# Elbow Plots ####
# Function to generate and save elbow plots
# Function to generate and save elbow plots
elbow_plots <- function(seurat_objects, output_path, image_number, 
                                 tiff_extension = ".tiff", pdf_extension = ".pdf") {
  # Print the names of the input objects for debugging
  message("Input object names: ", paste(names(seurat_objects), collapse = ", "))
  
  # Loop through each Seurat object by name
  for (obj_name in names(seurat_objects)) {
    message("Processing object: ", obj_name)  # Debugging message
    
    obj <- seurat_objects[[obj_name]]  # Retrieve the Seurat object
    
    # Check if the object is NULL or invalid
    if (is.null(obj)) {
      message("Warning: Object ", obj_name, " is NULL or invalid. Skipping.")
      next  # Skip to the next object
    }
    
    # Try generating and saving the ElbowPlot
    tryCatch({
      # Debugging: Print object class to confirm it's valid
      message("Object class: ", class(obj))
      
      # Generate the ElbowPlot
      plot <- ElbowPlot(obj)
      
      # Display the plot for debugging
      print(plot)
      
      # Define file paths
      filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), 
                         "_ElbowPlot_", obj_name)
      
      # Save plot in TIFF and PDF formats
      ggsave(paste0(filename, tiff_extension), plot, width = 8, height = 6, units = "in")
      ggsave(paste0(filename, pdf_extension), plot, width = 8, height = 6, units = "in")
      
      # Increment the global figure counter
      image_number <<- image_number + 1
      message("Saved ElbowPlot for ", obj_name, " as ", filename)  # Debugging message
      
    }, error = function(e) {
      message("Error processing ", obj_name, ": ", e$message)  # Capture errors
    })
  }
  
  message("Elbow plot generation completed.")  # Final debugging message
}


# Clustering and UMAP ####

cluster_and_umap <- function(seurat_objects, output_path, image_number, 
                             tiff_extension = ".tiff", pdf_extension = ".pdf", 
                             dims = 1:11, resolution = 0.5, components =2L) {
  # Create a list to store the updated Seurat objects
  updated_seurat_objects <- list()
  
  # Loop through each Seurat object by name
  for (obj_name in names(seurat_objects)) {
    obj <- seurat_objects[[obj_name]]  # Retrieve Seurat object by name
    
    tryCatch({
      # Perform clustering
      obj <- FindNeighbors(obj, dims = dims)
      obj <- FindClusters(obj, resolution = resolution)
      
      # Run UMAP
      obj <- RunUMAP(obj, dims = dims, n.components = components)
      
      # Create UMAP plot and add the object name as the title
      plot <- DimPlot(obj, reduction = "umap", label = TRUE) + 
        ggtitle(paste("UMAP Plot -", obj_name))
      print(plot)
      
      # Save UMAP plot to file
      filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), 
                         "_UMAP_", obj_name)
      ggsave(paste0(filename, tiff_extension), plot, width = 8, height = 6, units = "in")
      ggsave(paste0(filename, pdf_extension), plot, width = 8, height = 6, units = "in")
      
      # Increment the figure counter
      image_number <<- image_number + 1  # Update global figure counter
      
      # Add the updated object to the list
      updated_seurat_objects[[obj_name]] <- obj
      
    }, error = function(e) {
      message(paste("Error processing clustering/UMAP for", obj_name, ":", e$message))
    })
  }
  
  # Ensure the returned list retains the names of the input objects
  if (length(updated_seurat_objects) < length(seurat_objects)) {
    warning("Some objects failed during clustering/UMAP and are not included in the output.")
  }
  
  return(updated_seurat_objects)
}

# Data Integration ####

integrate_and_cluster_parallel <- function(seurat_objects, seurat_names, output_path, image_number, 
                                           dims = 1:11, resolution = 0.5, components =2L,
                                           tiff_extension = ".tiff", pdf_extension = ".pdf", gene) {
  
  # Plan for parallelization (adjust workers based on system)
  plan(multicore, workers = 4)
  
  # Find integration anchors in parallel with a random seed
  anchors <- future({
    FindIntegrationAnchors(object.list = lapply(seurat_names, function(name) seurat_objects[[name]]), anchor.features = 3000)
  }, seed = TRUE)  # Safe parallel seed
  
  anchors <- value(anchors)  # Wait for anchors
  # Asegúrate de incluir Dnajc5
  anchors <- unique(c(anchors, gene))
  # Perform data integration with a random seed
  data_int <- future({
    IntegrateData(anchorset = anchors)
  }, seed = TRUE)  # Safe parallel seed
  
  data_int <- value(data_int)  # Wait for integration
  
  # Process the integrated data with a random seed
  data_int <- future({
    data_int <- ScaleData(data_int) %>% 
      RunPCA() %>% 
      RunUMAP(dims = dims, n.components = components) %>% 
      FindNeighbors(dims = dims) %>% 
      FindClusters(resolution = resolution)
    data_int
  }, seed = TRUE)  # Safe parallel seed
  
  data_int <- value(data_int)  # Wait for processing
  
  # Create and save UMAP plots in parallel with random seed
  future_lapply(seurat_names, function(obj_name) {
    obj <- seurat_objects[[obj_name]]  # Retrieve the corresponding Seurat object
    
    # Generate a UMAP plot
    plot <- DimPlot(data_int, reduction = "umap", label = TRUE) + 
      ggtitle(paste("UMAP Plot -", obj_name))
    print(plot)
    
    # Save the plot
    filename <- paste0(output_path, "figures/", sprintf("%03d", image_number), "_Integrated_UMAP_", obj_name)
    ggsave(paste0(filename, tiff_extension), plot, width = 8, height = 6, units = "in")
    ggsave(paste0(filename, pdf_extension), plot, width = 8, height = 6, units = "in")
    
    # Increment the figure counter
    image_number <<- image_number + 1
  }, future.seed = TRUE, 
  future.globals = c("seurat_objects", "data_int", "output_path", 
                     "image_number", "tiff_extension", "pdf_extension"))
  
  return(data_int)
}




