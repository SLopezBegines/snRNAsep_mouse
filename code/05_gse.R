

# libraries ####

source("./code/00_packages.R")
# Global variables
tiff_extension <- ".tiff"
pdf_extension <- ".pdf" #Vectorial format

# Setting thresholds ####
# Setting threshold for p-value and Fold-Change
p_val <- 0.05
FC <- 0.1


#Load data results tables
excel_sheet_reader <- function(filename) {
  # Get the names of all the sheets in the Excel file
  sheets <- readxl::excel_sheets(filename)
  # Read in each sheet, rename the first column to "ID", and store it in a list
  sheet_data <- lapply(sheets, function(sheet_name) {
    df <- readxl::read_excel(filename, sheet = sheet_name) %>%
      dplyr::rename(ID = ...1) %>%  # Rename the first column to "ID"
      dplyr::select(ID, everything())  # Reorder the columns with "ID" as the first column
    return(df)
  })
  # Combine the list of data frames into a single named list
  result <- setNames(sheet_data, sheets)
  # Return the list
  return(result)
}
# Call the function with the filename as the argument
input_data <- excel_sheet_reader(paste0(output_path,"tables/","cluster_list.xlsx"))

#Filtering results
# Function to filter a single data frame
filter_dataframe <- function(df) {
  df %>%
    filter(abs(avg_log2FC) > FC, p_val_adj < 0.05)  %>%
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
}
filtered_list <- map(input_data, filter_dataframe)

# Function to split the filtered data frame into two based on avg_log2FC values
split_dataframe <- function(filtered_df) {
  up_df <- filtered_df %>%
    filter(avg_log2FC > 0) %>%
    mutate(direction = "UP")
  
  down_df <- filtered_df %>%
    filter(avg_log2FC < 0) %>%
    mutate(direction = "DOWN")
  
  return(list(UP = up_df, DOWN = down_df))
}
# Apply the splitting function to each filtered data frame
split_filtered_list <- map(filtered_list, split_dataframe)
split_filtered_list <- unlist(split_filtered_list, recursive = FALSE)# Remove one level in the list

# Define a function to process each element in the KEGG_input_dataframes list
process_gene_list <- function(df) {
  # Extract foldchange and ID columns
  original_gene_list <- df$avg_log2FC
  names(original_gene_list) <- df$ID
  
  # Remove NA values
  gene_list <- na.omit(original_gene_list)
  
  # Sort the gene_list in decreasing order
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}
# Apply the function to each element in the KEGG_input_dataframes list
GSE_gene_lists <- lapply(split_filtered_list, process_gene_list)
# Gene Ontology ####
# GO functions ####
### EnrichGO functions ####
# GO terms list
# Function to enrichGO for the dataframes of a list
perform_gseGO <- function(gene_set) {
  results <- clusterProfiler::gseGO(geneList = gene_set,
                                    ont = "ALL",
                                    keyType = "SYMBOL", 
                                    nPerm = 10000, 
                                    minGSSize = 3, 
                                    maxGSSize = 800, 
                                    pvalueCutoff = 0.05, 
                                    verbose = TRUE, 
                                    OrgDb = "org.Mm.eg.db", 
                                    pAdjustMethod = "none")
  return(results)
}



# Apply the enrichGO function to each gene set in differential_genes
results_list_gseGO_ALL <- lapply(GSE_gene_lists, perform_gseGO)


save(results_list_gseGO_ALL, file = paste0(output_path,"RData/results_list_gseGO_ALL.RData"))

# load results_list_go_ALL RData file
#load("./data/output/RData/SC/results_list_gseGO_ALL.RData")


## Plots ####

# GO Dotplots function
perform_dotplots <- function(x) {
  if (length(x[1][10]) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    #results <- enrichplot::dotplot(x, showCategory = 30, split=".sign") + facet_grid(ONTOLOGY ~ ., scale = "free")+ scale_color_viridis()
    results <- enrichplot::dotplot(x, showCategory = 10, split=".sign") + 
      facet_grid(.~.sign) + 
      scale_color_viridis() 
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}


# GO cnetplots function
perform_cnetplots <- function(x){
  ## remove redundent GO terms
  x2 <- simplify(x)
  #results <- enrichplot::cnetplot(x2)
  #For circular Gene Concept map
  results <- enrichplot::cnetplot(x2, circular = TRUE, colorEdge = TRUE)
  #results <- enrichplot::cnetplot(x2, foldChange=x@geneList, circular = TRUE, colorEdge = TRUE)
  return(results)
}

#HeatMap plot function
perform_heatmapplot <- function(x){
  results <- enrichplot::heatplot(x, foldChange=x@geneList)
  return(results)
}

# RidgePlots
perform_ridgeplots <- function(x) {
  if (length(x[1][1]) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    #results <- enrichplot::dotplot(x, showCategory = 30, split=".sign") + facet_grid(ONTOLOGY ~ ., scale = "free")+ scale_color_viridis()
    results <- enrichplot::ridgeplot(x) +  
      scale_fill_viridis()
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}
# PMCPlots
perform_pmcplots <- function(x) {
  if (length(x[1][1]) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    #results <- enrichplot::dotplot(x, showCategory = 30, split=".sign") + facet_grid(ONTOLOGY ~ ., scale = "free")+ scale_color_viridis()
    terms <- x@result$Description[1:3]
    results <- enrichplot::pmcplot(terms, 2010:2023, proportion=FALSE)
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}

image_number = 231

### Dotplot ####

dotplot_results <- lapply(results_list_gseGO_ALL, perform_dotplots)
# Iterate over dotplot_results and print each plot individually
for (i in seq_along(dotplot_results)) {
  if (nrow(dotplot_results[[i]][[1]][10]) == 0) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- dotplot_results[[i]] +
      ggplot2::ggtitle(paste0("Dot plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO/",image_number +i,"_GSE_Dotplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    #ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
    #ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  }
}
image_number <- image_number +i

### CNET plot ####

cnetplot_results <- lapply(results_list_gseGO_ALL, perform_cnetplots)

for (i in seq_along(cnetplot_results)) {
  if (is.list(cnetplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- cnetplot_results[[i]]+
      ggplot2::ggtitle(paste0("CNET plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO/",image_number+i,"_GSE_cnet_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  }
}
image_number <- image_number+i




### Heatmap plot ####
heatmapplot_results <- lapply(results_list_gseGO_ALL, perform_heatmapplot)

for (i in seq_along(heatmapplot_results)) {
  if (is.list(heatmapplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- heatmapplot_results[[i]]+
      ggplot2::ggtitle(paste0("HeatMap plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO/",image_number+i,"_GSE_heatmap_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

### Ridge plot ####
ridgeplot_results <- lapply(results_list_gseGO_ALL, perform_ridgeplots)

for (i in seq_along(ridgeplot_results)) {
  if (is.list(ridgeplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- ridgeplot_results[[i]]+
      ggplot2::ggtitle(paste0("Ridge plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO/",image_number+i,"_GSE_ridgeplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

### PMCplot ####
# Number of articles in pubmed with the input ges descriptions
pmcplot_results <- lapply(results_list_gseGO_ALL, perform_pmcplots)

for (i in seq_along(pmcplot_results)) {
  if (is.list(pmcplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- pmcplot_results[[i]]+
      ggplot2::ggtitle(paste0("PMC plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO/",image_number+i,"_GSE_pmcplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

