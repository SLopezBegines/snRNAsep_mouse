


# libraries ####

#source("./code/00_packages.R")
#source("./code/global_variables.R")

#Filtering results
# Function to filter a single data frame
filter_dataframe <- function(df) {
  df %>%
    filter(abs(avg_log2FC) > FC, p_val_adj < 0.05)  %>%
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
}
filtered_list <- map(cluster_list, filter_dataframe)

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
  names(original_gene_list) <- df$gene
  
  # Remove NA values
  gene_list <- na.omit(original_gene_list)
  
  # Sort the gene_list in decreasing order
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}
# Apply the function to each element in the KEGG_input_dataframes list
GSE_gene_lists <- lapply(split_filtered_list, process_gene_list)
GSE_gene_lists <- Filter(Negate(is.null), GSE_gene_lists)




#Use GSE_gene_list as input
convert_to_entrez <- function(gene_set) {
  # Convert gene symbols to ENTREZ IDs
  gene_entrez <- bitr(names(gene_set), fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Mm.eg.db)
  
  # Check if conversion returned any results
  if (nrow(gene_entrez) == 0) {
    warning("No genes were mapped to ENTREZ IDs. Returning NULL.")
    return(NULL)  # Avoid processing an empty set
  }
  
  # Filter out genes that were not mapped
  gene_list_entrez <- gene_set[names(gene_set) %in% gene_entrez$SYMBOL]
  
  # Rename to ENTREZ IDs
  gene_list_entrez <- gene_list_entrez[gene_entrez$SYMBOL]
  names(gene_list_entrez) <- gene_entrez$ENTREZID
  
  # Ensure values are sorted decreasingly
  gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
  
  return(gene_list_entrez)
}

#Run gseKEGG() with Local KEGG Data
#Modify run_gseKEGG() to use local KEGG pathway data.
run_gseKEGG <- function(gene_list_entrez) {
  # Ensure input is valid
  if (is.null(gene_list_entrez) || length(gene_list_entrez) == 0) {
    warning("No valid ENTREZ IDs. Skipping KEGG analysis.")
    return(NULL)
  }
  
  # Debugging: Print first few ENTREZ IDs
  print("Running gseKEGG() with the following gene list:")
  print(gene_list_entrez)
  
  # Run gseKEGG using online KEGG database
  results <- tryCatch({
    gseKEGG(geneList = gene_list_entrez,
            organism = "mmu",
            keyType = "ENTREZID",
            nPerm = 10000,
            minGSSize = 3,
            maxGSSize = 800,
            pvalueCutoff = 0.05,
            pAdjustMethod = "none",
            use_internal_data = FALSE)  # Force online KEGG
  }, error = function(e) {
    warning("gseKEGG() failed: ", e$message)
    return(NULL)
  })
  
  return(results)
}

# Convert gene symbols to ENTREZ IDs
GSE_gene_lists_ENTREZ <- lapply(GSE_gene_lists, convert_to_entrez)
GSE_gene_lists_ENTREZ <- Filter(Negate(is.null), GSE_gene_lists_ENTREZ)

# Run KEGG enrichment using local pathway data
results_list_gseKEGG_ALL <- lapply(GSE_gene_lists_ENTREZ, run_gseKEGG)

#Save results
save(results_list_gseKEGG_ALL, file = paste0(output_path,"RData/results_list_gseKEGG_ALL.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path,"RData/results_list_gseKEGG_ALL.RData"))


## Lolliplots ####

# Make dataframe list from GO results #
results_df_gseKEGG_ALL <- map(results_list_gseKEGG_ALL, ~ .x@result)
# Write each dataframe in the list to a separate sheet in an Excel file
write_xlsx(results_df_gseKEGG_ALL, path = paste0(output_path,"tables/results_df_gseKEGG_ALL.xlsx"))

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
    dplyr::mutate(Description = reorder(Description, enrichmentScore)) %>%
    top_n(15, enrichmentScore) %>%
    ggplot(aes(x = Description,
               y = enrichmentScore,
               colour = p.adjust)) +
    geom_segment(aes(x = Description,
                     xend = Description,
                     y = 0,
                     yend = enrichmentScore)) +
    geom_point(aes(size = setSize-1), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    #facet_wrap(ONTOLOGY ~ ., scale = "free", ncol=1)+
    coord_flip() +
    theme_minimal() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "enrichmentScore") + 
    ggtitle(paste("Plot for", plot_name))
  
  if (!is.null(file_prefix)) {
    tiff_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".tiff")
    pdf_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename = tiff_filename, plot = plot, device = "tiff")
    ggsave(filename = pdf_filename, plot = plot, device = "pdf")
  }
  print(plot)
  return(plot) 
}

# Apply lolliplot function
for (i in seq_along(results_df_gseKEGG_ALL)) {
  lolliplot_result <- lolliplot(i, results_df_gseKEGG_ALL, file_prefix = paste0(output_path,"figures/KEGG_GO/Lolliplot_0",i))
}




#Select manually the desired pathway
'
#KEGG PATHWAY
library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=GSE_gene_lists$names_CLN3_Lux1_vs_WT_UP, pathway.id="dre01200", species = kegg_organism)

print(dme)
# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=GSE_gene_lists$names_CLN3_Lux1_vs_WT_UP, pathway.id="dre00620", species = kegg_organism, kegg.native = T)
knitr::include_graphics("dre00620.pathview.png")
'