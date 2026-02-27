

# libraries ####

#source("./code/00_packages.R")
# Setting thresholds ####
# Setting threshold for p-value and Fold-Change
p_val <- 0.05
FC <- 0.5

# Global variables
tiff_extension <- ".tiff"
pdf_extension <- ".pdf" #Vectorial format


# Load data ####
import_excel_as_list <- function(output_path) {
  # Obtener los nombres de todas las hojas en el archivo Excel
  sheet_names <- excel_sheets(output_path)
  
  # Filtrar solo las hojas cuyo nombre empieza con "cluster"
  sheet_names_filtered <- sheet_names[str_starts(sheet_names, "cluster")]
  
  # Leer solo las hojas filtradas
  sheet_list <- sheet_names_filtered %>%
    set_names() %>%
    map(~read_excel(output_path, sheet = .x)) %>%
    compact() # Eliminar NULL (en caso de que alguna hoja esté vacía)
  
  return(sheet_list)
}
input_data <- import_excel_sheets(paste0(output_path, "tables/", "cluster_list.xlsx"))

# Filtrado de Resultados ####
filter_dataframe <- function(df) {
  df %>%
    dplyr::filter(abs(avg_log2FC) > FC, p_val_adj < p_val) %>%
    dplyr::mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
}

# Aplicar filtrado a todas las tablas
filtered_list <- purrr::map(input_data, filter_dataframe)
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

# Apply the splitting function to each filtered data frame
split_filtered_list <- map(filtered_list, split_dataframe)
split_filtered_list <- unlist(split_filtered_list, recursive = FALSE)# Remove one level in the list

id_list <- lapply(split_filtered_list, function(df) data.frame(gene = df$gene))

# Strings ####
# load string database
#9606 for Human, 10090 for mouse, 7955 for zebrafish
string_db <- STRINGdb$new(version = "11.5", species = 10090, score_threshold = 200, input_directory="")
class(string_db)


string_function <- function(x) {
  string_example <- string_db$map(x, "gene", removeUnmappedRows = TRUE)
  dimension <- dim(string_example)[1]
  hits <- string_example$STRING_id[1:dimension]
  link <- as.character(string_db$get_link(hits))
  return(link)
}

STRING_plotnames <- lapply(names(id_list), function(name) {
  link <- string_function(id_list[[name]])
  data.frame(Name = name, Link = link, stringsAsFactors = FALSE)
})


### Save String plot links ####
STRING_plotnames_df <- do.call(rbind, STRING_plotnames)
STRING_plotnames_table <- STRING_plotnames_df %>%
  kable("html") %>%
  kable_styling()

STRING_plotnames_table

STRING_plotnames_df %>% 
  write_xlsx(paste0(output_path, "tables/STRING_plotnames_table.xlsx"))

