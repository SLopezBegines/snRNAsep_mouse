

# libraries ####

#source("./code/00_packages.R")
source("./code/global_variables.R")

# Leer datos desde Excel
excel_sheet_reader <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  sheet_data <- lapply(sheets, function(sheet_name) {
    readxl::read_excel(filename, sheet = sheet_name) %>%
      dplyr::select(gene, everything())
  })
  setNames(sheet_data, sheets)
}

# Filter dataframes
#filtration variable will filter for Fold-Change and adjusted p-value or it will take all values without filter
#
filter_dataframe <- function(df, filtration, FC, pval = 0.05) {
  if (filtration == 1) {
    df <- df %>%
      filter(abs(avg_log2FC) > FC, p_val_adj < pval, p_val_adj >0) %>%
      mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
  } else {
    df <- df %>%
      mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
  }
  return(df)
}

# Generate volcano plots function

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

perform_vulcano <- function(df, title) {
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), col = direction, label = gene)) +
    geom_point() +
    theme_DEP1() +
    ggtitle(title) +
    labs(x = "Fold-Change", y = expression(-log[10] * "(p-value)")) +
    geom_vline(xintercept = c(-FC, FC), col = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(p_val), col = "red", linetype = "dashed") +
    scale_colour_manual(values = mycolors) +
    geom_text_repel(max.overlaps = 10) 
}








input_data <- excel_sheet_reader(input_file)

# Filtrar resultados
filtered_list <- map(input_data, ~filter_dataframe(.x, filtration = filtration, FC = FC))

# Generar gráficos volcán
for (i in seq_along(filtered_list)) {
  plot_name <- names(filtered_list)[i] %>%
    gsub("^names_", "", .) %>%
    gsub("\\.name$", "", .) %>%
    gsub("_", " ", .)
  
  p <- perform_vulcano(filtered_list[[i]], paste("Vulcano plot for", plot_name))
  print(p)
  
  # Guardar los gráficos
  filename <- paste0(output_dir, image_number + i, "_vulcano_plot_", plot_name)
  ggsave(paste0(filename, pdf_extension), p, width = 8, height = 6, units = "in", dpi = 300)
  ggsave(paste0(filename, tiff_extension), p, width = 8, height = 6, units = "in", dpi = 300)
}

image_number <- image_number + length(filtered_list)


