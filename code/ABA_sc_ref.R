library(hdf5r)
library(data.table)
library(Matrix)

# Load metadata
reference_metadata <- fread("./rawdata/metadata.csv")

# Load expression matrix
h5_file <- H5File$new("./rawdata/expression_matrix.hdf5", mode = "r")
h5_file
# Load the counts matrix as sparse
counts_data <- h5_file[["data/counts"]][, ]
counts_sparse <- Matrix(counts_data, sparse = TRUE)


allen_expr <- h5_file[["data"]]
allen_expr <- h5_file[["data/counts"]][, ]

h5_file$close_all()

neuronal_metadata <- metadata[metadata$class_label == "Neuron"]



# Load the counts matrix as sparse
counts_data <- h5_file[["data/counts"]][]
counts_sparse <- Matrix(counts_data, sparse = TRUE)