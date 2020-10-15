require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)
 
# peak-bc matrix
mex_dir_path <- "/opt/sample345/outs/filtered_peak_bc_matrix"

mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')
 
features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
 
mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

# tf-bc matrix
mex_dir_path <- "/opt/sample345/outs/filtered_tf_bc_matrix"

mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "motifs.tsv", sep = '/')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')
 
features <- readr::read_tsv(feature_path, col_names = c('feature', 'common_name'))
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
 
mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

