library(reticulate)
ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad("scanpy_object.h5ad")
# access normalized data matrix
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
# access meta data
meta.data <- py_to_r(ad_object$obs)
identity <- meta.data
