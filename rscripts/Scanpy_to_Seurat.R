library(argparse)
library(stringr)
library(Seurat)
library(Matrix)
library(mindr)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(reticulate)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--adata",
                    type="character",
                    default="")


parser$add_argument("--outdir",
                    type="character",
                    default="results")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("# Convert into R object")
options(stringsAsFactors = FALSE)
ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad(args$adata)
if(py_has_attr(ad_object,"raw")){
	raw=ad_object$raw$to_adata()
}else{
	raw=ad_object
}

# access normalized data matrix
data.input <- t(py_to_r(raw$X))
rownames(data.input) <- rownames(py_to_r(raw$var))
colnames(data.input) <- rownames(py_to_r(raw$obs))
# access meta data
meta.data <- py_to_r(ad_object$obs)

print("# Create Seurat")
seurat=CreateSeuratObject(data.input,meta.data=meta.data)

print("# Normalize")
seurat=NormalizeData(seurat,normalization.method = "LogNormalize",verbose=TRUE)

print("# Add Embeddings")
obsm=py_to_r(ad_object$obsm)
m=obsm$get("X_tsne")
colnames(m)<-paste("tSNE_",1:ncol(m),sep = "")
seurat[["tsne"]]<-CreateDimReducObject(embeddings =m,
                                  key = "tSNE_",
                                  assay = DefaultAssay(seurat))

m=obsm$get("X_umap")
colnames(m)<-paste("UMAP_",1:ncol(m),sep = "")
seurat[["umap"]]<-CreateDimReducObject(embeddings =m,
                                  key = "UMAP_",
                                  assay = DefaultAssay(seurat))


print("# Save ")
saveRDS(seurat,file.path(args$outdir,"adata.seurat.rds"))
