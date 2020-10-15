library(Seurat)
library(argparse)
library(monocle3)
library(stringr)
library(future)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="the path of seurat has been created")

parser$add_argument("--column",
                    type="character",
                    default="orig.Cluster",
                    help="which column to use subset")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default="",
                    help="the subset")

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading Dataset")
object<-readRDS(args$seurat)
metadata=object@meta.data
Idents(object)=metadata[[args$column]]
object<-subset(object,idents=args$subset,invert=args$invert)
print(paste0("after subset,number of cells : ",ncol(object)))

object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)


counts<-GetAssayData(object,"counts")
print("#### gene meta data")
pd<-data.frame("orig.Cluster"=object$orig.Cluster,"ident"=object$ident)
rownames(pd)<-colnames(counts)
fd <- data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
print("#### new cell data set")
cds<-new_cell_data_set(counts,cell_metadata=pd,gene_metadata=fd)

print("### preprocess ")
cds<-detect_genes(cds)

#############

cds <- preprocess_cds(cds,
                      num_dim = 50,
                      method="PCA",
                      norm_method="log")
                      #residual_model_formula_str="~Size_Factor+num_genes_expressed",
                      #alignment_group="sample")


print("### Align")
cds <- align_cds(cds,
                 preprocess_method="PCA",
                 alignment_k=20,
                 residual_model_formula_str="~Size_Factor+num_genes_expressed",
                 alignment_group="ident")


print("### reduce dimension")
cds <- reduce_dimension(cds,reduction_method="tSNE",preprocess_method="Aligned",cores=8)
cds <- reduce_dimension(cds,reduction_method="UMAP",preprocess_method="Aligned",cores=8)

print("### cluster")
cds<-cluster_cells(cds,
                   reduction_method="UMAP",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

cds<-cluster_cells(cds,
                   reduction_method="tSNE",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

print("### learn graph")
cds<-learn_graph(cds,
                 use_partition=TRUE,
                 close_loop=TRUE)

print("### Save monocle")
saveRDS(cds,file.path(args$outdir,"monocle.rds"))

print(paste0("Size of counts after monocle3 selection [ ",nrow(counts),",",ncol(counts)," ]"))
print("### Normalize Data")
object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = FALSE)

print("### Scale Data")
object<-ScaleData(object,features=VariableFeatures(object),model.use = "linear",
               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)
#object=SCTransform(object,vars.to.regress="nFeature_RNA",verbose = FALSE)
#########################
print("### Create ReducedDim from monocle and add clusters")
tSNE_clusters<-clusters(cds,reduction_method="tSNE")
UMAP_clusters<-clusters(cds,reduction_method="UMAP")

tSNE_partitions<-partitions(cds,reduction_method="tSNE")
UMAP_partitions<-partitions(cds,reduction_method="UMAP")

monocle_meta<-data.frame("tSNE_clusters"=tSNE_clusters,
                         "UMAP_clusters"=UMAP_clusters,
                         "tSNE_partitions"=tSNE_partitions,
                         "UMAP_partitions"=UMAP_partitions,
                         row.names=names(tSNE_clusters))

print("### Add MetaData")
object<-AddMetaData(object,metadata=monocle_meta)

print("### Add reducedDims")
print("#### Add tSNE")
mat<-reducedDims(cds)[["tSNE"]]
colnames(mat)<-paste("tSNE_",1:ncol(mat),sep = "")
object[["tsne"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "tSNE_",
                                  assay = DefaultAssay(object))
print("#### Add UMAP")
mat<-reducedDims(cds)[["UMAP"]]
colnames(mat)<-paste("UMAP_",1:ncol(mat),sep = "")
object[["umap"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "UMAP_",
                                  assay = DefaultAssay(object))

print("#### Add PCA")
mat<-reducedDims(cds)[["PCA"]]
colnames(mat)<-paste("PCA_",1:ncol(mat),sep = "")
object[["pca"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "PCA_",
                                  assay = DefaultAssay(object))

print("#### Add Aligned")
mat<-reducedDims(cds)[["Aligned"]]
colnames(mat)<-paste("Aligned_",1:ncol(mat),sep = "")
object[["aligned"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "Aligned_",
                                  assay = DefaultAssay(object))

saveRDS(object,file.path(args$outdir,"seurat.rds"))

mat<-table(object$UMAP_clusters)
write.table(mat,file.path(args$outdir,"cluster_number.csv"),sep=",",quote=F,row.names=F)


mat<-table(object$ident,object$UMAP_clusters)
mat<-as.matrix(mat)
write.table(mat,file.path(args$outdir,"ident_cluster_number.csv"),sep=",",quote=F)

