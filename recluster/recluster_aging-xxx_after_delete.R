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

args <- parser$parse_args()
dataset<-args$outdir
model.dir<-file.path(dataset,"model")
plot.dir<-file.path(dataset,"plot")
if(!dir.exists(model.dir)){
        dir.create(model.dir,recursive=TRUE)
}

if(!dir.exists(plot.dir)){
        dir.create(plot.dir,recursive=TRUE)
}


print("### Loading Dataset")
path="/home/ye/Work/BioAligment/10X/output/20200407_aging-XXX/outs"
counts=Read10X(file.path(path,"filtered_feature_bc_matrix"))
cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"))
rownames(cluster)<-cluster$Barcode


cluster<-subset(cluster,select=Cluster)
cell=rownames(cluster)

counts<-counts[,cell]
rownames(cluster)=colnames(counts)
print(paste0("Size of counts  [ ",nrow(counts),",",ncol(counts)," ]"))
print("### Create Seurat object")
genes<-rownames(counts)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
counts=counts[keep_genes,]

object<-CreateSeuratObject(counts= counts,
                       assay = "RNA",
                       project ="scRNA",
                       names.delim="_",
                       min.cells=0,
                       min.features=0)

object<-AddMetaData(object,metadata=cluster,col.name="orig.Cluster")
cells<-colnames(object)

ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident
status<-ifelse(ident%in%c(1:5),"AS","YS")
print(table(status))
object$status<-status
print(table(object$orig.Cluster))
Idents(object)<-object$orig.Cluster

object<-subset(object,idents=c("2","14","18","36"),invert=TRUE)
print(paste0("number of cells : ",ncol(object)))
print(table(object$orig.Cluster))

#object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT-")
#object[["percent.rpl"]] <- PercentageFeatureSet(object,pattern = "^RPL")
#object[["percent.rps"]] <- PercentageFeatureSet(object,pattern = "^RPS")

object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)

counts<-GetAssayData(object,"counts")
print("#### gene meta data")
pd<-data.frame("orig.Cluster"=object$orig.Cluster,"status"=object$status,"ident"=object$ident)
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
saveRDS(cds,file.path(model.dir,"monocle.rds"))

print(paste0("Size of counts after monocle3 selection [ ",nrow(counts),",",ncol(counts)," ]"))
print("### Normalize Data")
object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = FALSE)

print("### Scale Data")
object<-ScaleData(object,features=VariableFeatures(object),model.use = "linear",
               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)
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

saveRDS(object,file.path(model.dir,"seurat.rds"))

mat<-table(object$UMAP_clusters)
write.table(mat,file.path(model.dir,"cluster_number.csv"),sep=",",quote=F,row.names=F)


mat<-table(object$ident,object$UMAP_clusters)
mat<-as.matrix(mat)
write.table(mat,file.path(model.dir,"ident_cluster_number.csv"),sep=",",quote=F)

