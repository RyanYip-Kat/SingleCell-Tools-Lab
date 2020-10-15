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

parser$add_argument("--path",
                    type="character",
                    default=NULL,
                    help="the path of dataset")

parser$add_argument("--subset",
                    nargs="+",
                    type="double",
                    default=NULL)

parser$add_argument("--column",
                    type="character",
                    default=NULL)

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

args <- parser$parse_args()
dataset<-args$outdir
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print("### Loading Dataset")
path=args$path
counts=Read10X(file.path(path,"filtered_feature_bc_matrix"))
gp_cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"))
km_cluster<-read.csv(file.path(path,"analysis/clustering/kmeans_10_clusters/clusters.csv"))
aggregation=read.csv(file.path(path,"outs/aggregation.csv"))

rownames(gp_cluster)<-gp_cluster$Barcode
rownames(km_cluster)<-km_cluster$Barcode


gp_cluster<-subset(gp_cluster,select=Cluster)
km_cluster<-subset(km_cluster,select=Cluster)
print(paste0("Size of orig counts  [ ",nrow(counts),",",ncol(counts)," ]"))
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

object<-AddMetaData(object,metadata=gp_cluster,col.name="GraphCluster")
object<-AddMetaData(object,metadata=km_cluster,col.name="KmeansCluster")
cells<-colnames(object)

ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident

object$barcode=Cells(object)
############################
print("### Add message from aggregation.csv file")
aggregation$ident=as.character(1:nrow(aggregation))
ident=aggregation$ident
library_id=aggregation$library_id
status=aggregation$status

names(library_id)=ident
names(status)=ident

seurat_sample=dplyr::recode(object$ident,!!!library_id)
seurat_status=dplyr::recode(object$ident,!!!status)
object$sample=seurat_sample
object$status=seurat_status
new_cells=paste0(object$sample,"_",object$status,"_",Cells(object))
object=RenameCells(object,new.names=new_cells)

metadata=object@meta.data
if(!is.null(args$subset)){
        if(is.null(args$column)| (!args$column%in%colnames(metadata))){
                stop("Please provide valid column in metadata")
        }else{
                column=metadata[[args$column]]
                Idents(object)<-column
                object<-subset(object,idents=args$subset,invert=args$invert)
		print(paste0("Size of object after subset [ ",nrow(object,",",ncol(object," ]"))i))
        }
}


#object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT-")
#object[["percent.rpl"]] <- PercentageFeatureSet(object,pattern = "^RPL")
#object[["percent.rps"]] <- PercentageFeatureSet(object,pattern = "^RPS")

object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)

counts<-GetAssayData(object,"counts")
print("#### gene meta data")
pd<-data.frame("GraphCluster"=object$GraphCluster,"KmeansCluster"=object$KmeansCluster,"ident"=object$ident)
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

#object=SCTransform(object,vars.to.regress="nFeature_RNA",verbose = FALSE)
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

saveRDS(object,file.path(args$outdir,"seurat.rds"))

mat<-table(object$UMAP_clusters)
write.table(mat,file.path(args$outdir,"cluster_number.csv"),sep=",",quote=F,row.names=F)


mat<-table(object$ident,object$UMAP_clusters)
mat<-as.matrix(mat)
write.table(mat,file.path(args$outdir,"ident_cluster_number.csv"),sep=",",quote=F)

