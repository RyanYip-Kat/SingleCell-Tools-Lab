library(Seurat)
library(argparse)
library(DoubletFinder)
library(stringr)
library(future)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--path",
                    type="character",
                    default="")

parser$add_argument("--annotation",
                    type="character",
                    default="orig.Cluster")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading Dataset")
counts=Read10X(file.path(args$path,"filtered_feature_bc_matrix"))
cluster<-read.csv(file.path(args$path,"analysis/clustering/graphclust/clusters.csv"))
rownames(cluster)<-cluster$Barcode

cluster<-subset(cluster,select=Cluster)
cell=rownames(cluster)

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
                       min.cells=2,
                       min.features=3)

cluster=cluster[colnames(object),]
object<-AddMetaData(object,metadata=cluster,col.name="orig.Cluster")
cells<-colnames(object)

ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident
print(paste0("number of cells : ",ncol(object)))
print(table(object$orig.Cluster))


object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)

print("### Normalize Data")
object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = FALSE)

print("### Scale Data")
object<-ScaleData(object,features=VariableFeatures(object),model.use = "linear",
               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)

print("### RunPCA")
object <- RunPCA(object, features = VariableFeatures(object = object),verbose=FALSE)
mat=read.csv(file.path(args$path,"analysis/pca/10_components/projection.csv"))
rownames(mat)=mat$Barcode
mat=mat[,-1]
mat=mat[colnames(object),]
mat=as.matrix(mat)
colnames(mat)<-paste("oPCA_",1:ncol(mat),sep = "")
object[["opca"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "oPCA_",
                                  assay = DefaultAssay(object))

print("### RunUMAP")
object <- RunUMAP(object, dims = 1:10)

print("### RunTSNE")
#object <- RunTSNE(object, dims = 1:10)
mat=read.csv(file.path(args$path,"analysis/tsne/2_components/projection.csv"))
rownames(mat)=mat$Barcode
mat=mat[,-1]
mat=mat[colnames(object),]
mat=as.matrix(mat)
colnames(mat)<-paste("tSNE_",1:ncol(mat),sep = "")
object[["tsne"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "tSNE_",
                                  assay = DefaultAssay(object))

print("### Clusters")
object <- FindNeighbors(object, dims = 1:10)
object <- FindClusters(object, resolution = 0.8)


#saveRDS(object,file.path(model.dir,"seurat.rds"))
mat<-table(object$seurat_clusters)
write.table(mat,file.path(args$outdir,"cluster_number.csv"),sep=",",quote=F,row.names=F)


print("### DetectDoublets")
## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(object, PCs = 1:10, sct = FALSE)
#gt.calls <- object@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"]
#sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
#annotations <- object$seurat_clusters
metadata=object@meta.data
if(!args$annotation%in%colnames(metadata)){
        stop("Invalid columns in metadata!!!")
}
annotations=metadata[[args$annotation]]
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- object@meta.data$ClusteringResults
nExp_poi <- round(0.075*ncol(object))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
object <- doubletFinder_v3(object, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

saveRDS(object,file.path(args$outdir,"seurat_doublet.rds"))

print("### Extract DF result")
metadata=object@meta.data
Names=colnames(metadata)
columns=Names[str_detect(Names,"DF.classifications|pANN")]
DATA=metadata[,columns]
colnames(DATA)=c("pANN_score","DF_predict")
DATA$barcode=rownames(DATA)
DATA=DATA[,c(3,1,2)]
write.table(DATA,file.path(args$outdir,"DoubletFinder.csv"),sep=",",quote=FALSE,row.names=FALSE)
