library(Seurat)
library(SeuratWrappers)
library(harmony)
library(argparse)
library(monocle3)
library(stringr)
library(DoubletFinder)
library(future)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat_list",
                    type="character",
                    default=NULL,
                    help="the txt file for seurat.rds list")

parser$add_argument("--batch_key",
                    type="character",
                    default=NULL,
                    help="barcode file use for subsetting")

parser$add_argument("--batch_correct",
                    action='store_true', default=FALSE)

parser$add_argument("--addDoublets",
                    action='store_true', default=FALSE)

parser$add_argument("--method",
                    type="character",
                    default="harmony",
		    choices=c("fastmnn","harmony"),
                    help="batch correct method")

parser$add_argument("--annotation",
                    type="character",
                    default="seurat_clusters",
                    help="annotation groups for DoubletFinder")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading dataset")
seurat_lists=read.table(args$seurat_list,stringsAsFactors=FALSE)
seurat_files=seurat_lists$V2
seurat_samples=seurat_lists$V1
seurat_status=seurat_lists$V3

seurats<-list()
for(i in 1:length(seurat_files)){
	print(paste0("#### Loading Sample : ",seurat_samples[i]," from : ",seurat_files[i]))
	seurat=readRDS(seurat_files[i])
	seurat$Status=seurat_status[i]
	seurat$Sample=seurat_samples[i]
	print(paste0("#### The Sample's size ---  #features : ",nrow(seurat), " #cells : ",ncol(seurat)))
	seurats[[i]]=seurat
}

object<-Reduce(merge,seurats)
print("### NormalizeData")
object<-NormalizeData(object)
print("### FindVariableFeatures")
object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures =5000,verbose = FALSE)

print("### Run ScaleData")
object<-ScaleData(object,features=VariableFeatures(object),model.use = "linear",
               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)

print("### Run PCA")
object<-RunPCA(object,npcs = 50,verbose = FALSE)

use_rep="pca"
if(args$batch_correct){
	vars=ifelse(!is.null(args$batch_key),args$batch_key,"Sample")
	stopifnot(vars%in%colnames(object@meta.data))
	if(args$method=="harmony"){
		use_rep="harmony"
		print(paste0("### Run Harmony with :",vars))
		object=RunHarmony(object, group.by.vars =vars,
				  reduction = "pca",dims.use=1:30)
	}else if(args$method=="fastmnn"){
		object_list <-object.list = SplitObject(object, split.by = vars)
		object<-RunFastMNN(object_list,features=3000)
		use_rep="mnn"
	}else{
		stop("Invalid batch correct method !")
	}
}


print("### Run UMAP")
object <- RunUMAP(object, reduction =use_rep, dims = 1:30)
print("### Run TSNE")
object=RunTSNE(object,reduction=use_rep)

print("### Find Clusters")
object <- FindNeighbors(object, reduction =use_rep, dims = 1:30)
object<-FindClusters(object,resolution=0.8)


if(args$addDoublets & !is.null(args$annotation)){
	print("### DetectDoublets")
        ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
        sweep.res.list_kidney <- paramSweep_v3(object, PCs = 1:30, sct = FALSE)
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
        object <- doubletFinder_v3(object, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
}


#################################
print("### use monocel3 method")
counts<-GetAssayData(object,"counts")
#print("#### gene meta data")
#pd<-data.frame("orig.Cluster"=object$orig.Cluster,"ident"=object$ident)
pd=object@meta.data
rownames(pd)<-colnames(counts)
fd <- data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
print("#### new cell data set")
cds<-new_cell_data_set(counts,cell_metadata=pd,gene_metadata=fd)

#print("### preprocess ")
cds<-detect_genes(cds)

#############
#print("### Convert into monocle3 object")
#cds <- as.cell_data_set(object)
cds<-detect_genes(cds)

print("### preprocess cds")
cds <- preprocess_cds(cds,
                      num_dim = 50,
                      method="PCA",
                      norm_method="log")


print("### Align")
vars=ifelse(!is.null(args$batch_key),args$batch_key,"Sample")
stopifnot(vars%in%colnames(object@meta.data))
cds <- align_cds(cds,
                 preprocess_method="PCA",
                 alignment_k=20,
                 residual_model_formula_str="~Size_Factor+num_genes_expressed",
                 alignment_group=vars)


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

#########################
print("### Create ReducedDim from monocle and add clusters")
tSNE_clusters<-clusters(cds,reduction_method="tSNE")
UMAP_clusters<-clusters(cds,reduction_method="UMAP")

tSNE_partitions<-partitions(cds,reduction_method="tSNE")
UMAP_partitions<-partitions(cds,reduction_method="UMAP")

monocle_meta<-data.frame("MtSNE_clusters"=tSNE_clusters,
                         "MUMAP_clusters"=UMAP_clusters,
                         "MtSNE_partitions"=tSNE_partitions,
                         "MUMAP_partitions"=UMAP_partitions,
                         row.names=names(tSNE_clusters))

print("### Add MetaData")
object<-AddMetaData(object,metadata=monocle_meta)

print("### Add reducedDims")
print("#### Add tSNE")
mat<-reducedDims(cds)[["tSNE"]]
colnames(mat)<-paste("MtSNE_",1:ncol(mat),sep = "")
object[["motsne"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "MtSNE_",
                                  assay = DefaultAssay(object))
print("#### Add UMAP")
mat<-reducedDims(cds)[["UMAP"]]
colnames(mat)<-paste("MUMAP_",1:ncol(mat),sep = "")
object[["moumap"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "MUMAP_",
                                  assay = DefaultAssay(object))

print("#### Add PCA")
mat<-reducedDims(cds)[["PCA"]]
colnames(mat)<-paste("MPCA_",1:ncol(mat),sep = "")
object[["mopca"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "MPCA_",
                                  assay = DefaultAssay(object))

print("#### Add Aligned")
mat<-reducedDims(cds)[["Aligned"]]
colnames(mat)<-paste("MAligned_",1:ncol(mat),sep = "")
object[["moaligned"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "MAligned_",
                                  assay = DefaultAssay(object))

print("### Save Seurat object")
saveRDS(object,file.path(args$outdir,"seurat.rds"))


