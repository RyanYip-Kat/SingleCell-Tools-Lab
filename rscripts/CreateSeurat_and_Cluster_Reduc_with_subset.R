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

parser$add_argument("--path",
                    type="character",
                    default="",
                    help="the path to of cellranger counts or aggr")

parser$add_argument("--column",
                    type="character",
                    default="gp_cluster",
                    help="which column to use subset")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset")

parser$add_argument("--barcode",
                    type="character",
                    default=NULL,
                    help="barcode file use for subsetting")

parser$add_argument("--batch_key",
                    type="character",
                    default=NULL,
                    help="barcode file use for subsetting")

parser$add_argument("--batch_correct",
                    action='store_true', default=FALSE)

parser$add_argument("--method",
                    type="character",
                    default="harmony",
		    choices=c("fastmnn","harmony"),
                    help="batch correct method")

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

parser$add_argument("--soupX",
                    action='store_true', default=FALSE)


parser$add_argument("--annotations",
                    type="character",
                    default="seurat_clusters",
                    help="annotations use for doublet finder")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading dataset")
path=args$path
toc=Read10X(file.path(path,"filtered_feature_bc_matrix"))
print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))
gp_cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"),stringsAsFactors=FALSE)
barcode=gp_cluster$Barcode
rownames(gp_cluster)<-gp_cluster$Barcode

km_cluster<-read.csv(file.path(path,"analysis/clustering/kmeans_10_clusters/clusters.csv"),stringsAsFactors=FALSE)
rownames(km_cluster)<-km_cluster$Barcode

gp_cluster<-subset(gp_cluster,select=Cluster)
km_cluster<-subset(km_cluster,select=Cluster)
colnames(gp_cluster)="gp_cluster"
colnames(km_cluster)="km_cluster"
print(head(gp_cluster))

Data=cbind(gp_cluster,km_cluster)
Data$gp_cluster=as.character(Data$gp_cluster)
Data$km_cluster=as.character(Data$km_cluster)
rownames(Data)=rownames(gp_cluster)
Data$barcode=rownames(gp_cluster)

ident=unlist(lapply(Data$barcode,function(cell){return(str_split(cell,"-")[[1]][2])}))
Data$ident<-ident

if(!is.null(args$column) & !is.null(args$subset) & is.null(args$barcode)){
	stopifnot(args$column%in%colnames(Data))
	mat=Data[Data[[args$column]]%in%args$subset,,drop=FALSE]
	if(args$invert){
		mat=Data[!Data[[args$column]]%in%args$subset,,drop=FALSE]
	}

	print(paste0("After subset,Size of cells :",nrow(mat)))
	print(head(mat))
}else if(!is.null(args$barcode)){
	barcode=read.csv(args$barcode,stringsAsFactors=FALSE)
	colnames(barcode)=paste("V",1:ncol(barcode),sep="")
	cell=as.character(barcode$V1)
	mat=Data[Data$barcode%in%cell,,drop=FALSE]
	if(args$invert){
		mat=Data[!Data$barcode%in%cell,,drop=FALSE]
	}
	print(paste0("After subset,Size of cells :",nrow(mat)))
}

toc=toc[,mat$barcode]
print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))
genes<-rownames(toc)
keep_genes<-genes[!str_detect(genes,"^mt-|^Rpl|^Rps")]
keep_genes<-keep_genes[!str_detect(keep_genes,"\\.")]
toc=toc[keep_genes,]
print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))
if(args$soupX){
        library(SoupX)
        print("### Use SoupX method")
        tod=Read10X(file.path(path,"raw_feature_bc_matrix"))
        tod=tod[keep_genes,]
        print(paste0("Size of tod  [ ",nrow(tod),",",ncol(tod)," ]"))
        sc = SoupChannel(tod, toc)
        sc = setClusters(sc, setNames(Data[[args$column]], rownames(Data)))
        sc = autoEstCont(sc)
        toc = adjustCounts(sc)
}

print("### Create Object")
object<-CreateSeuratObject(counts= toc,
                       assay = "RNA",
                       project ="scRNA",
                       names.delim="_",
                       min.cells=0,
                       min.features=0)

print("### AddMetadata from oirginal clusters")
object<-AddMetaData(object,metadata=gp_cluster,col.name="gp_cluster")
object<-AddMetaData(object,metadata=km_cluster,col.name="km_cluster")
cells<-colnames(object)

print("### Get Donor message from barcode")
ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident

############################################
aggregation_dir=file.path(path,"aggregation.csv")
if(file.exists(aggregation_dir)){
	print("### Add meassga from aggregation.csv")
	aggregation=read.csv(aggregation_dir,stringsAsFactors=FALSE)
	extract_cols=c("library_id","molecule_h5")
        add_cols=setdiff(colnames(aggregation),extract_cols)

	aggregation$ident=as.character(1:nrow(aggregation))
	aggregation=aggregation[aggregation$ident%in%unique(Data$ident),]
        library_id=aggregation$library_id

	ident=aggregation$ident
	print(ident)
	names(library_id)=ident
	seurat_sample=dplyr::recode(object$ident,!!!library_id)
	object$Sample=seurat_sample
	for(col in add_cols){
		print(paste0("Add : ",col," into metadata"))
		x<-as.character(aggregation[[col]])
		names(x)=ident
		xx<-dplyr::recode(object$ident,!!!x)
		object@meta.data[[str_to_title(col)]]=xx
	}
}

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
	vars=ifelse(!is.null(args$batch_key),args$batch_key,"ident")
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
vars=ifelse(!is.null(args$batch_key),args$batch_key,"ident")
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


