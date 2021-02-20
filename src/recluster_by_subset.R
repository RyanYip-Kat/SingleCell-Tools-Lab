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

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="seurat object rds ")

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


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

#################################
print("### Loading Dataset")
object<-readRDS(args$seurat)
if(!is.null(args$column) & !is.null(args$subset)){
        metadata=object@meta.data
        Idents(object)=metadata[[args$column]]
        object<-subset(object,idents=args$subset,invert=args$invert)
        print(paste0("after subset,number of cells : ",ncol(object)))
}

DefaultAssay(object)="RNA"
message("INFO :  NormalizeData")
object=NormalizeData(object)

message("INFO :  FindVariableFeatures")
object =  FindVariableFeatures(object, selection.method = "vst",
                            nfeatures =2000,verbose = FALSE)


message("INFO : Cell-Cycle Scoring and Regression ")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

object = CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
message("INFO :  Run ScaleData")
object=ScaleData(object,features=VariableFeatures(object),model.use = "linear",
               vars.to.regress = c("S.Score", "G2M.Score"),verbose =FALSE)
# object = SCTransform(object, vars.to.regress =c("S.Score", "G2M.Score"), verbose = FALSE)


message("INFO :  Run PCA")
object = RunPCA(object,features=gex_genes,npcs = 50,verbose = FALSE)

use_rep="pca"
if(args$batch_correct){
	vars=ifelse(!is.null(args$batch_key),args$batch_key,"Sample")
	stopifnot(vars%in%colnames(object@meta.data))
	if(args$method=="harmony"){
		use_rep="harmony"
		cat(sprintf("INFO : Run Harmony with : [ %s ]\n",vars))
		object=RunHarmony(object, group.by.vars =vars,
				  reduction = "pca",dims.use=1:30)
	}else if(args$method=="fastmnn"){
		object_list <-object.list = SplitObject(object, split.by = vars)
		object<-RunFastMNN(object_list,features=3000)
		use_rep="mnn"
	}else{
		stop("INFO : Invalid batch correct method !")
	}
}


message("INFO : Run UMAP ")
object = RunUMAP(object, reduction =use_rep, dims = 1:20)

message("INFO : Run TSNE")
object = RunTSNE(object,reduction=use_rep,dims=1:20)

message("INFO: Find Clusters")
object = FindNeighbors(object, reduction =use_rep, dims = 1:20)
object = FindClusters(object,resolution=0.8,algorithm=4)

DefaultAssay(object)="ADT"
object = NormalizeData(object, assay = "ADT", normalization.method = "CLR")
object = ScaleData(object, assay = "ADT")

DefaultAssay(object)="RNA"
#################################
message("INFO :  use monocel3 method")
counts=GetAssayData(object,slot="counts",assay="RNA")
pd=object@meta.data
rownames(pd)<-colnames(counts)
fd <- data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
message("INFO : create  new cell data set")
cds<-new_cell_data_set(counts,cell_metadata=pd,gene_metadata=fd)
cds<-detect_genes(cds)

###################################
#cds <- as.cell_data_set(object)

message("INFO :  preprocess cds")
cds = preprocess_cds(cds,
                      num_dim = 50,
                      method="PCA",
                      norm_method="log")


message("INFO :  Align")
vars=ifelse(!is.null(args$batch_key),args$batch_key,"Sample")
stopifnot(vars%in%colnames(object@meta.data))
cds = align_cds(cds,
                 preprocess_method="PCA",
                 alignment_k=20,
                 residual_model_formula_str="~Size_Factor+num_genes_expressed",
                 alignment_group=vars)


message("INFO :  reduce dimension")
cds = reduce_dimension(cds,reduction_method="tSNE",preprocess_method="Aligned",cores=8)
cds = reduce_dimension(cds,reduction_method="UMAP",preprocess_method="Aligned",cores=8)

message("INFO : cluster ")
cds=cluster_cells(cds,
                   reduction_method="UMAP",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

cds=cluster_cells(cds,
                   reduction_method="tSNE",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

message("INFO : learn graph")
cds=learn_graph(cds,
                 use_partition=TRUE,
                 close_loop=TRUE)

message("INFO : Save monocle")
saveRDS(cds,file.path(args$outdir,"monocle.rds"))

#########################
message("INFO : Create ReducedDim from monocle and add clusters")
tSNE_clusters=clusters(cds,reduction_method="tSNE")
UMAP_clusters=clusters(cds,reduction_method="UMAP")

tSNE_partitions=partitions(cds,reduction_method="tSNE")
UMAP_partitions=partitions(cds,reduction_method="UMAP")

monocle_meta=data.frame("MTSNE_clusters"=tSNE_clusters,
                         "MUMAP_clusters"=UMAP_clusters,
                         "MTSNE_partitions"=tSNE_partitions,
                         "MUMAP_partitions"=UMAP_partitions,
                         row.names=names(tSNE_clusters))

message("INFO :  Add MetaData")
object=AddMetaData(object,metadata=monocle_meta)

message("INFO : Add reducedDims")
message("INFO :  Add tSNE")
mat=reducedDims(cds)[["tSNE"]]
colnames(mat)=paste("MTSNE_",1:ncol(mat),sep = "")
object[["motsne"]]=CreateDimReducObject(embeddings =mat,
                                  key = "MTSNE_",
                                  assay = DefaultAssay(object))

message("INFO : Add UMAP")
mat=reducedDims(cds)[["UMAP"]]
colnames(mat)=paste("MUMAP_",1:ncol(mat),sep = "")
object[["moumap"]]=CreateDimReducObject(embeddings =mat,
                                  key = "MUMAP_",
                                  assay = DefaultAssay(object))

message("INFO : Add PCA")
mat=reducedDims(cds)[["PCA"]]
colnames(mat)=paste("MPCA_",1:ncol(mat),sep = "")
object[["mopca"]]=CreateDimReducObject(embeddings =mat,
                                  key = "MPCA_",
                                  assay = DefaultAssay(object))

message("INFO : Add Aligned")
mat=reducedDims(cds)[["Aligned"]]
colnames(mat)=paste("MAligned_",1:ncol(mat),sep = "")
object[["moaligned"]]=CreateDimReducObject(embeddings =mat,
                                  key = "MAligned_",
                                  assay = DefaultAssay(object))

message("INFO : Save Seurat object")
saveRDS(object,file.path(args$outdir,"seurat.rds"))

message("INFO : Done!")

