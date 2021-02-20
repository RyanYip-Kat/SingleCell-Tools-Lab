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

parser$add_argument("--seurat_file",
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

#################################
message("INFO : Loading dataset ")
seurat_lists=read.table(args$seurat_file,sep="\t",stringsAsFactors=FALSE)
seurat_samples=seurat_lists$V1
seurat_files=seurat_lists$V2
seurat_status=seurat_lists$V3

seurats=list()
gex_genes_list=list()
for(i in seq_along(seurat_files)){
	cat(sprintf("INFO :  Loading Sample  :  [ %d of %d ]  --- [ %s from %s ]\n",i,length(seurat_files),seurat_samples[i],seurat_files[i]))
	seurat=readRDS(seurat_files[i])
	seurat$Status=seurat_status[i]
	seurat$Sample=seurat_samples[i]
	cat(sprintf("INFO :  The Sample's size --- [ %s # cells %d  ]\n ",seurat_samples[i],ncol(seurat)))
	seurats[[i]]=seurat
	gex_genes_list[[i]]=seurat@misc[["gex_genes"]]
}

object=Reduce(merge,seurats)
gex_genes=Reduce(intersect,gex_genes_list)
object@misc[["gex_genes"]]=gex_genes

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
object = NormalizeData(object, assay = "ADT", normalization.method = "CLR", margin = 2)
object = ScaleData(object, assay = "ADT")
message("INFO : Cluster directly on protein levels")
object <- RunPCA(object, features = rownames(object), reduction.name = "pca_adt", reduction.key = "pca_adt_",
    verbose = FALSE)

# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
object[["rnaClusterID"]] <- Idents(object)

# Now, we rerun tSNE using the PCA only on ADT (protein) levels.
object <- RunTSNE(object, dims = 1:10, reduction = "pca_adt", reduction.key = "adtTSNE_", reduction.name = "tsne_adt")
object <- RunUMAP(object, dims = 1:10, reduction = "pca_adt", reduction.key = "adtUMAP_", reduction.name = "umap_adt")
object <- FindNeighbors(object, features = rownames(object), dims = NULL)
object <- FindClusters(object, resolution = 0.2, graph.name = "ADT_snn")
object[["adtClusterID"]] <- Idents(object)
###########################################
message("INFO : Cluster directly on RNA levels")
DefaultAssay(object)="RNA"
#################################
if(args$addDoublets & !is.null(args$annotation)){
	message("INFO : DetectDoublets")
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
            stop("INFO : Invalid columns in metadata!!!")
        }
        annotations=metadata[[args$annotation]]
        homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- object@meta.data$ClusteringResults
        nExp_poi <- round(0.075*ncol(object))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

        ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
        object <- doubletFinder_v3(object, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
}


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

