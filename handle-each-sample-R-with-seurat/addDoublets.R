library(Seurat)
library(SeuratWrappers)
library(harmony)
library(argparse)
library(stringr)
library(DoubletFinder)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help=" seurat.rds ")

parser$add_argument("--annotation",
                    type="character",
                    default="seurat_clusters",
                    help="annotation groups for DoubletFinder")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading dataset")
object=readRDS(args$seurat)

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
print("### Finish doubletFinder")

print("### Save object")
saveRDS(object,file.path(args$outdir,"seurat_doublet.rds"))


