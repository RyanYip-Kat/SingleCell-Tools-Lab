library(Seurat)
library(argparse)
library(DoubletFinder)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--annotation",
                    type="character",
                    default="seurat_clusters")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading Dataset")
object=readRDS(args$seurat)

print("### DetectDoublets")
## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(object, PCs = 1:10, sct = FALSE,num.cores=4)
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
object <- doubletFinder_v3(object, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

saveRDS(object,file.path(args$outdir,"seurat_doublet.rds"))

print("### Extract DF result")
metadata=seurat@meta.data
Names=colnames(metadata)
columns=Names[str_detect(Names,"DF.classifications|pANN")]
DATA=metadata[,columns]
colnames(DATA)=c("pANN_score","DF_predict")
DATA$barcode=rownames(DATA)
DATA=DATA[,c(3,1,2)]
write.table(DATA,file.path(args$outdir,"DoubletFinder.csv"),sep=",",quote=FALSE,row.names=FALSE)
