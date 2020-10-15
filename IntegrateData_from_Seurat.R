library(Seurat)
library(argparse)
library(stringr)
library(DropletUtils)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the path of seurat has been created")

parser$add_argument("--split_by",
                    type="character",
                    default="time_point",
                    help="the column used for split by ")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
seurat_obj<-readRDS(args$seurat)

print("### Split object")
seurat.list <- SplitObject(seurat_obj, split.by =args$split_by)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
seurat.combined <- IntegrateData(anchorset = seurat.anchors, dims = 1:20)
DefaultAssay(seurat.combined) <- "integrated"


#metadata=seurat_obj@meta.data
#meta=metadata[,c("item","label")]
#print(table(meta$label))
#write.table(meta,file.path(args$outdir,"label.csv"),sep=",",quote=FALSE,row.names=FALSE)
write10xCounts(x =GetAssayData(seurat_obj,assay="integrated",slot="data"), path=file.path(args$outdir,"matrixi_Integrate"),overwrite=TRUE)





