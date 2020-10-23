library(Seurat)
library(argparse)
library(stringr)
library(DropletUtils)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--path",
                    type="character",
                    default=NULL,
                    help="the path of seurat has been created")

parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="the path of seurat has been created")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

counts=Read10X(args$path)
seurat <- CreateSeuratObject(counts =counts)
meta=read.csv(args$meta,stringsAsFactors=FALSE)
colnames(meta)=c("barcode","cluster")
rownames(meta)=meta$barcode
meta=meta[colnames(counts),]

seurat=AddMetaData(seurat,meta)
seurat=NormalizeData(seurat)
print("### Write data")
write.csv(100 * (exp(as.matrix(GetAssayData(object = seurat, assay = "RNA", slot = "data"))) - 1), file.path(args$outdir,"em.csv"), row.names = T)
write.csv(meta,file.path(args$outdir,"metadata.csv"), row.names=F,quote=F)
