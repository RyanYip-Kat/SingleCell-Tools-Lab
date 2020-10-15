library(argparse)
library(DropletUtils)
library(Seurat)
library(stringr)


parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default="the seurat data")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat=readRDS(args$seurat)
counts=GetAssayData(seurat,"counts")
#counts=t(as.matrix(counts))
metadata=seurat@meta.data


print("# Write data into csv")
write10xCounts(x =GetAssayData(seurat,"counts"), path=file.path(args$outdir,"matrix"),overwrite=TRUE)
#write.table(counts,file.path(args$outdir,"counts.csv"),sep=",",quote=FALSE)
write.table(metadata,file.path(args$outdir,"cell_meta.csv"),sep=",",quote=FALSE)

print(" Get Seurat Embeddings")
for(method in Reductions(seurat)){
	print(paste0("## Get ",method," Embedding"))
	emb=Embeddings(seurat,method)
	emb=as.data.frame(emb)
	barcode=data.frame(Barcode=rownames(emb),row.names=rownames(emb))
	m=cbind(barcode,emb)
	write.table(m,file.path(args$outdir,paste0(method,"_projection.csv")),sep=",",quote=FALSE,row.names=FALSE)
}


