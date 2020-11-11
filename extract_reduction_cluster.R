library(Seurat)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
		    type="character",
		    default="")

parser$add_argument("--column",
                    default="Clusters",
                    type="character")

parser$add_argument("--reduction",
		    type="character",
		    default="umap")
args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
mat<-Embeddings(seurat,args$reduction)
mat<-as.data.frame(mat)
mat$Barcode<-rownames(mat)
mat<-mat[,c(3,1,2)]
write.table(mat,file.path(args$outdir,paste0(args$reduction,"_projection.csv")),sep=",",quote=FALSE,row.names=FALSE)

mat<-seurat@meta.data[,args$column,drop=FALSE]
mat$barcode=rownames(mat)
mat<-mat[,c(2,1)]
write.table(mat,file.path(args$outdir,paste0(args$reduction,"_clusters.csv")),sep=",",quote=FALSE,row.names=FALSE)


if("new.Clusters"%in%colnames(seurat@meta.data)){
	mat<-table(seurat$new.Clusters)
	write.table(mat,file.path(args$outdir,"new.Cluster_number.csv"),sep=",",quote=F,row.names=F)

	mat<-table(seurat$ident,seurat$new.Clusters)
	mat<-as.matrix(mat)
	write.table(mat,file.path(args$outdir,"ident_new.Cluster_number.csv"),sep=",",quote=F)
}


