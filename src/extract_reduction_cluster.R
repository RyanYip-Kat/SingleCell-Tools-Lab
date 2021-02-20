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
		    nargs="+",
                    default=NULL,
                    type="character")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
slots=c("tsne","umap","motsne","moumap")
for(slot in slots){
	mat<-Embeddings(seurat,slot)
        mat<-as.data.frame(mat)
        mat$Barcode<-rownames(mat)
        mat<-mat[,c(3,1,2)]
        write.table(mat,file.path(args$outdir,paste0(slot,"_projection.csv")),sep=",",quote=FALSE,row.names=FALSE)
}

columns=args$column
metadata=seurat@meta.data
if(is.null(columns)){
	columns=c("seurat_clusters","MtSNE_clusters","MUMAP_clusters","Status")
	columns=columns[columns%in%colnames(metadata)]
}
for(column in columns){
	mat<-metadata[,column,drop=FALSE]
        mat$barcode=rownames(mat)
        mat<-mat[,c(2,1)]
        write.table(mat,file.path(args$outdir,paste0(column,"_table.csv")),sep=",",quote=FALSE,row.names=FALSE)
}




