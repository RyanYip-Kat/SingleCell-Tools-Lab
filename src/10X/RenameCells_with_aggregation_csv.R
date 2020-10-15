library(Seurat)
library(argparse)
library(monocle3)
library(stringr)
library(future)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the path of dataset")

parser$add_argument("--aggregation_csv",
                    type="character",
                    default=NULL,
                    help="the path of dataset")



args <- parser$parse_args()
dataset<-args$outdir
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print("### Loading Dataset")
aggregation=read.csv(args$aggregation_csv)
object=readRDS(args$seurat)

############################
print("### Add message from aggregation.csv file")
aggregation$ident=as.character(1:nrow(aggregation))
ident=aggregation$ident
library_id=aggregation$library_id
status=aggregation$status

names(library_id)=ident
names(status)=ident

seurat_sample=dplyr::recode(object$ident,!!!library_id)
seurat_status=dplyr::recode(object$ident,!!!status)
object$sample=seurat_sample
object$status=seurat_status
new_cells=paste0(object$sample,"_",object$status,"_",Cells(object))
object=RenameCells(object,new.names=new_cells)

saveRDS(object,file.path(args$outdir,"seurat.rds"))

