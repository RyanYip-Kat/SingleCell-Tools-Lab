library(argparse)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--cluster",
                    default=NULL,
                    type="character")

parser$add_argument("--start_cell",
                    default=NULL,
                    type="character")

parser$add_argument("--reduction",
		    default="character",
                    default="tsne",choices=c("tsne","umap"))

parser$add_argument("--column1",
                    type="character",
                    default=NULL)

parser$add_argument("--subset1",
                    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--column2",
                    type="character",
                    default=NULL)

parser$add_argument("--subset2",
                    nargs="+",
                    type="character",
                    default=NULL)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat_obj<-readRDS(args$seurat)
genes=rownames(seurat_obj)
point_genes=genes[!str_detect(genes,"\\.")]
seurat_obj=subset(seurat_obj,features=point_genes)
print(dim(seurat_obj))
DefaultAssay(seurat_obj)="RNA"


metadata=seurat_obj@meta.data
if(!is.null(args$column1)){
        if(!args$column1%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset1,collapse=",")
                print(paste0("Subset :",args$column1," with :",s))
                Idents(seurat_obj)=metadata[[args$column1]]
                seurat_obj<-subset(seurat_obj,idents=args$subset1)
        }
        metadata=seurat_obj@meta.data
        if(!is.null(args$column2)){
                if(!args$column2%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset2,collapse=",")
                print(paste0("Subset :",args$column2," with :",s))
                Idents(seurat_obj)=metadata[[args$column2]]
                seurat_obj<-subset(seurat_obj,idents=args$subset2)
                }
        }
}

metadata=seurat_obj@meta.data
stopifnot(args$column%in%colnames(metadata))
print("# Run sligslot")
sds <- slingshot(Embeddings(seurat_obj, args$reduction), clusterLabels =metadata[[args$cluster]], 
                 start.clus =args$start_cell, stretch = 0)


print("# Save sligslot")
saveRDS(sds,file.path(args$outdir,"slingshot.rds"))

