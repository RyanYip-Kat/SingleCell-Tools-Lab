library(argparse)
library(ggplot2)
library(stringr)
library(Seurat)
library(Cairo)
library(Matrix)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")

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


parser$add_argument("--outdir",
                    type="character",
                    default="vision_results")

parser$add_argument("--genes",
                    type="character",
                    default="")

parser$add_argument("--groupby",
                    type="character",
                    default="status")

parser$add_argument("--level",
                    nargs="+",
                    type="character",
                    default=NULL)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("Loading data")
seurat_obj<-readRDS(args$seurat)

metadata=seurat_obj@meta.data
meta=metadata[,args$group,drop=F]
rownames(meta)=rownames(metadata)
if(!is.null(args$column1)){
	if(!args$column1%in%colnames(metadata)){
		stop("Invaild columns in metadata!")
	}else{
		s<-paste(args$subset1,collapse=",")
		print(paste0("Subset :",args$column1," with :",s))
		Idents(seurat_obj)=metadata[[args$column1]]
		seurat_obj<-subset(seurat_obj,idents=args$subset1)
	}
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
features<-read.csv(args$genes,sep=",",stringsAsFactors = F,header=FALSE)$V1
seurat_obj<-NormalizeData(seurat_obj,normalization.method = "LogNormalize",verbose = FALSE)

metadata=seurat_obj@meta.data
stopifnot(args$groupby%in%colnames(metadata))
metadata[[args$groupby]]=factor(metadata[[args$groupby]],levels=args$level)
#seurat_obj$status<-factor(seurat_obj$status,levels=c("NPDR","DME"))
DotPlot(seurat_obj,features=features,group.by=args$groupby, dot.scale=20,scale.by="size")+
	theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y =element_text(size=20,face = "bold"),
          axis.text.x = element_text(size=20,face = "bold"),
          axis.title = element_blank())+coord_flip()
filename=file.path(args$outdir,"dotplot.pdf")
ggsave(filename = filename,width = 12,height = 14,limitsize = F,device = cairo_pdf)
