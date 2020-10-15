library(Seurat)
library(ggplot2)
library(argparse)
library(stringr)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--column",
                    type="character",
                    default=NULL)

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="vision_results")

parser$add_argument("--groupby",
                    type="character",
                    default="label_fine")

parser$add_argument("--features",
                    nargs="+",
                    type="character",
                    default="")

parser$add_argument("--pt_size",
                    type="double",
                    default="1.0")

parser$add_argument("--level",
                    nargs="+",
                    type="character",
                    default=NULL)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=25),
                  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())


seurat_obj<-readRDS(args$seurat)
metadata=seurat@meta.data
if(!is.null(args$level)){
	metadata[[args$groupby]]=factor(metadata[[args$groupby]],levels=args$level)
	seurat_obj@meta.data[[args$groupby]]=metadata[[args$groupby]]
}

if(!is.null(args$column)){
        if(!args$column%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset,collapse=",")
                print(paste0("Subset :",args$column," with :",s))
                Idents(seurat_obj)=metadata[[args$column]]
                seurat_obj<-subset(seurat_obj,idents=args$subset)
        }
}

for(feature in args$features){
	filename=paste0(feature,"_violin.pdf")
	pdf(file.path(args$outdir,filename),width=16,height=12)
	p<-VlnPlot(seurat_obj,features=feature,group.by=args$groupby,pt.size=args$pt_size)+my_theme+NoGrid()
        print(p)
	dev.off()
}
