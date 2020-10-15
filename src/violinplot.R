library(Seurat)
library(ggplot2)
library(argparse)
library(ggpubr)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the path to save result")

#parser$add_argument("--genes",
#		    nargs="+",
#                    type="character",
#                    default=NULL)

parser$add_argument("--groupby",
		    type="character",
                    default="condition")

parser$add_argument("--pt_size",
                    type="double",
                    default=0)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat=readRDS(args$seurat)
#genes=args$genes
genes=rownames(seurat)
metadata=seurat@meta.data
if(!args$groupby%in%colnames(metadata)){
	stop("Invalid group use !!!")
}

my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=25),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

print(paste0("Violin Plot with : ",length(genes)," genes"))
for(gene in genes){
	pdf(file.path(args$outdir,paste0(gene,"_vilonplot.pdf")),width=16,height=12)
	p=VlnPlot(seurat,features=gene,group.by=args$groupby,pt.size=args$pt_size)+
		stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                       method = "t.test")+my_theme
	print(p)
	dev.off()
}
