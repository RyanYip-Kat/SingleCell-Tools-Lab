library(Seurat)
library(argparse)
library(stringr)
library(pheatmap)
library(dplyr)
library(ggplot2)
source("../pheatmap.R")

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


object<-readRDS("../20200416/output/aging-covid/model/seurat.rds")
genes<-rownames(object)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

Idents(object)<-object$UMAP_clusters
object<-subset(object,idents=c("2","5","14"))

status=object$status
x<-str_replace(status,"OHC","ACR")
x<-str_replace(x,"YHC","YCR")
x<-str_replace(x,"OH","AH")

object$status<-x
Idents(object)=object$status
print(table(Idents(object)))


genes<-c("FOS","JUN","IL1B","CCL3","CXCL8","S100A8","S100A12","KLF4","ZFP36","IFITM1","IFITM2","IFITM3","ISG15","IRF1","IRF7","DUSP1","DUSP2","LAMTOR1","IRAK1","PSMB8","PSME1")
my_theme<-theme(axis.title.x = element_blank(),
                axis.text.x = element_text(size=18,angle=90,face = "bold"),
                axis.text.y = element_text(size=18,face="bold"),
                axis.title.y = element_blank(),
                legend.title = element_text(size=18),
                plot.title=element_text(size=22),
                legend.text = element_text(size=18),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())


pdf(file.path(args$outdir,"dotplot1.pdf"),width=8,height=16)
DotPlot(object,features=genes,dot.scale = 10)+coord_flip()+my_theme
dev.off()

pdf(file.path(args$outdir,"dotplot2.pdf"),width=10,height=12)
DotPlot(object,features=genes,dot.scale=10)+coord_flip()+my_theme
dev.off()

pdf(file.path(args$outdir,"dotplot3.pdf"),width=16,height=8)
DotPlot(object,features=genes,dot.scale=10)+coord_flip()+my_theme
dev.off()
