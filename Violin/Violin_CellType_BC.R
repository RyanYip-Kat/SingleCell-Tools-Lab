suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

print("*** Configure Parameters ***")
parser <- ArgumentParser(description='Process some tasks')

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")


args<-parser$parse_args()
if(!dir.exists(args$outdir)){
	dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat=readRDS("../20200409/output/aging-xxx/BC/model/seurat.rds")

print("Subset seurat")
Idents(seurat)<-seurat$new.Clusters
print(table(Idents(seurat)))
seurat<-RenameIdents(seurat,"1"="Naive","4"="Naive","5"="Naive",
		     "2"="Memory","3"="Memory","6"="ASC","7"="ASC","ASC"="ASC",
		     "ABC"="ABC")

seurat<-subset(seurat,idents=c("Naive","Memory","ASC","ABC"))
print(table(Idents(seurat)))
seurat$celltype<-Idents(seurat)
seurat$status<-factor(seurat$status,levels=c("YS","AS"))

genes<-c("IL4R","TCL1A","CD27","CD38","MZB1","IGHG1","IGHA1","FAS","SSR4","CXCR4","ITGAX","TNFRSF13B","FCRL3","FAS","ZEB2")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, cols=c("yellow","blue"),features =gene, split.by = "status", group.by = "celltype",pt.size=0.75)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}
