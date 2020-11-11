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
seurat=readRDS("../20200411/output/aging-xxx/TC/CD4/model/seurat.rds")

print("Subset seurat")
Idents(seurat)<-seurat$new.Clusters
print(table(Idents(seurat)))
seurat<-RenameIdents(seurat,"1"="Naive","5"="Naive","6"="Naive","9"="Naive",
		     "3"="Tcm","10"="Tcm","4"="Exhausted",
		     "11"="Cytotoxic","12"="Cytotoxic",
		     "2"="Tem","7"="Tem","8"="Tem","Th17"="Tem")

seurat<-subset(seurat,idents=c("Naive","Tcm","Exhausted","Cytotoxic","Tem"))
print(table(Idents(seurat)))
seurat$celltype<-Idents(seurat)
seurat$status<-factor(seurat$status,levels=c("YS","AS"))
genes<-c("CCR7","LEF1","AQP3","CD69","CCR6","CCL5","CXCR6","CTLA4","RORA","RORC","TBX21","GATA3","FOXP3","PDCD1","EOMES","Tnfsf8","HAVCR2","LAG3", "GZMB","GZMK","GNLY")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, cols=c("yellow","blue"),features =gene, split.by = "status", group.by = "celltype",pt.size=0.75)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}

jpeg(file.path(args$outdir,"cluster1.jpeg"),width=1024,height=1024)
DimPlot(seurat,reduction="tsne",label=TRUE,label.size=7.5)
dev.off()

jpeg(file.path(args$outdir,"cluster2.jpeg"),width=1024,height=1024)
DimPlot(seurat,reduction="umap",label=TRUE,label.size=7.5)
dev.off()

Idents(seurat)<-seurat$new.Clusters
seurat<-RenameIdents(seurat,"1"="Naive","5"="Naive","6"="Naive","9"="Naive",
                     "3"="Tcm","10"="Tcm","4"="Exhausted",
                     "11"="Cytotoxic","12"="Cytotoxic",
                     "2"="Th2","7"="Th1","8"="Treg","Th17"="Th17")

jpeg(file.path(args$outdir,"cluster3.jpeg"),width=1024,height=1024)
DimPlot(seurat,reduction="tsne",label=TRUE,label.size=7.5)
dev.off()

jpeg(file.path(args$outdir,"cluster4.jpeg"),width=1024,height=1024)
DimPlot(seurat,reduction="umap",label=TRUE,label.size=7.5)
dev.off()

