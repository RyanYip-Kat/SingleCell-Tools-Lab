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
seurat=readRDS("../20200418/output/aging-xxx/TC/model/seurat.rds")

print("Subset seurat")
Idents(seurat)<-seurat$new.Clusters
print(table(Idents(seurat)))
seurat<-RenameIdents(seurat,"10"="CD4","11"="CD4","15"="CD4","16"="CD4","17"="CD4","19"="CD4",
		     "4"="CD4","6"="CD4","8"="CD4","9"="CD4","1"="CD8","12"="CD8","13"="CD8",
		     "14"="CD8","18"="CD8","22"="CD8","3"="CD8","5"="CD8","CD8-1"="CD8","CD8-2"="CD8",
		     "2"="NKT","7"="NKT","20"="T-Pro","21"="T-Pro")

seurat<-subset(seurat,idents=c("CD4","CD8","NKT","T-Pro"))
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
