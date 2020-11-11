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
seurat=readRDS("../20200410/output/aging-xxx/Mono/model/seurat.rds")

print("Subset seurat")
Idents(seurat)<-seurat$new.Clusters
print(table(Idents(seurat)))
seurat<-RenameIdents(seurat,"1"="CD14","2"="CD14","3"="CD14",
		     "5"="CD14","6"="CD14","7"="CD14","8"="CD14","9"="CD14",
		     "4"="CD16","CD16"="CD16","MED"="MED")

seurat<-subset(seurat,idents=c("CD14","CD16","MED"))
print(table(Idents(seurat)))
seurat$celltype<-Idents(seurat)
seurat$status<-factor(seurat$status,levels=c("YS","AS"))

genes<-c("IFITM1","IFITM2","IFITM3","IFITM5","CCL4", "CCL5", "IL1B", "DUSP2","JUNB", "FOS", "JUN","TGFBR2","IRF7","IRAK1","MAPK1","MYD88","NLRP3")
genes<-str_to_upper(genes)
for(gene in genes){
	p<-VlnPlot(seurat, cols=c("yellow","blue"),features =gene, split.by = "status", group.by = "celltype",pt.size=0.75)
	jpeg(file.path(args$outdir,paste0(gene,".jpeg")),width=1024,height=1024)
	print(p)
	dev.off()
}
