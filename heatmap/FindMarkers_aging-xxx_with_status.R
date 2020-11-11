library(Seurat)
library(argparse)
library(monocle3)
library(stringr)
library(pheatmap)
library(dplyr)
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


object<-readRDS("../20200429/output/aging-xxx/Mono/model/mono_delete.rds")
genes<-rownames(object)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
Idents(object)=object$Clusters
object=subset(object,idents=c("1","2","3","5","7","9"))
print(table(Idents(object)))

Idents(object)=object$status
markers <- FindAllMarkers(object, only.pos = FALSE,
                          features = keep_genes,
                          test.use = "wilcox",
                          logfc.threshold = 0.1,
                          min.pct = 0.2,
			  pseudocount.use =0.5 )
print("### Save")
write.table(markers,file.path(args$outdir,"AS_YS_markers.csv"),sep=",",quote=FALSE,row.names=FALSE)

object$status<-factor(object$status,levels=c("YS","AS"))
topn<-markers%>%group_by(cluster)%>%top_n(30,wt=avg_logFC)
genes=topn$gene
pdf(file.path(args$outdir,"heatmap1.pdf"),width=8,height=16)
pheatmap.idents(object,genes,"data")
dev.off()

pdf(file.path(args$outdir,"heatmap2.pdf"),width=10,height=10)
pheatmap.idents(object,genes,"data")
dev.off()

pdf(file.path(args$outdir,"heatmap3.pdf"),width=12,height=16)
pheatmap.idents(object,genes,"data")
dev.off()
