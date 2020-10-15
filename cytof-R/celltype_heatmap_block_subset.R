library(Seurat)
library(RColorBrewer)
library(pheatmap)
library(argparse)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat=readRDS("../20200623/recluster/CC20/Total/seurat.rds")
Idents(seurat)=seurat$label.fine

markers=c("CD3","CD4","CD8A","CD28","CD127","CD27","CD45RO","CD25","CD161","TCR-RD","CCR7","CD45RA","CDR6","CD19","CD20","IGD","CXCR5","CD56","CD57","HLA-DR","CD294","CD38","CD11C","CD123","CXCR3","CCR4","CD14","CD16","CD66B")
cts=AverageExpression(seurat,features=markers)[[DefaultAssay(seurat)]]
#cts=as.data.frame(cts)
level=c("CD4 Naive","CD4 Tcm","CD4 Tem","CD4 Treg","CD4 CTL","CD8 Naive","CD8 Tem","CD8 CTL","CD4-CD8-","CD4+CD8+",
	"NKT","GDT","Naive BC","Memory BC","Plasma BC","ABC","NK1","NK2","NK3","cDC","pDC","CD14","CD16","Transitional")
cts=cts[,level]
#cts=as.matrix(cts)
color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
pdf(file.path(args$outdir,"subset_block_heatmap1.pdf"),width=10,height=16)
p=pheatmap(cts,color=color,
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
            scale="column",
            cluster_rows=F,
            cluster_cols=F,fontsize=10)

print(p)
dev.off()


pdf(file.path(args$outdir,"subset_block_heatmap2.pdf"),width=16,height=10)
p=pheatmap(t(cts),color=color,
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
            scale="column",
            cluster_rows=F,
            cluster_cols=F,fontsize=10)

print(p)
dev.off()
