library(ggplot2)
library(RColorBrewer)
library(argparse)
library(Seurat)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
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

# cols = c("#191970","yellow","#FF3030")
seurat<-readRDS(args$seurat)
print(table(seurat$celltype))
if(!is.null(args$level)){
	seurat$celltype=factor(seurat$celltype,levels=args$level)
}

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
n=length(unique(seurat$celltype))
pdf(file.path(args$outdir,"cluster.pdf"),width=12,height=12)
DimPlot(seurat,reduction="tsne",group.by="celltype",pt.size=args$pt_size,cols=getPalette(n))+my_theme
	#scale_fill_manual(values = getPalette(n))+my_theme
dev.off()

pdf(file.path(args$outdir,"condition_cluster.pdf"),width=24,height=8)
DimPlot(seurat,reduction="tsne",group.by="celltype",split.by="condition",pt.size=args$pt_size,cols=getPalette(n))+my_theme
        #scale_fill_manual(values = getPalette(n))+my_theme
dev.off()

mat=as.matrix(table(seurat$sample_id,seurat$celltype))
write.table(mat,file.path(args$outdir,"ident_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$sample_id,seurat$seurat_clusters))
write.table(mat,file.path(args$outdir,"ident_clusters_number.csv"),sep=",",quote=F)


mat=as.matrix(table(seurat$condition,seurat$celltype))
write.table(mat,file.path(args$outdir,"condition_celltype_number.csv"),sep=",",quote=F)

mat=as.matrix(table(seurat$celltype))
write.table(mat,file.path(args$outdir,"celltype_number.csv"),sep=",",quote=F)

