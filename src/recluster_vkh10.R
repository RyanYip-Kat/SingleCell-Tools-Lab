library(Seurat)
library(argparse)
library(stringr)
library(ggplot2)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="TC")

parser$add_argument("--resolution",
                    type="double",
                    default=1.0)
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=15),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())


seurat<-readRDS("../20200712/output/VKH-12/Total/seurat.rds")
seurat_x=readRDS(args$seurat)
Idents(seurat)=seurat$seurat_clusters
print(dim(seurat))
print(table(Idents(seurat)))
seurat<-subset(seurat,cells=Cells(seurat_x))
print(dim(seurat))

print("### Normalize Data")
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)

print("### Run TSNE")
seurat <- RunTSNE(object = seurat, reduction = "harmony", dims = 1:20)
print("Run FindNeighbors")

print("Run FindNeighbors")
seurat <- FindNeighbors(object = seurat,reduction = "harmony", dims = 1:20)

print("Run Find Clusters")
seurat <- FindClusters(object = seurat,resolution=args$resolution, verbose =TRUE)


saveRDS(seurat,file.path(args$outdir,"seurat.rds"))

jpeg(file.path(args$outdir,"cluster.jpeg"),width=1024,height=1024)
p<-DimPlot(seurat,reduction="tsne",pt.size=2.0,label=T,label.size=7.5)+my_theme
print(p)
dev.off()

mat=table(seurat$seurat_clusters)
write.table(mat,file.path(args$outdir,"ident_clusters.csv"),sep=",",quote=F)

