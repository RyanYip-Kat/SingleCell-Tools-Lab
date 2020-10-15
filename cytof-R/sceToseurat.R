library(scater)
library(SingleCellExperiment)
library(Seurat)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--hm",
                      type="character",
                      default="")

parser$add_argument("--sce",
		    type="character",
		    default="")
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

sce=readRDS(args$sce)
data<-assay(sce,"exprs")

seurat<-CreateSeuratObject(counts=data,
                       assay = "exprs",
                       project ="cytof",
                       names.delim="_",
                       min.cells=0,
                       min.features=1)

meta.data=as.data.frame(colData(sce))
seurat<-AddMetaData(seurat,meta.data)

print("### Normalize Data")
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)

print("### Add harmony")
mat<-readRDS(args$hm)
mat<-as.matrix(mat)
colnames(mat)<-paste("Harmony_",1:ncol(mat),sep = "")
seurat[["harmony"]]<-CreateDimReducObject(embeddings =mat,
                                    key = "Harmony_",
                                    assay = DefaultAssay(seurat))

print("### Run UMAP")
seurat <- RunUMAP(object = seurat, reduction = "harmony", dims = 1:20)
print("### Run TSNE")
seurat <- RunTSNE(object = seurat, reduction = "harmony", dims = 1:20)
print("Run FindNeighbors")
seurat <- FindNeighbors(object = seurat,reduction = "harmony", dims = 1:20)

print("Run Find Clusters")
seurat <- FindClusters(object = seurat,resolution=0.8, verbose =TRUE)


saveRDS(seurat,file.path(args$outdir,"seurat.rds"))

jpeg(file.path(args$outdir,"cluster1.jpeg"),width=1024,height=1024)
p<-DimPlot(seurat,reduction="tsne",label=TRUE,label.size=7.5)
print(p)
dev.off()

jpeg(file.path(args$outdir,"cluster2.jpeg"),width=1024,height=1024)
p<-DimPlot(seurat,reduction="umap",label=TRUE,label.size=7.5)
print(p)
dev.off()
