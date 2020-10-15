library(Seurat)
library(stringr)
library(argparse)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--path",
                    type="character",
                    default="",
                    help="the path to of cellranger counts or aggr")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading dataset")
path=args$path
toc=Read10X(file.path(path,"filtered_feature_bc_matrix"))
rownames(toc)=str_to_upper(rownames(toc))
print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))
gp_cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"),stringsAsFactors=FALSE)
barcode=gp_cluster$Barcode
rownames(gp_cluster)<-gp_cluster$Barcode

km_cluster<-read.csv(file.path(path,"analysis/clustering/kmeans_10_clusters/clusters.csv"),stringsAsFactors=FALSE)
rownames(km_cluster)<-km_cluster$Barcode

gp_cluster<-subset(gp_cluster,select=Cluster)
km_cluster<-subset(km_cluster,select=Cluster)
colnames(gp_cluster)="gp_cluster"
colnames(km_cluster)="km_cluster"
print(head(gp_cluster))

Data=cbind(gp_cluster,km_cluster)
Data$gp_cluster=as.character(Data$gp_cluster)
Data$km_cluster=as.character(Data$km_cluster)
rownames(Data)=rownames(gp_cluster)
Data$barcode=rownames(gp_cluster)

ident=unlist(lapply(Data$barcode,function(cell){return(str_split(cell,"-")[[1]][2])}))
Data$ident<-ident

seurat <- CreateSeuratObject(counts =toc,
			   min.cells = 1, min.features =1,
			   meta.data=Data)
print(head(rownames(seurat)))
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rpl"]] <- PercentageFeatureSet(seurat, pattern = "^RPL")
seurat[["percent.rps"]] <- PercentageFeatureSet(seurat, pattern = "^RPS-")
#keep_genes<-keep_genes[!str_detect(keep_genes,"\\.")]
print(head(seurat@meta.data))
# Visualize QC metrics as a violin plot

metrics=c("nFeature_RNA","nCount_RNA","percent.mt","percent.rpl","percent.rps")
plot_list<-lapply(metrics,function(metric){
			  p<-VlnPlot(seurat, features =metric,pt.size=0.1)
			  return(p)
			   })
MyplotPDF(plot_list, name = "dataset_qc_plots.pdf", outpath=args$outdir, addDOC = FALSE, width = 8, height = 8)


