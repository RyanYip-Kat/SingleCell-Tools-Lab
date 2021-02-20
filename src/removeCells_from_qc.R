library(Seurat)
library(SeuratWrappers)
library(harmony)
library(argparse)
library(stringr)
library(DoubletFinder)
library(future)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--path",
                    type="character",
                    default="",
                    help="the path to of cellranger counts or aggr")


parser$add_argument("--upper",
                    type="double",
                    default=NULL,
		    help="the nFeature_count upper limit")

parser$add_argument("--lower",
                    type="double",
                    default=NULL,
                    help="the nFeature_count lower limit")


parser$add_argument("--mt",
                    type="double",
                    default=NULL,
                    help="the percent.mt")

parser$add_argument("--soupX",
                    action='store_true', default=FALSE)


parser$add_argument("--suffix",
                    type="character",
                    default="1",
                    help="the suffix in the barcode,like : AATCGTACGTC-1")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading dataset")
path=args$path
toc=Read10X(file.path(path,"filtered_feature_bc_matrix"))
rownames(toc)=str_to_upper(rownames(toc))
print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))

if(!is.null(args$suffix)){
	barcode=colnames(toc)
	cells=unlist(lapply(barcode,function(cell){return(str_split(cell,"-")[[1]][1])}))
	cells=paste(cells,args$suffix,sep="-")
	colnames(toc)=cells
}

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

if(!is.null(args$suffix)){
	cells=unlist(lapply(rownames(gp_cluster),function(cell){return(str_split(cell,"-")[[1]][1])}))
	cells=paste(cells,args$suffix,sep="-")
	rownames(Data)=cells
	Data$barcode=cells
}

print("### Create Seurat")
#seurat <- CreateSeuratObject(counts =toc,
#                           min.cells = 1, min.features =1,
#                           meta.data=Data)
seurat <- CreateSeuratObject(counts =toc,
                           min.cells = 0, min.features =0,
			   meta.data=Data)
        
print(head(rownames(seurat)))
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rpl"]] <- PercentageFeatureSet(seurat, pattern = "^RPL")
seurat[["percent.rps"]] <- PercentageFeatureSet(seurat, pattern = "^RPS-")

print(paste0("Before subset,size  of dataset  [ ",nrow(seurat),",",ncol(seurat)," ]"))
print("### Subset ")
seurat <- subset(seurat, subset = nFeature_RNA > args$lower & nFeature_RNA < args$upper & percent.mt < args$mt)
print(paste0("After subset,size  of dataset  [ ",nrow(seurat),",",ncol(seurat)," ]"))

toc=GetAssayData(seurat,"counts")
genes<-rownames(toc)
keep_genes<-genes[!str_detect(genes,"^MT-|^RRL|^RPS")]
keep_genes<-keep_genes[!str_detect(keep_genes,"\\.")]

seurat<-subset(seurat,features=keep_genes)
toc=toc[keep_genes,]

metadata=seurat@meta.data
print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))
if(args$soupX){
        library(SoupX)
        print("### Use SoupX method")
        tod=Read10X(file.path(path,"raw_feature_bc_matrix"))
	rownames(tod)=str_to_upper(rownames(tod))
        tod=tod[keep_genes,]
        print(paste0("Size of tod  [ ",nrow(tod),",",ncol(tod)," ]"))
        sc = SoupChannel(tod, toc)
        sc = setClusters(sc, setNames(metadata[["gp_cluster"]], rownames(metadata)))
        sc = autoEstCont(sc)
        toc = adjustCounts(sc)
        metadata<-metadata[colnames(toc),]
	seurat<-CreateSeuratObject(counts= toc,meta.data=metadata)
}

print("### NormalizeData")
seurat<-NormalizeData(seurat)

print("### Save seurat")
saveRDS(seurat,file.path(args$outdir,"seurat.rds"))
