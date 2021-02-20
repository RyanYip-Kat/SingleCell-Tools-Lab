library(Seurat)
library(argparse)
library(stringr)
library(ggplot2)
parser <- ArgumentParser(description="--- Program to create Gex ADT Seurat Object ---")
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--path",
                    type="character",
                    default=NULL,
                    help="the path to of cellranger counts or aggr")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="which column to use subset")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset")

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

parser$add_argument("--suffix",
		    type="character",
		    default=NULL,
		    help="suffix for renaming cells(for seurats merge and show in loupe,eg(1,2,3,4,5 ...)")
		    

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

###########################
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
message("INFO : Loading dataset ...")
path=args$path
if(is.null(args$suffix)){
	suffix="-1" # not replace
}else{
	suffix=paste0("-",args$suffix)
}

countsList=Read10X(file.path(path,"filtered_feature_bc_matrix"))
counts=countsList[["Gene Expression"]]
adt=countsList[["Antibody Capture"]]
cat(sprintf("INFO : raw data size [ %d , %d ]\n",nrow(counts),ncol(counts)))
newBarcode=str_replace(colnames(counts),"-1",suffix)
colnames(counts)=newBarcode
colnames(adt)=newBarcode

message("INFO : Loading graphclust ...")
gphcDF<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"),stringsAsFactors=FALSE)
rownames(gphcDF)=str_replace(gphcDF$Barcode,"-1",suffix)
gphcDF$Cluster=as.character(gphcDF$Cluster)
gphcDFs=subset(gphcDF,select=Cluster)
colnames(gphcDFs)="graphclust"

message("INFO : Loading kmeans clusters  ...")
kmcList=list()
km_clusters=list.files(file.path(path,"analysis/clustering"),"kmeans")
for(i in seq_along(km_clusters)){
	kmc=km_clusters[i]
	kmcp=file.path(path,"analysis/clustering",kmc,"clusters.csv")
	cat(sprintf("INFO : loaing [ %d of %s ]  from [ %s ]\n",i,kmc,kmcp))
	
	kmcDF=read.csv(kmcp,stringsAsFactors=FALSE)
	kmcDF$Cluster=as.character(kmcDF$Cluster)

        rownames(kmcDF)=str_replace(kmcDF$Barcode,"-1",suffix)
	kmcDFs=subset(kmcDF,select=Cluster)
	colnames(kmcDFs)=kmc
	kmcList[[kmc]]=kmcDFs
}
kmcDFL=do.call(cbind,kmcList)
orgMeta=cbind(gphcDFs,kmcDFL)

rownames(orgMeta)=rownames(gphcDFs)
orgMeta$barcode=rownames(gphcDFs)

#############################  filter
if(!is.null(args$column) & !is.null(args$subset)){
	stopifnot(args$column%in%colnames(orgMeta))
	orgMeta=orgMeta[orgMeta[[args$column]]%in%args$subset,,drop=FALSE]
	if(args$invert){
		orgMeta=orgMeta[!orgMeta[[args$column]]%in%args$subset,,drop=FALSE]
	}
}

counts=counts[,orgMeta$barcode]
adt=adt[,orgMeta$barcode]
cat(sprintf("INFO : raw data size [ %d , %d ]\n",nrow(counts),ncol(counts)))

#####################################  clean gex genes,will store in seurat misc slot
gex_genes=rownames(counts)
gex_genes=gex_genes[!str_detect(gex_genes,"^MT-|^RPL|^RPS|^Mt-|^mt-|^Rps|^Rpl")]
gex_genes=gex_genes[!str_detect(gex_genes,"\\.")]

######################################
message("IFNO :  Create Seurat Object")
seurat=CreateSeuratObject(counts=counts,
                       assay = "RNA",
                       project ="GexAdt",
                       names.delim="_",
                       min.cells=0,
                       min.features=0)

message("INFO : AddMetadata from oirginal clusters")
orgMeta=orgMeta[colnames(seurat),]
seurat=AddMetaData(seurat,metadata=orgMeta)

#######################################
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
message("INFO : caculate QC metrics ")
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rpl"]] <- PercentageFeatureSet(seurat, pattern = "^RPL")
seurat[["percent.rps"]] <- PercentageFeatureSet(seurat, pattern = "^RPS-")

message("INFO : add ADT Slot ")
seurat[["ADT"]]=CreateAssayObject(adt)
seurat = NormalizeData(seurat, assay = "ADT", normalization.method = "CLR", margin = 2)
seurat <- ScaleData(seurat, assay = "ADT")
seurat@misc[["gex_genes"]]=gex_genes
seurat@misc[["adt_genes"]]=rownames(adt)

message("INFO : Save Object")
saveRDS(seurat,file.path(args$outdir,"seuratRAW.rds"))

message("INFO : plot QC metric ")
metrics=c("nFeature_RNA","nCount_RNA","percent.mt","percent.rpl","percent.rps")
plot_list<-lapply(metrics,function(metric){
                          p=VlnPlot(seurat, features =metric,pt.size=0.1)
                          return(p)
                           })
MyplotPDF(plot_list, name = "QCplots.pdf", outpath=args$outdir, addDOC = FALSE, width = 8, height = 8)


