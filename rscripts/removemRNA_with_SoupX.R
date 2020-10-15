library(Seurat)
library(monocle3)
library(argparse)
library(stringr)
library(SoupX)
library(ggplot2)
library(Cairo)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

#parser$add_argument("--barcode",
#                    type="character",
#                    default="")

parser$add_argument("--subset",
                    nargs="+",
                    type="double",
                    default=NULL)

parser$add_argument("--column",
                    type="character",
                    default=NULL)

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

parser$add_argument("--path",
                    type="character",
                    default="cellranger aggr or counts 'outs' path")
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("### Loading Dataset from : ",args$path))
sc=load10X(args$path)
metaData=sc$metaData
metaData=subset(metaData,select=c(nUMIs,clusters,clustersFine))
colnames(metaData)=c("nUMIs","GraphCluster","KmeansCluster")

sc = autoEstCont(sc)
out = adjustCounts(sc)
object<-CreateSeuratObject(counts= counts,
                       assay = "RNA",
                       project ="scRNA",
                       names.delim="_",
                       min.cells=0,
                       min.features=0)
object<-AddMetaData(object,metadata=metaData)
object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT-")
object[["percent.rpl"]] <- PercentageFeatureSet(object,pattern = "^RPL")
object[["percent.rps"]] <- PercentageFeatureSet(object,pattern = "^RPS")

print("### Get QC Metric")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","nUMIs"), ncol = 2)
ggsave(filename =file.path(args$outdir,"QCs.pdf"),width = 12,height =16 ,limitsize = F,device = cairo_pdf)
dev.off()

genes<-rownames(object)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
object=subset(object,features=keep_genes)

cells<-colnames(object)

ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident

metadata=object@meta.data
if(!is.null(args$subset)){
        if(is.null(args$column)| (!args$column%in%colnames(metadata))){
                stop("Please provide valid column in metadata")
        }else{
                column=metadata[[args$column]]
                Idents(object)<-column
                object<-subset(object,idents=args$subset,invert=args$invert)
                print(paste0("Size of object after subset [ ",nrow(object,",",ncol(object," ]"))i))
        }
}

