suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(pathview))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Rgraphviz))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Seurat))


set.seed(7777)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--slot",
                    type="character",
                    default="data",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="the dataset txt file,2 columns with \t")

parser$add_argument("--status",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="")
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}


DATA=read.table(args$meta,sep="\t",header=FALSE,stringsAsFactors=FALSE)
files=as.character(DATA$V1)
meta=as.character(DATA$V2)


seurat_list<-c()
for(i in 1:length(files)){
	seurat=readRDS(files[i])
	seurat$celltype=meta[i]
	seurat_list<-c(seurat_list,seurat)
}

print("### Merge dataset")
seurat<-merge(x=seurat_list[1],y=seurat_list[2:length(seurat_list)],merge.data=TRUE)
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)

print(table(seurat$celltype))

if(!is.null(args$status)){
        Idents(seurat)<-seurat$status
        seurat<-subset(seurat,idents=args$status)
}

print("### Get metrix data")
DATA<-GetAssayData(seurat,args$slot)

counts<-as.data.frame(as.matrix(DATA))
genes<-rownames(DATA)
symbols=mapIds(x=org.Hs.eg.db,keys=genes,keytype="SYMBOL",column="ENSEMBL")
Genes<-data.frame(Gene=as.character(symbols))
counts<-cbind(Genes,counts)
counts<-na.omit(counts)

celltypes<-as.character(Idents(seurat))
cells<-colnames(seurat)
meta.data<-data.frame(Cell=cells,cell_type=celltypes,stringsAsFactors=FALSE)

write.table(counts,"cat_counts.txt",sep="\t",row.names = FALSE,quote=F)
write.table(meta.data,"cat_meta.txt",sep="\t",row.names = FALSE,quote=F)




