#suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(AnnotationDbi))
#suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(Rgraphviz))
#suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stringr))
library(Seurat)


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

parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="the dataset txt file,2 columns with \t")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="")

parser$add_argument("--nfeatures",
                    type="integer",
                    default=10000,
                    help="the dataset  to be used")

parser$add_argument("--vst",
                    action="store_true",
                    default=FALSE)


args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}


DATA=read.table(args$meta,sep=",",header=FALSE,stringsAsFactors=FALSE)
files=as.character(DATA$V1)
meta=as.character(DATA$V2)


seurat_list<-c()
for(i in 1:length(files)){
	x=readRDS(files[i])
	x$celltype=meta[i]
	print(dim(x))
	if(i==1){
		seurat=x
	}else{
		seurat<-merge(x=seurat,y=x,merge.data=TRUE)
	}
	print(dim(seurat))
}

print("### Merge dataset")
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)
genes=rownames(seurat)
point_genes=genes[!str_detect(genes,"\\.")]
seurat=subset(seurat,features=point_genes)

print(dim(seurat))
print(table(seurat$celltype))

if(!is.null(args$column)){
        #Idents(seurat)<-seurat$status
        #seurat<-subset(seurat,idents=args$status)
	metadata=seurat@meta.data
	Idents(seurat)=metadata[[args$column]]
        seurat<-subset(seurat,idents=args$subset)
}

if(args$vst){
	seurat<-FindVariableFeatures(seurat,nfeatures=args$nfeatures)
        seurat<-subset(seurat,features=VariableFeatures(seurat))
}

print("### Get metrix data")
DATA<-GetAssayData(seurat,args$slot)

counts<-as.data.frame(as.matrix(DATA))
genes<-rownames(DATA)
symbols=mapIds(x=org.Hs.eg.db,keys=genes,keytype="SYMBOL",column="ENSEMBL")
Genes<-data.frame(Gene=as.character(symbols))
counts<-cbind(Genes,counts)
counts<-na.omit(counts)

celltypes<-seurat$celltype
cells<-colnames(seurat)
meta.data<-data.frame(Cell=cells,cell_type=celltypes,stringsAsFactors=FALSE)

write.table(counts,file.path(args$outdir,"cell_counts.txt"),sep="\t",row.names = FALSE,quote=F)
write.table(meta.data,file.path(args$outdir,"cell_meta.txt"),sep="\t",row.names = FALSE,quote=F)




