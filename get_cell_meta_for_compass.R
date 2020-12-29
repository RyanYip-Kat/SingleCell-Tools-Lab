library(Seurat)
library(stringr)
library(argparse)


set.seed(7777)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="seurat rds file")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")


parser$add_argument("--column1",
                    type="character",
                    default=NULL,
                    help="the column1 as subset1 column")

parser$add_argument("--subset1",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset1")

parser$add_argument("--column2",
                    type="character",
                    default=NULL,
                    help="the column2 as subset2 column")

parser$add_argument("--subset2",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset2")


parser$add_argument("--nFeatures",
                    type="integer",
                    default=7000,
                    help="the number of vst genes to use")

parser$add_argument("--vst",
                    action="store_true",
                    default=FALSE)


parser$add_argument("--sample_n",
		    type="integer",
                    default=NULL)

parser$add_argument("--sample_col",
                    type="character",
                    default=NULL,
                    help="sample column")


args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}


message("INFO : loading dataset ...")
seurat=readRDS(args$seurat)

###########################
if(!is.null(args$column1) & !is.null(args$subset1)){
	metadata=seurat@meta.data
	stopifnot(args$column1%in%colnames(metadata))
	column1=metadata[[args$column1]]
	stopifnot(sum(args$subset1%in%unique(column1))==length(args$subset1))
        Idents(seurat)<-column1
        seurat<-subset(seurat,idents=args$subset1)

	if(!is.null(args$column2) & !is.null(args$subset2)){
		   metadata=seurat@meta.data
                   stopifnot(args$column2%in%colnames(metadata))
                   column2=metadata[[args$column2]]
                   stopifnot(sum(args$subset2%in%unique(column2))==length(args$subset2))
                   Idents(seurat)<-column2
                   seurat<-subset(seurat,idents=args$subset2)
	}
}

###########################
if(args$vst){
	nfeatures=ifelse(!is.null(args$nFeatures),args$nFeatures,5000)
	seurat<-FindVariableFeatures(seurat,nfeatures=nfeatures)
        seurat<-subset(seurat,features=VariableFeatures(seurat))
}


############################
if(!is.null(args$sample_n) & !is.null(args$sample_col)){
	metadata=seurat@meta.data
        metadata$cells=rownames(metadata)
	stopifnot(args$sample_col%in%colnames(metadata))
	N=args$sample_n
	d_count=table(metadata[[args$sample_col]])
	d_min=min(d_count)
	if(N> d_min){
		N=d_min
	}
	column=metadata[[args$sample_col]]
	cells=unlist(lapply(unique(column),function(x){
			     cell=metadata$cells[which(column==x)]
			     cell=sample(cell,size=N,replace=FALSE)}))
	seurat=subset(seurat,cells=cells)
}

###############################
message("INFO : NormalizeData ...")
seurat<-NormalizeData(seurat)

###############################
metadata=seurat@meta.data
metadata$cells=rownames(metadata)
cells=paste("SRR",1:nrow(metadata),sep="")
metadata$cell_id=cells

###########################
message("INFO : write cell metadata ...")
filename=file.path(args$outdir,"cell_metadata.csv")
write.table(metadata,filename,sep=",",row.names=FALSE,quote=FALSE)

message("INFO : write gene metadata ...")
filename=file.path(args$outdir,"gene_metadata.csv")
gene_meta=data.frame("symbol"=rownames(seurat))
write.table(gene_meta,filename,sep=",",row.names=FALSE,quote=FALSE)

message("INFO : write matrix ...")
DF=GetAssayData(seurat,slot="data")
colnames(DF)=cells
filename=file.path(args$outdir,"linear_gene_expression_matrix.tsv")
write.table(DF,filename,sep="\t",quote=FALSE)

message("INFO : Done!")

