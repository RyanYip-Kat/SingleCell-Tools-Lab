library(Seurat)
library(monocle)
library(argparse)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="the dataset  to be used")
parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="")

parser$add_argument("--nfeatures",
                    type="integer",
                    default=7000,
                    help="the dataset  to be used")

parser$add_argument("--vst",
                    action="store_true",
                    default=FALSE)
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
metadata=seurat@meta.data
if(!is.null(args$subset)){
	if(is.null(args$column)| (!args$column%in%colnames(metadata))){
		stop("Please provide valid column in metadata")
	}else{
		column=metadata[,args$column]
		Idents(seurat)<-column
		seurat<-subset(seurat,idents=args$subset)
	}
}
seurat<-FindVariableFeatures(seurat,nfeatures=args$nfeatures)
if(args$vst){
	seurat<-subset(seurat,features=VariableFeatures(seurat))
}

cds<-as.CellDataSet(seurat)
print("### preprocess")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

print("### filter")
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

print("### DDRtree")
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = 'DDRTree')

cds <- orderCells(cds)
print("### Saving")
saveRDS(cds,file.path(args$outdir,"cds.rds"))
