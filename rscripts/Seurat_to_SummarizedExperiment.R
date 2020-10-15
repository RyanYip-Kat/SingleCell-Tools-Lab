library(Seurat)
library(SummarizedExperiment)
library(stringr)
library(argparse)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="seurat path")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--name",
                    type="character",
                    default="pbmc",
                    help="the dataset  to be used")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
metadata<-seurat@meta.data

logcount=GetAssayData(seurat,"data")
count=GetAssayData(seurat,"counts")

print("Convert into SummarizedExperiment")
se=SummarizedExperiment(assays=list(counts=count,logcounts=logcount),colData=metadata)

print("Saving")
saveRDS(se,file.path(args$outdir,paste0(args$name,"_Summarized-Experiment.rds")))
