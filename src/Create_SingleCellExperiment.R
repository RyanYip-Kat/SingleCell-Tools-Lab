library(cowplot)
library(scater)
library(SingleCellExperiment)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--count",
                    type="character",
                    default="")

parser$add_argument("--tsne",
                      type="character",
                      default=NULL)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
(cores <- parallel::detectCores())
(bpParam <- safeBPParam(cores))

counts<-readRDS(args$count)
counts<-t(counts)
print(dim(counts))
sce=SingleCellExperiment(list(exprs=counts))
Name<-colnames(counts)
#ident=unlist(lapply(Name,function(name){return(str_split(name,"_")[[1]][1])}))
ident=unlist(lapply(Name,function(name){return(str_split(name,"_\\d+")[[1]][1])}))
print(table(ident))
experiment_info<-as.data.frame(as.matrix(table(ident)))
print(head(experiment_info))
colnames(experiment_info)<-"n_cells"
experiment_info$sample_id<-rownames(experiment_info)
experiment_info$patient_id<-rownames(experiment_info)
experiment_info$condition<-unlist(lapply(rownames(experiment_info),function(x){return(str_split(x,"-")[[1]][1])}))
metadata(sce)$experiment_info=experiment_info
print(sce)

colData(sce)$sample_id=ident
colData(sce)$patient_id=ident
colData(sce)$condition<-unlist(lapply(Name,function(x){return(str_split(x,"-")[[1]][1])}))
print(table(colData(sce)$condition))
rowData(sce)$antigen=rownames(sce)

if(!is.null(args$tsne)){
	   DATA<-readRDS(args$tsne)
           DATA<-DATA[,c(1,2)]
           colnames(DATA)=c("X1","X2")
           DATA<-DATA[colnames(sce),]

           print("###  Update tSNE")
           reducedDim(sce,"TSNE")=DATA
           colData(sce)$metacluster=DATA$cluster
}

saveRDS(sce,file.path(args$outdir,"sce.rds"))
