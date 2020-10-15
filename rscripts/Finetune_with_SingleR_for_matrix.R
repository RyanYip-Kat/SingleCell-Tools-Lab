library(SingleR)
library(SingleCellExperiment)
library(scater)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
#parser$add_argument("--outdir",
#                    type="character",
#                    default="output",
#                    help="the path to save result")

parser$add_argument("--matrix",
                    type="character",
                    default=NULL,
                    help="the test dataset")

parser$add_argument("--ref",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the reference dataset in singleR")

parser$add_argument("--label",
                    type="character",
                    default="label.main",
                    help="the reference dataset'label (label.main,label.fine)")
#parser$add_argument("--invert",
#                    action='store_true', default=FALSE)

args <- parser$parse_args()
#if(!dir.exists(args$outdir)){
#        dir.create(args$outdir,recursive=TRUE)
#}

####################
safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
(cores <- parallel::detectCores())
(bpParam <- safeBPParam(8))

#####################
print("### Loading dataset")
test_input=readRDS(args$matrix)

outdir=dirname(args$matrix)

print("### Loading reference")
if(length(args$ref)>1){
	ref_list<-lapply(args$ref,function(ref){return(readRDS(ref))})
        label_list<-lapply(ref_list,function(ref){return(ref[[args$label]])})
}else{
	ref_list=readRDS(args$ref)
	label_list=ref_list[[args$label]]
}

print("### Fintune")
predictor <- SingleR(test = test_input, 
    ref = ref_list, 
    labels =label_list,
    BPPARAM=bpParam)


print(table(predictor$labels))

SingleR.labels <- predictor$labels
DATA=data.frame("barcode"=colnames(test_input),"SingleR.labels"=SingleR.labels)
saveRDS(DATA,file.path(outdir,"SinglerPrediction.rds"))
write.table(DATA,file.path(outdir,paste0(args$label,"_SinglerPrediction.csv")),sep=",",quote=TRUE,row.names=F)

print("### Update seurat")
saveRDS(predictor,file.path(outdir,paste0(args$label,"_singleR-predictor.rds")))

