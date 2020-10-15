library(SingleR)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
#parser$add_argument("--outdir",
#                    type="character",
#                    default="output",
#                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the test dataset")

parser$add_argument("--ref",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the reference dataset in singleR")

parser$add_argument("--slot",
                    type="character",
                    default="data",
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
seurat=readRDS(args$seurat)
test_input=GetAssayData(seurat,args$slot)

outdir=dirname(args$seurat)

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

seurat[["SingleR.labels"]] <- predictor$labels
metadata=seurat@meta.data
metadata$barcode=rownames(metadata)
DATA=metadata[,c("barcode","SingleR.labels")]
write.table(DATA,file.path(outdir,"SinglerPrediction.csv"),sep=",",quote=F,row.names=F)

print("### Update seurat")
saveRDS(seurat,args$seurat)
saveRDS(predictor,file.path(outdir,"singleR-predictor.rds"))

