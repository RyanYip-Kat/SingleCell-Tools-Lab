library(argparse)
library(stringr)
library(Seurat)
library(VISION)
library(Matrix)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--column",
                    type="character",
                    default=NULL)

parser$add_argument("--subset",
		    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="vision_results")

parser$add_argument("--group",
                    type="character",
                    default="vision_results")

parser$add_argument("--sig",
		    nargs="+",
                    type="character",
                    default="")

parser$add_argument("--nfeatures",
                    type="integer",
                    default=NULL,
                    help="the dataset  to be used")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


seurat_obj<-readRDS(args$seurat)
genes=rownames(seurat_obj)
point_genes=genes[!str_detect(genes,"\\.")]
seurat_obj=subset(seurat_obj,features=point_genes)
DefaultAssay(seurat_obj)="RNA"

metadata=seurat_obj@meta.data
if(!is.null(args$column)){
        if(!args$column%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset,collapse=",")
                print(paste0("Subset :",args$column," with :",s))
                Idents(seurat_obj)=metadata[[args$column]]
                seurat_obj<-subset(seurat_obj,idents=args$subset)
        }
}

if(!is.null(args$nfeatures)){
        seurat_obj<-FindVariableFeatures(seurat_obj,nfeatures=args$nfeatures)
        keep_genes=VariableFeatures(seurat_obj)
}else{
        keep_genes=rownames(seurat_obj)
}

metadata=seurat_obj@meta.data
#meta=metadata[,args$group,drop=F]
print(paste0("Using signatures from :",args$sig))
signatures <-args$sig

print(" Analysis ")
#args<-list()
exprData=GetAssayData(seurat_obj,"counts")[keep_genes,]
totals <- colSums(exprData)
scalefactor <- 10000
exprData <- t(t(exprData) / totals * scalefactor)

#args[["data"]]<-exprData
#args[["meta"]] <- meta
#args[["signatures"]]<-signatures

#vision.obj <- do.call(Vision, args)
print(" Vision ")
vision.obj<-Vision(exprData,signatures=signatures,meta = metadata)

print(" Add Reductions from seurat ")
for(method in Reductions(seurat_obj)){
	name=paste0("Seurat_",method)
	print(name)
	coordinates=Embeddings(seurat_obj,method)
	vision.obj=addProjection(vision.obj,name,coordinates)
}

print(" Analyze running ")
options(mc.cores=10)
vision.obj=analyze(vision.obj)
str(vision.obj)

print(" Save ")
saveRDS(vision.obj, file.path(args$outdir,'vision_results.rds'))

score=getSignatureDifferential(vision.obj)
name=names(score[[args$group]])
score_df=score[[args$group]]
DATA=lapply(name,function(n){
		    df=score_df[[n]]
		    df$pathway=rownames(df)
		    df$cluster=n
		    return(df)})

m=do.call(rbind,DATA)
write.table(m, file.path(args$outdir,'SignatureDifferential.csv'),sep=",",quote=F,row.names=F)


saveRDS(vision.obj@SigScores,file.path(args$outdir,'SigScores.rds'))
#print("View result")
#viewResults(vision.obj,host="10.100.110.103",port="1234",browser=FALSE)



