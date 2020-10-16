library(limma)
library(GSVA)
library(GSEABase)
library(limma)
library(argparse)
library(stringr)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default="",
                    help="expression data,csv,tsv")


parser$add_argument("--metadata",
                    type="character",
                    default="",
                    help="cell metadata")


parser$add_argument("--column",
                    type="character",
                    default="status",
                    help="column in metadata for group")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")


parser$add_argument("--gmt",
                    type="character",
                    default="",
                    help="gmt file,eg,c2.cp.kegg.v7.0.entrez.gmt")

#parser$add_argument("--nfeatures",
#                    type="integer",
#                    default=7000,
#                    help="the dataset  to be used")

#parser$add_argument("--vst",
#                    action="store_true",
#                    default=FALSE)
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

print("### Loading data")
file=basename(args$data)
prefix=str_split(file,"\\.")[[1]][2]
#sep=ifelse(prefix%in%c("txt","tsv"),",","\t")
Data=read.csv(args$data,sep="\t",row.names=1)  #also can be use for seurat slot data
print(head(Data[,1:10]))
metadata=read.csv(args$metadata,sep=",")   # also can be from seurat metadata
stopifnot(args$column%in%colnames(metadata))

print("### Make group ")
grouP=metadata[[args$column]] %>% as.factor()
desigN <- model.matrix(~ grouP + 0)
head(desigN)
saveRDS(desigN,file.path(args$outdir,"desigN.rds"))
colnames(desigN)=c("grouPA","grouPB")
rownames(desigN) <- colnames(Data)
comparE <- makeContrasts(grouPB - grouPA, levels=desigN)

print("### Loading gmt")
geneSet <- getGmt(args$gmt)
print("### Run gsva")
gEs <- gsva(expr=as.matrix(Data), gset.idx.list=geneSet,kcdf="Gaussian", parallel.sz=12)
# pheatmap(gEs) # do pheatmap

print("### Run with gsva result for limma method")
fiT <- lmFit(gEs, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
Diff <- topTable(fiT3, coef=1, number=500)
saveRDS(fiT3,file.path(args$outdir,"gsvaEbayes_fit.rds"))
write.table(Diff,file.path(args$outdir,"gsva_difftable.csv"),sep=",",quote=F)

