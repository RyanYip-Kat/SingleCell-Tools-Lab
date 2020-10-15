library(Seurat)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--path",
                    type="character",
                    default="")
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


path=args$path
print(paste0("### Loading Dataset from : ",path))
counts=Read10X(file.path(path,"filtered_feature_bc_matrix"))

filename=file.path(args$outdir,"expression.tsv")
print(paste0("### Saving expression matrix into : ",filename))
write.table(counts,filename,sep="\t",quote=FALSE) # df=read.table(filename,header=T,row.names=1)
