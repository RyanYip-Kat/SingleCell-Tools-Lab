library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--matrix",
                    type="character",
                    default="")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

df=read.csv(args$matrix,stringsAsFactors=F)
mat=as.matrix(table(df[,1],df[,2]))
write.table(mat,file.path(args$outdir,"ident_cluster_number.csv"),sep=",",quote=F)

