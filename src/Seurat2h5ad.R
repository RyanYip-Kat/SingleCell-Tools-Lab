library(SeuratDisk)
library(SeuratWrappers)
library(argparse)
library(stringr)
####################################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

parser <- ArgumentParser(description='Program to handle monocle3 pstime')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="monocle3 object")

parser$add_argument("--outdir",
                    type="character",
                    default="./outdir")

parser$add_argument("--prefix",
                    type="character",
                    default="seurat",
		    help="prefix name for convert object")

args <- parser$parse_args()

outDir=args$outdir
prefix=args$prefix
makedir(outDir)

message("INFO : Loading dataset ...")
seurat=readRDS(args$seurat)

message("INFO : SaveH5Seurat ...")
dest=file.path(outDir,paste0(prefix,".h5Seurat"))
SaveH5Seurat(seurat, filename = dest)

message("INFO : Convert into h5ad ...")
Convert(dest, dest = "h5ad")

message("INFO : Done!")

