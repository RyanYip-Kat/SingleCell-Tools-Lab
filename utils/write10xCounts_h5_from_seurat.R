library(argparse)
library(stringr)
library(Seurat)
library(DropletUtils)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default=NULL)

parser$add_argument("--rename",
                    action="store_true",
                    default=FALSE)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat=readRDS(args$seurat)

counts=GetAssayData(seurat,"counts")
rownames(counts)=str_to_title(rownames(counts))
if(args$rename){
	cellNames=unlist(lapply(colnames(counts),function(cell)return(str_split(cell,"-")[[1]][1])))
	cellNames=paste(cellNames,"1",sep="-")
	colnames(counts)=cellNames
}
print(head(rownames(counts)))
print(head(colnames(counts)))
filename=file.path(args$outdir,"filtered_feature_bc_matrix.h5")
print(paste0("### Write counts matrix into : ",filename))
write10xCounts(x =counts, path=filename,type="HDF5",overwrite=TRUE,version="3")
filename=file.path(args$outdir,"matrix")
print(paste0("### Write counts matrix into : ",filename))
write10xCounts(x =counts, path=filename,overwrite=TRUE,version="3")

