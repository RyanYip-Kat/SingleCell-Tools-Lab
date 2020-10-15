library(Matrix)
library(Seurat)
library(data.table)
library(magrittr)
library(argparse)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--mtx",
                    type="character",
                    default=NULL,
                    help="the path to save result")

parser$add_argument("--barcode",
                    type="character",
                    default=NULL,
                    help="the path of barcode file")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

################
read_mtx_scATACpro <- function(mtx_path){
  #mtx_path <- paste0(dirt, "matrix.mtx")
  mtx.dir = dirname(mtx_path)
  feature_path <- paste0(mtx.dir, "/features.txt")
  barcode_path <- paste0(mtx.dir, "/barcodes.txt")


  features <- fread(feature_path, header = F)
  barcodes <- fread(barcode_path, header = F)

  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$V1)%>%
    magrittr::set_colnames(barcodes$V1)

  return(mtx)
}
###############

write_mtx<-function(mtx,barcode=NULL,outdir="outdir"){
	#counts=GetAssayData(seurat,"counts")
	print("### Loadind mtx")
	counts=read_mtx_scATACpro(mtx)
	if(!is.null(barcode)){
		print(paste0("Subset : ",length(barcode)," barcodes"))
		counts=counts[,barcode]
	}
	features=rownames(counts)
	barcodes=colnames(counts)
	writeMM(obj = counts, file=file.path(outdir,"matrix.mtx"))
	# save genes and cells names
        write(x = features, file = file.path(outdir,"features.txt"))
        write(x = barcodes, file =file.path(outdir, "barcodes.txt"))
}

print("loading data")
DATA=read.csv(args$barcode,header=T,stringsAsFactors=F)
barcode=DATA$Barcode

print("write data")
write_mtx(args$mtx,barcode,args$outdir)



