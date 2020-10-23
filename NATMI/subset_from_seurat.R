library(Seurat)
library(argparse)
library(stringr)
library(DropletUtils)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the path of seurat has been created")

parser$add_argument("--column1",
                    type="character",
                    default=NULL,
                    help="which column to use subset")

parser$add_argument("--subset1",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset")

parser$add_argument("--column2",
                    type="character",
                    default=NULL,
                    help="which column to use subset")

parser$add_argument("--subset2",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset")

parser$add_argument("--nfeatures",
                    type="integer",
                    default=NULL)

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
seurat_obj<-readRDS(args$seurat)
if(!is.null(args$column1) & !is.null(args$subset1)){
	metadata=seurat_obj@meta.data
        if(!args$column1%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset1,collapse=",")
                print(paste0("Subset :",args$column1," with :",s))
                Idents(seurat_obj)=metadata[[args$column1]]
                seurat_obj<-subset(seurat_obj,idents=args$subset1)

		if(!is.null(args$column2) & !is.null(args$subset2)){
			metadata=seurat_obj@meta.data
			if(!args$column2%in%colnames(metadata)){
                                stop("Invaild columns in metadata!")
			}else{
				s<-paste(args$subset2,collapse=",")
				print(paste0("Subset :",args$column2," with :",s))
				Idents(seurat_obj)=metadata[[args$column2]]
				seurat_obj<-subset(seurat_obj,idents=args$subset2)
			}
		}
	}
}

if(!is.null(args$nfeatures)){
        seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures =args$nfeatures)
        seurat_obj <-subset(seurat_obj,features=VariableFeatures(seurat_obj))
}

#seurat_obj=SCTransform(seurat_obj)
metadata=seurat_obj@meta.data
meta=metadata[,c("item","label")]
print(table(meta$label))
write.table(meta,file.path(args$outdir,"label.csv"),sep=",",quote=FALSE,row.names=FALSE)
write10xCounts(x =GetAssayData(seurat_obj,assay="RNA",slot="counts"), path=file.path(args$outdir,"matrix"),overwrite=TRUE)





