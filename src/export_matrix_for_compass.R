library(argparse)
library(stringr)
library(Seurat)


set.seed(7777)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--slot",
                    type="character",
                    default="data",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the dataset")

parser$add_argument("--column1",
                    type="character",
                    default=NULL)

parser$add_argument("--subset1",
                    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--column2",
                    type="character",
                    default=NULL)

parser$add_argument("--subset2",
                    nargs="+",
                    type="character",
                    default=NULL)


parser$add_argument("--nfeatures",
                    type="integer",
                    default=10000,
                    help="the dataset  to be used")

parser$add_argument("--num_sample",
                    type="integer",
                    default=100,
                    help="the number sample for each idents")


parser$add_argument("--vst",
                    action="store_true",
                    default=FALSE)

parser$add_argument("--genes",
                    type="character",
                    default=NULL,
		    help="csv file of genes file")

args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

print("# Loading dataset") 
seurat_obj=readRDS(args$seurat)
genes=rownames(seurat_obj)
point_genes=genes[!str_detect(genes,"\\.")]

keep_genes=point_genes[!str_detect(point_genes,"^RPL|^RPS|^MT-")]
seurat_obj=subset(seurat_obj,features=keep_genes)
print(dim(seurat_obj))

metadata=seurat_obj@meta.data
if(!is.null(args$column1)){
        if(!args$column1%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset1,collapse=",")
                print(paste0("Subset :",args$column1," with :",s))
                Idents(seurat_obj)=metadata[[args$column1]]
                seurat_obj<-subset(seurat_obj,idents=args$subset1)
        }
        metadata=seurat_obj@meta.data
        if(!is.null(args$column2)){
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

if(!is.null(args$num_sample)){
	cells=Cells(seurat_obj)
	samples=unique(seurat_obj$sample)
	cells_list=unlist(lapply(samples,function(s){
				  cell=cells[seurat_obj$sample==s]
				  sample_cell=sample(cell,size=args$num_sample)
				  return(sample_cell)}))
	seurat_obj=subset(seurat_obj,cells=cells_list)
}

if(!is.null(args$genes)){
	features=read.csv(args$genes,stringsAsFactors=F,header=F)$V1
	vst=FALSE
}else{
	vst=args$vst
}

if(vst){
	seurat_obj<-FindVariableFeatures(seurat_obj,nfeatures=args$nfeatures)
	features=VariableFeatures(seurat_obj)
}

seurat_obj<-subset(seurat_obj,features=features)
print(dim(seurat_obj))
print("### Get metrix data")
seurat_obj<-NormalizeData(seurat_obj,normalization.method = "LogNormalize",verbose = FALSE)
DATA<-GetAssayData(seurat_obj,args$slot)
cells=colnames(DATA)
cells=unlist(lapply(1:length(cells),function(i){return(paste0("SRR",i))}))
colnames(DATA)=cells
print(dim(DATA))

print("### Get Meta")
seurat_obj$cell_id=cells
metadata=seurat_obj@meta.data
meta.data=metadata[,c("cell_id","status","condition","label_fine")]

print("## Save cell meta")
meta_filename=file.path(args$outdir,"cell_metadata.csv")
write.table(meta.data,meta_filename,sep=",",row.names = FALSE,quote=F)

print("## Save linear_gene_expression_matrix")
matrix_filename=file.path(args$outdir,"linear_gene_expression_matrix.tsv")
write.table(DATA,matrix_filename,sep="\t",row.names =TRUE,quote=F)
