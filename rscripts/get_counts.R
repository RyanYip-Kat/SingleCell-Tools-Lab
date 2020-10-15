library(argparse)
library(stringr)
library(Seurat)
library(DropletUtils)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--meta",
                    type="character",
                    default="")

parser$add_argument("--method",
                    type="character",
                    default="union",
		    choices=c("union","intersect"))

parser$add_argument("--min_cells",
                    type="integer",
                    default=0)

parser$add_argument("--nfeatures",
                    type="integer",
                    default=NULL)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

meta=read.csv(args$meta,header=TRUE,stringsAsFactors=FALSE) # Path,Status,Sample
Path=meta$Path
Status=meta$Status
Sample=meta$Sample
time_point=meta$Time
N=nrow(meta)

count_list<-list()
meta_list<-list()
gene_list<-list()
for(i in 1:N){
	count=Read10X(Path[i])
	cell_name=paste(Sample[i],colnames(count),sep="#")
	colnames(count)=cell_name
	print(dim(count))
	status=rep(Status[i],ncol(count))
	sam=rep(Sample[i],ncol(count))
	time=rep(time_point[i],ncol(count))
	metadata<-data.frame("cell_id"=colnames(count),"Sample"=sam,"Status"=status,"Time"=time)
	rownames(metadata)=colnames(count)
	count_list[[i]]=count
	meta_list[[i]]=metadata
	gene_list[[i]]=rownames(count)
}

Metadata=do.call(rbind,meta_list)
#Features=unique(unlist(gene_list))
if(args$method=="union"){
	Features=Reduce(union,gene_list)
	counts=lapply(count_list,function(read){
                     feature=rownames(read)
                     impute_feature=setdiff(Features,feature)
                     impute_matrix=matrix(0,nrow=length(impute_feature),ncol=ncol(read))
                     rownames(impute_matrix)=impute_feature
                     colnames(impute_matrix)=colnames(read)
                     impute_matrix=as(impute_matrix, "dgCMatrix")
                     m=rbind(read,impute_matrix)
                     return(m)})
}else{
	Features=Reduce(intersect,gene_list)
	counts=lapply(count_list,function(read)return(read[Features,]))
}

mat<-do.call(cbind,counts)
Metadata=Metadata[colnames(mat),]
write.table(Metadata,file.path(args$outdir,"metadata.csv"),sep=",",quote=F,row.names=F)
seurat=CreateSeuratObject(mat,min.cells=args$min_cells)
if(!is.null(args$nfeatures)){
	seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
	seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures =args$nfeatures)
	seurat<-subset(seurat,features=VariableFeatures(seurat))
}

write10xCounts(x =GetAssayData(seurat,"counts"), path=file.path(args$outdir,"matrix"),overwrite=TRUE)

