library(Seurat)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--barcode",
                    type="character",
                    default="")

parser$add_argument("--name",
                    type="character",
                    default="celltype")

parser$add_argument("--reduction",
                    type="character",
                    default=NULL,help="tsne,umap")

parser$add_argument("--path",
                    type="character",
                    default="")
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print("### Loading Dataset")
#path="/home/ye/Work/BioAligment/10X/output/vkh_513/outs"
path=args$path
counts=Read10X(file.path(path,"filtered_feature_bc_matrix"))
cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"))
rownames(cluster)<-cluster$Barcode


cluster<-subset(cluster,select=Cluster)
cell=rownames(cluster)

rownames(cluster)=colnames(counts)
print(paste0("Size of counts  [ ",nrow(counts),",",ncol(counts)," ]"))
print("### Create Seurat object")
genes<-rownames(counts)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
keep_genes=keep_genes[!str_detect(keep_genes,"\\.")]
counts=counts[keep_genes,]
print(dim(counts))
object<-CreateSeuratObject(counts= counts,
                       assay = "RNA",
                       project ="scRNA",
                       names.delim="_",
                       min.cells=0,
                       min.features=0)

object<-AddMetaData(object,metadata=cluster,col.name="orig.Cluster")
cells<-colnames(object)

ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
object$ident<-ident
#status<-ifelse(ident%in%c(1:8),"AA","YA")
#print(table(status))
#object$status<-status


print("### Subset")
DATA=read.csv(args$barcode,stringsAsFactors=F)
colnames(DATA)=c("barcode","celltype")
rownames(DATA)<-DATA$barcode


DATA<-subset(DATA,select=celltype)
cell=rownames(DATA)
object<-subset(object,cells=cell)
object<-AddMetaData(object,metadata=DATA,col.name=args$name)

######################

object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)

print(paste0("Size of object after selection [ ",nrow(object),",",ncol(object)," ]"))
print("### Normalize Data")
object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = FALSE)
#object=SCTransform(object,vars.to.regress="nFeature_RNA",verbose = FALSE)

#print("### Scale Data")
#object<-ScaleData(object,features=VariableFeatures(object),model.use = "linear",
#               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)
#########################
print("#### Add reduction")
mat=read.csv(args$reduction,stringsAsFactors=F)
rownames(mat)=mat$Barcode
mat=mat[colnames(object),]
mat=mat[,-1]
pcs=colnames(mat)
mat=as.matrix(mat)
#colnames(mat)<-paste(paste0(args$reduction,"_"),1:ncol(mat),sep = "")
key=str_split(colnames(mat)[1],"_")[[1]][1]
r=str_to_lower(key)
object[[r]]<-CreateDimReducObject(embeddings =mat,
                                  key =paste0(key,"_"),
                                  assay = DefaultAssay(object))


m<-Embeddings(object,r)
print(dim(m))
saveRDS(object,file.path(args$outdir,"seurat.rds"))
#cds<-as.CellDataSet(object)
#saveRDS(cds,file.path(args$outdir,"cds.rds")

