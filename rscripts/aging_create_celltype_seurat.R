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
                    default=NULL)


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

cluster=read.csv(args$barcode,stringsAsFactors=F)
seurat=readRDS("/home/ye/Work/R/SingleCell/Su/human/VDJ/Project/VKH/20200509/output/aging-11/model/seurat.rds")

rownames(cluster)<-cluster$barcode


cluster<-subset(cluster,select=celltype)
cell=rownames(cluster)

genes<-rownames(seurat)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

seurat<-subset(seurat,cells=cell,features=keep_genes)

seurat<-AddMetaData(seurat,metadata=cluster,col.name="celltype")

seurat <- FindVariableFeatures(seurat, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)
print("### Normalize Data")
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)

print("### Scale Data")
seurat<-ScaleData(seurat,features=VariableFeatures(seurat),model.use = "linear",
               vars.to.regress = c("nFeature_RNA"),verbose =FALSE)


print("### saving")
saveRDS(seurat,file.path(args$outdir,"seurat.rds"))
