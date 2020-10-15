library(SingleCellSignalR)
library(Seurat)

library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--column",
                    type="character",
                    default="label_main")

parser$add_argument("--gene",
                    type="character",
                    default="S1PR1")

parser$add_argument("--cluster",
                    type="character",
                    default="TC")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


seurat=readRDS(args$seurat)
seurat <- NormalizeData(seurat,scale.factor = 10000)

genes<-rownames(seurat)
all.genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]

# Retreiving the results of the preprocessing from the Seurat object
Idents(seurat)=seurat@meta.data[[args$column]]
cluster = as.character(Idents(seurat))
data = data.frame(seurat[["RNA"]]@data)

print("# Ligand/Receptor analysis using SingleCellSignalR")
signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)

print("# Visualization")
pdf(file.path(args$outdir,paste0("signal_",args$gene,"_",args$cluster,".pdf"),width=12,height=16)
visualize(signal)
dev.off()

pdf(file.path(args$outdir,paste0("intra_network_",args$gene,"_",args$cluster,".pdf"),width=12,height=16)
intra = intra_network(args$gene,data,all.genes,cluster,args$cluster,signal = signal)
dev.off()
