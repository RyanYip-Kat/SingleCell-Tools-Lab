library(Seurat)
library(argparse)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the path to save result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat=readRDS(args$seurat)
metadata=seurat@meta.data

#mat=table(seurat$sample_id,seurat$celltype)
mat=table(seurat$sample_id,seurat$seurat_clusters)
write.table(as.matrix(mat),file.path(args$outdir,"ident_cluster_numbers.csv"),sep=",",quote=F)

mat=table(seurat$seurat_clusters)
write.table(as.matrix(mat),file.path(args$outdir,"cluster_numbers.csv"),sep=",",quote=F)

