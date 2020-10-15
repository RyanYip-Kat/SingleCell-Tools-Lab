library(compassR)
library(ggrepel)
library(stringr)
library(argparse)
library(Seurat)
library(ggplot2)
library(pheatmap)
source("compass_function.R")

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default=NULL,
                    help="")


parser$add_argument("--outdir",
                    type="character",
                    default="compass_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

compass_data=readRDS(args$data)
mat=compass_data$reaction_consistencies
seurat=CreateSeuratObject(mat)

idents=unlist(lapply(Cells(seurat),function(cell){return(str_split(cell,"-")[[1]][2])}))
print(table(idents))
seurat$idents=idents
Idents(seurat)=seurat$idents

cts=AverageExpression(seurat,features=rownames(seurat))[[DefaultAssay(seurat)]]
pdf(file.path(args$outdir,"block_heatmap.pdf"),width=8,height=64)
p=pheatmap(cts,,
            border_color=NA,
            show_rownames=T,
            show_colnames=T,
	    scale="column",
            cluster_rows=T,
            cluster_cols=F,fontsize=7)

print(p)
dev.off()

