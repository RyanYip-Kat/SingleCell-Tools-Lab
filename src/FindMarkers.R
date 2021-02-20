library(Seurat)
library(argparse)
library(stringr)
library(dplyr)
library(ggplot2)
library(future)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")


parser$add_argument("--outdir",
                    type="character",
                    default="./Results")


parser$add_argument("--assay",
		    type="character",
		    default="RNA",
		    help="which assay to be used")

parser$add_argument("--groupby",
		    type="character",
		    default="seurat_clusters",
		    help="which column used to be caculating DE markers")

parser$add_argument("--n_genes",
		    type="integer",
		    default=10,
		    help="the number of top genes")

parser$add_argument("--heatmap",
		    action="store_true",
		    default=FALSE,
		    help="whether to plot heatmap")
args <- parser$parse_args()

############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
makedir(outDir)

###############################
message("INFO : loading dataset ...")
seurat=readRDS(args$seurat)
DefaultAssay(seurat)=args$assay
metadata=seurat@meta.data

assertthat::assert_that(args$groupby%in%colnames(metadata))
Idents(seurat)=metadata[[args$groupby]]

###############################
plan("multiprocess", workers = 16)
mtFeatures=rownames(seurat)[str_detect(rownames(seurat),"^mt-|^MT-")]
markers=FindAllMarkers(object=seurat,
		       assay=args$assay,
		       features=mtFeatures,
		       logfc.threshold = 0.25,
		       test.use = "wilcox",
		       slot="data")

write.table(markers,file.path(outDir,"All.markers.csv"),sep=",",quote=FALSE,row.names=F)
if(args$heatmap){
	message("INFO : Do Heatmap ...")
	topN <- markers %>% group_by(cluster) %>% top_n(n = as.integer(args$n_genes), wt = avg_logFC)
	write.table(topN,file.path(outDir,paste0("top",args$n_genes,"_markers.csv")),sep=",",quote=FALSE,row.names=F)
        p=DoHeatmap(seurat, features = topN$gene) #+ NoLegend()
	ggsave(file.path(outDir,"heatmap.pdf"),plot=p,width=12,height=16)
}
message("INFO : Done!")
