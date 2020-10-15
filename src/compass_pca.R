# Reproduces Fig. 2(e) from our paper. (This script works out-of-the-box with the Th17 cell data set included in the package.)

library(compassR)
library(ggrepel)
library(tidyverse)
library(ggplot2)
library(stringr)
library(argparse)
library(Cairo)
library(factoextra)
library(FactoMineR)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default=NULL,
                    help="directory should include files named :cell_metadata.csv, reactions.tsv, and linear_gene_expression_matrix.tsv")


parser$add_argument("--outdir",
                    type="character",
                    default="compass_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


compass_data=readRDS(args$data)
cell_meta=as.data.frame(compass_data$cell_metadata)
rownames(cell_meta)=cell_meta$cell_id
compass_score=compass_data$reaction_consistencies
compass_score=compass_score[,cell_meta$cell_id]

res.pca <- PCA(t(compass_score), graph = FALSE)
coord=res.pca$ind$coord

eig.val <- get_eigenvalue(res.pca)
coord.vars<-eig.val[colnames(coord),]

#Data=cbind(as.data.frame(coord),as.data.frame(coord.vars))
Data=as.data.frame(coord)
Data=Data[cell_meta$cell_id,]
Data$Status=cell_meta$status

ggplot(data=Data, aes(x=Dim.1,y=Dim.2)) +
  geom_point(aes(color = Status)) +
  geom_smooth( method = "lm",color = "black")+theme_bw()

ggsave(file.path(args$outdir,"compass_pca_1_2.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()

ggplot(data=Data, aes(x=Dim.1,y=Dim.3)) +
  geom_point(aes(color = Status)) +
  geom_smooth( method = "lm",color = "black")+theme_bw()

ggsave(file.path(args$outdir,"compass_pca_1_3.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()


ggplot(data=Data, aes(x=Dim.2,y=Dim.3)) +
  geom_point(aes(color = Status)) +
  geom_smooth( method = "lm",color = "black")+theme_bw()

ggsave(file.path(args$outdir,"compass_pca_2_3.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()

#fviz_pca_ind(res.pca,
#     geom.ind = "point", # show points only (nbut not "text")
#     col.ind = cell_meta$status, # color by groups
#     palette = c("#00AFBB", "#E7B800"),
#     addEllipses = TRUE, # Concentration ellipses
#     legend.title = "Status"
#     )
