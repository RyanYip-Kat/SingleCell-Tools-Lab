library(argparse)
library(pheatmap)
library(ggplot2)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(RColorBrewer)
library(viridis)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--loom",
                    type="character",
                    default=NULL,
                    help="the path of pysceni result")


parser$add_argument("--outdir",
                    type="character",
                    default=NULL)
args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

loom <- open_loom(args$loom, mode="r")

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulonsAuc(loom)
regulonsAucThresholds <- get_regulonThresholds(loom)
embeddings <- get_embeddings(loom)

#exprMat <- get_dgem(loom)
cellInfo <- get_cellAnnotation(loom)
clusterings <- get_clusterings_withName(loom)

close_loom(loom)

auc_mtx=regulonsAUC@assays@data@listData$AUC
pheatmap.idents<-function(auc_mtx,features,clusters){
        cts <- auc_mtx
        features<-features[features%in%rownames(cts)]
        new_cluster <- sort(clusters)
        cts <- as.matrix(cts[features, names(new_cluster)])
        ac=data.frame(cluster=new_cluster)
        #color<-c("#40E0D0","#B0171F")
        p<-pheatmap(cts,color = viridis(100),show_colnames =F,show_rownames = F,
                    cluster_rows = F,
                    cluster_cols = F,
                    fontsize_row=10,
                    fontsize=8,
                    annotation_col=ac)
        return(p)
}

clusters=as.character(cellInfo$CellType)
names(clusters)=rownames(cellInfo)

pdf(file.path(args$outdir,"rss_auc_heatmap.pdf"),width=10,height=16)
pheatmap.idents(auc_mtx,rownames(auc_mtx),clusters)
dev.off()

jpeg(file.path(args$outdir,"rss_auc_heatmap.jpeg"),width=1024,height=1024)
pheatmap.idents(auc_mtx,rownames(auc_mtx),clusters)
dev.off()
