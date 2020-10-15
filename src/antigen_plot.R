library(ggplot2)
library(argparse)
library(Seurat)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--pt_size",
                    type="double",
                    default="1.0")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=25),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

# cols = c("#191970","yellow","#FF3030")
seurat<-readRDS(args$seurat)
antigens<-rownames(seurat)
for(antigen in antigens){
        jpeg(file.path(args$outdir,paste0(antigen,"_tsne_condition_seurat.jpeg")),width=1024,height=1024)
        p=FeaturePlot(seurat,reduction="tsne",pt.size=args$pt_size,features=antigen,split.by="condition")+
                my_theme
        print(p)
        dev.off()

        jpeg(file.path(args$outdir,paste0(antigen,"_tsne_seurat.jpeg")),width=1024,height=1024)
        p=FeaturePlot(seurat,reduction="tsne",pt.size=args$pt_size,features=antigen)+my_theme
        print(p)
        dev.off()


}

jpeg(file.path(args$outdir,"tsne_cluster.jpeg"),width=1024,height=1024)
DimPlot(seurat,reduction="tsne",pt.size=args$pt_size,label=TRUE)+my_theme
dev.off()

