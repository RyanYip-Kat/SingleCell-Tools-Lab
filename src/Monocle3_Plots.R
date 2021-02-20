library(monocle3)
library(Seurat)
library(argparse)
library(stringr)
library(ggplot2)
library(Cairo)
library(plotly)
####################################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

parser <- ArgumentParser(description='Program to handle monocle3 pstime')
parser$add_argument("--outdir",
                    type="character",
                    default="./outdir")

args <- parser$parse_args()
outDir=args$outdir

Width=12
Height=10

makedir(outDir)
message("INFO : Loading dataset")
cds=readRDS("ips-result-subset1/Pseudotime-cds.rds")
cds_3d=readRDS("ips-result-subset1/Pseudotime-cds_3d.rds")

p=plot_cells(cds,
           color_cells_by = "Labels",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

ggsave(file.path(outDir,"Plot1.pdf"),plot=p,width=Width,height=Height,device =cairo_pdf)


p=plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(file.path(outDir,"Plot3.pdf"),plot=p,width=Width,height=Height,device =cairo_pdf)

p=plot_cells(cds,
           color_cells_by = "Labels",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

ggsave(file.path(outDir,"Plot2.pdf"),plot=p,width=Width,height=Height,device =cairo_pdf)


p=plot_cells_3d(cds_3d,
		reduction_method="UMAP",
		color_cells_by="Labels",
		color_palette="Paired")

htmlwidgets::saveWidget(p, "Plot1.html")

p=plot_cells_3d(cds_3d,
                reduction_method="UMAP",
                color_cells_by="pseudotime")

htmlwidgets::saveWidget(p, "Plot2.html")
