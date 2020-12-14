library(stringr)
library(ggplot2)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(argparse)
library(future)
#############################
parser <- ArgumentParser(description='A Program to Perform NicheNet analysis starting from a Seurat object ')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat object")


parser$add_argument("--outdir",
                    type="character",
                    default="Result")


parser$add_argument("--condition",
		    type="character",
		    default=NULL,
		    help="condition colname in Seurat metadata to  explain differential expression between two conditions,length must be two!!!")

parser$add_argument("--condition_io",
                    type="character",
                    default=NULL,
		    help="condition of interest name")

parser$add_argument("--condition_reference",
                    type="character",
                    default=NULL,
                    help="condition of reference name")


parser$add_argument("--receiver",
                    type="character",
                    default=NULL,
		    help="receiver celltype name,like : CD8 T,CD4 T,...")

parser$add_argument("--sender",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="sender celltype name,like : CD8 T,CD4 T,...")


args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

##############################  Paramenters
outDir=args$outdir
condition=args$condition
condition_io=args$condition_io
condition_reference=args$condition_reference
receiver=args$receiver
sender=args$sender

makedir(outDir)
############################## Configure files
message("INFO : Loading Nichenetr Configure ...")
ligand_target_matrix=readRDS("/home/ye/Work/R/SingleCell/Nichenetr/ligand_target_matrix.rds")
weighted_networks=readRDS("/home/ye/Work/R/SingleCell/Nichenetr/weighted_networks.rds")
lr_network=readRDS("/home/ye/Work/R/SingleCell/Nichenetr/lr_network.rds")

##############################
message("INFO : Loading dataset ...")
seuratObj=readRDS(args$seurat)
seuratObj$celltype=seuratObj$label_fine
Idents(seuratObj)=seuratObj$label_fine
message("INFO : Perform the NicheNet analysis ...")
cat(sprintf("INFO : Start [ %s ] ... \n",Sys.time()))
plan("multiprocess", workers = 16)
nichenet_output = nichenet_seuratobj_aggregate(seuratObj,
					       receiver=receiver,
					       condition_colname=condition,
					       condition_oi =condition_io, 
					       condition_reference =condition_reference,
					       sender=sender,
					       ligand_target_matrix = ligand_target_matrix, 
					       lr_network = lr_network, weighted_networks = weighted_networks, 
					       organism = "human")

cat(sprintf("INFO : End [ %s ] ... \n",Sys.time()))
message("INFO : Save Result ...")
saveRDS(nichenet_output,file.path(outDir,"nichenet_output.rds"))

message("INFO : Plot Results ...")
p=DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
ggplot(file.path(outDir,"top_ligands.pdf"),plot=p,height=12,width=18)

p=DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu",split.by=condition) + RotatedAxis()
ggplot(file.path(outDir,"top_ligands_condition.pdf"),plot=p,height=12,width=18)


p=nichenet_output$ligand_target_heatmap + 
	scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
	xlab("anti-LCMV response genes in CD8 T cells") + ylab("Prioritized immmune cell ligands")

ggplot(file.path(outDir,"ligand_target_heatmap.pdf"),plot=p,height=12,width=18)



