library(argparse)
library(stringr)
library(Seurat)
library(Matrix)
library(mindr)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
#############################
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--cellchat",
                    type="character",
                    default="")

parser$add_argument("--np",
                    type="integer",
                    default=5)

parser$add_argument("--signaling",
                    type="character",
                    default="MIF")

parser$add_argument("--outdir",
                    type="character",
                    default="output")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

options(stringsAsFactors = FALSE)
cellchat=readRDS(args$cellchat)
groupSize <- as.numeric(table(cellchat@idents))
nPatterns=args$np
pathways.show=args$signaling

pdf(file.path(args$outdir,paste0(pathways.show,"_netVisual.pdf")),width=16,height=16)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
dev.off()

cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") 

pdf(file.path(args$outdir,"outgoing_Patterns.pdf"),width=8,height=16)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,width=8,height=16)
dev.off()

pdf(file.path(args$outdir,"incoming_Patterns.pdf"),width=8,height=16)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width=8,height=16)
dev.off()


pdf(file.path(args$outdir,"outgoing_netAnalysis_river.pdf"),width=16,height=16)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf(file.path(args$outdir,"outgoing_netAnalysis_dot.pdf"),width=16,height=16)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

pdf(file.path(args$outdir,"incoming_netAnalysis_river.pdf"),width=16,height=16)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

pdf(file.path(args$outdir,"incoming_netAnalysis_dot.pdf"),width=16,height=16)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()





