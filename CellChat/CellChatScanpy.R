library(argparse)
library(stringr)
library(Seurat)
library(Matrix)
library(mindr)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(reticulate)
source("myfun.R")
#############################
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--adata",
                    type="character",
                    default="")


parser$add_argument("--outdir",
                    type="character",
                    default="results")

parser$add_argument("--group",
                    type="character",
                    default="label_main")

parser$add_argument("--search",
                    type="character",
                    default="Secreted Signaling",
		    choices=c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

options(stringsAsFactors = FALSE)
ad <- import("anndata", convert = FALSE)
ad_object <- ad$read_h5ad(args$adata)
# access normalized data matrix
data.input <- t(py_to_r(ad_object$X))
rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))
# access meta data
meta.data <- py_to_r(ad_object$obs)

print("# Create CellData")
cellchat <- createCellChat(data = data.input)
cellchat <- addMeta(cellchat, meta = meta.data)
cellchat <- setIdent(cellchat, ident.use =args$group) # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

print("# Loadin cellphone database")
CellChatDB <- CellChatDB.human

CellChatDB.use <- subsetDB(CellChatDB, search =args$search) # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

print("# Preprocess cellchat")
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

print("# Run Model")
#mycomputeCommunProb <-edit(computeCommunProb)  # computeCommunProb内部似乎有一些bug，同一套数据在window10上没事，到了Linux上有报错。发现是computeExpr_antagonist这个函数有问题，(matrix(1, nrow = 1, ncol = length((group))))，中应为(matrix(1, nrow = 1, ncol = length(unique(group))))？ 不然矩阵返回的不对。de了它。
#environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- computeCommunProb(cellchat)  # 这儿是我de过的。

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat,file.path(args$outdir,"cellchat.rds"))
