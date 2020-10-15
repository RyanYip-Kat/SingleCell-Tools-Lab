library(argparse)
library(stringr)
library(Seurat)
library(Matrix)
library(mindr)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
#source("myfun.R")
#############################
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--column1",
                    type="character",
                    default=NULL)

parser$add_argument("--subset1",
                    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--column2",
                    type="character",
                    default=NULL)

parser$add_argument("--subset2",
                    nargs="+",
                    type="character",
                    default=NULL)

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



parser$add_argument("--nfeatures",
                    type="integer",
                    default=NULL,
                    help="the dataset  to be used")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

options(stringsAsFactors = FALSE)
seurat_obj<-readRDS(args$seurat)
genes=rownames(seurat_obj)
keep_genes<-genes[!str_detect(genes,"^MT-|^RPL|^RPS")]
keep_genes=keep_genes[!str_detect(keep_genes,"\\.")]

seurat_obj=subset(seurat_obj,features=keep_genes)
print(dim(seurat_obj))
DefaultAssay(seurat_obj)="RNA"


metadata=seurat_obj@meta.data
if(!is.null(args$column1)){
        if(!args$column1%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset1,collapse=",")
                print(paste0("Subset :",args$column1," with :",s))
                Idents(seurat_obj)=metadata[[args$column1]]
                seurat_obj<-subset(seurat_obj,idents=args$subset1)
        }
	metadata=seurat_obj@meta.data
        if(!is.null(args$column2)){
                if(!args$column2%in%colnames(metadata)){
                stop("Invaild columns in metadata!")
        }else{
                s<-paste(args$subset2,collapse=",")
                print(paste0("Subset :",args$column2," with :",s))
                Idents(seurat_obj)=metadata[[args$column2]]
                seurat_obj<-subset(seurat_obj,idents=args$subset2)
                }
        }
}

if(!is.null(args$nfeatures)){
	seurat_obj<-FindVariableFeatures(seurat_obj,nfeatures=args$nfeatures)
	keep_genes=VariableFeatures(seurat_obj)
}else{
	keep_genes=rownames(seurat_obj)
}

print("# Analysis ")
#args<-list()
metadata=seurat_obj@meta.data
seurat_obj<-NormalizeData(seurat_obj,normalization.method = "LogNormalize",verbose = FALSE)
data.input=GetAssayData(seurat_obj,"data")[keep_genes,]

print("# Create CellData")
cellchat <- createCellChat(data = data.input)
cellchat <- addMeta(cellchat, meta = metadata)
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
