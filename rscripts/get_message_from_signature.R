library(argparse)
library(stringr)
library(ggplot2)
library(VISION)
library(Seurat)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default="vision_seurat.rds")

#parser$add_argument("--reduction",
#                    type="character",
#                    default="tsne",choices=c("tsne","umap"))

parser$add_argument("--sig",
		    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

vis=readRDS(args$data)
sigscores=getSignatureScores(vis)
object<-CreateSeuratObject(t(sigscores),assay="sigature")
print(head(rownames(object)))
embedding=getProjections(vis)

#pro=paste0("Seurat_",args$reduction)
#projection=embedding[[pro]]

mat=embedding[["Seurat_tsne"]]
colnames(mat)<-paste(paste0("tSNE","_"),1:ncol(mat),sep = "")
object[["tsne"]]<-CreateDimReducObject(embeddings =mat,
                                  key =paste0("tSNE","_"),
                                  assay = DefaultAssay(object))

mat=embedding[["Seurat_umap"]]
colnames(mat)<-paste(paste0("UMAP","_"),1:ncol(mat),sep = "")
object[["umap"]]<-CreateDimReducObject(embeddings =mat,
                                  key =paste0("UMAP","_"),
                                  assay = DefaultAssay(object))


saveRDS(object,file.path(args$outdir,"seurat_signature.rds"))
if(!is.null(args$sig)){
	for(sig in args$sig){
		filename=paste0(sig,"_tSNE.pdf")
		pdf(file.path(args$outdir,filename),width=16,height=12)
		print(FeaturePlot(object,reduction="tsne",features=sig,pt.size=2.0)+theme_bw()+NoGrid())
		dev.off()

		filename=paste0(sig,"_UMAP.pdf")
                pdf(file.path(args$outdir,filename),width=16,height=12)
                print(FeaturePlot(object,reduction="umap",features=sig,pt.size=2.0)+theme_bw()+NoGrid())
                dev.off()
         }
}





