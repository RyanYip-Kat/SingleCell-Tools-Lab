library(Seurat)
library(pheatmap)
#library(NMF)
library(dplyr)
library(viridis)
#library(hrbrthemes)
#library(igraph)
library(RColorBrewer)
library(stringr)
library(argparse)

set.seed(7777)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--slot",
                    type="character",
                    default="data",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--status",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="")
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}
combn.ident<-function(idents){
  df<-expand.grid(idents,idents)
  idx<-apply(df,1,function(z){
    if(z[1]==z[2]){
      return(FALSE)
    }else{
      return(TRUE)
    }
  })
  df<-df[idx,]
  rownames(df)<-1:nrow(df)
  return(t(df))
}
Seuratinfo<-function(seurat){
  value<-seurat
  n<-ncol(seurat) # number of cells
  c<-colnames(seurat) # Cells of seurat
  g<-rownames(seurat) # Genes of seurat
  return(list("v"=value,"n"=n,"c"=c,"g"=g))
}

InteractionScore<-function(sL,sR,Lg,Rg){
  # L : Ligand
  # R : Receptor
  vL<-sL[["v"]]
  nL<-sL[["n"]]
  Lc<-sL[["c"]]
  
  vR<-sR[["v"]]
  nR<-sR[["n"]]
  Rc<-sR[["c"]]
  
  Le<-vL[Lg,Lc]
  Re<-vR[Rg,Rc]
  
  p1<-sum(Le)
  p2<-sum(Re)
  
  I<-p1*p2/(nL*nR)
  return(I)
}

Score<-function(sL,sR,interaction_csv){
  score<-c()
  row<-c()
  pairs<-read.csv(interaction_csv,header=TRUE)
  colnames(pairs)<-c("Ligand","Receptor")
  n<-nrow(pairs)
  ligand<-as.character(pairs$Ligand)
  receptor<-as.character(pairs$Receptor)
  eps<-1e-6
  
  Lgs<-sL[["g"]]
  Rgs<-sR[["g"]]
  for(i in 1:n){
    Lg=ligand[i]
    Rg=receptor[i]
    if(Lg%in%Lgs & Rg%in%Rgs){
      I<-InteractionScore(sL,sR,Lg,Rg)
      score<-c(score,I)
      name<-paste0(Lg,"-",Rg)
      row<-c(row,name)
    }
  }
  names(score)<-row
  #score<-as.data.frame(score)
  #rownames(score)<-row
  #colnames(score)<-col.name
  return(score)
}

GetSeuratInteractionScore<-function(seurat,idents,interaction_csv,slot){
  if(length(idents)<2){
    stop("Input idents must more than 2")
  }
  eps<-1e-4
  scores<-list()
  
  mat<-GetAssayData(seurat,slot)
  idents.combn<-combn.ident(idents)
  
  n.combn<-ncol(idents.combn)
  for(i in 1:n.combn){
    sL<-mat[,WhichCells(seurat,idents=idents.combn[1,i])]
    sR<-mat[,WhichCells(seurat,idents=idents.combn[2,i])]
    sLinfo<-Seuratinfo(sL)
    sRinfo<-Seuratinfo(sR)
    
    I1<-Score(sLinfo,sRinfo,interaction_csv)
    I2<-Score(sRinfo,sLinfo,interaction_csv)
    
    s<-data.frame(log2((I1+eps)/(I2+eps)))
    colnames(s)<-paste0(idents.combn[1,i],"_",idents.combn[2,i])
    scores[[i]]<-s
  }
  return(scores)
}


interaction_csv<-"interaction.csv"
bc<-readRDS("../rscript/output/aging-11/BC/seurat.rds")
dc<-readRDS("../rscript/output/aging-11/DC/seurat.rds")
nk<-readRDS("../rscript/output/aging-11/NK/seurat.rds")
cd8<-readRDS("../rscript/output/aging-11/CD8/seurat.rds")
tc<-readRDS("../rscript/output/aging-11/TC/seurat.rds")


bc$celltype<-"BC"  # BC
dc$celltype<-"DC" # TC
nk$celltype<-"NK"
cd8$celltype="CD8"
tc$celltype="TC"

seurat<-merge(x=bc,y=c(dc,nk,cd8,tc),merge.data=TRUE)
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)
print(table(seurat$celltype))
if(!is.null(args$status)){
	Idents(seurat)<-seurat$status
	seurat<-subset(seurat,idents=args$status)
}
Idents(seurat)<-seurat$celltype

idents<-c("CD8","BC","DC","NK","TC")
scores<-GetSeuratInteractionScore(seurat = seurat,
                                  idents=idents ,
                                  interaction_csv = interaction_csv,
                                  slot="data")

scores<-do.call(cbind,scores)
saveRDS(scores,file.path(args$outdir,"InteractionScores.rds"))
#idx<-apply(scores,1,function(score){
#  if(all(score==0)){
#    return(TRUE)
#  }else{
#    return(FALSE)
#  }
#})
#scores<-scores[!idx,]
mat=as.matrix(scores)
cell2cell<-c("BC_DC","CD8_DC","TC_DC","NK_DC")
mat<-mat[,cell2cell]
write.table(mat,file.path(args$outdir,"scores.csv"),sep=",",quote=F)
pdf(file.path(args$outdir,"InteractionScore-1.pdf"),width=8,height=16)
pheatmap(mat[1:84,],
         border_color=NA,
         cluster_rows =F,
         cluster_cols = F,
         fontsize_row = 10,
         angle_col=90,
         fontsize =10)
dev.off()

pdf(file.path(args$outdir,"InteractionScore-2.pdf"),width=8,height=16)
pheatmap(mat[85:nrow(mat),],
           border_color=NA,
           cluster_rows =F,
           cluster_cols = F,
           fontsize_row = 10,
           angle_col=90,
           fontsize =10)
dev.off()

pdf(file.path(args$outdir,"InteractionScore-3.pdf"),width=8,height=32)
pheatmap(mat,
             border_color=NA,
             cluster_rows =F,
             cluster_cols = F,
             fontsize_row = 10,
             angle_col=90,
             fontsize_col =12)
dev.off()
