library(argparse)
library(Seurat)
library(ggplot2)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--genes",
		    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--ident",
                    type="character",
                    default="celltype")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
if(is.null(args$genes)){
	gene_set<-c("DDIT4",
            "DUSP2",
            "ISG20",
            "S100A10",
            "COX5A",
            "PDCD5",
            "CASP4",
            "PSMB9",
            "CRIP1",
            "TMEM219",
            "PSMB6",
            "COX8A",
            "PSMA7",
            "ZFP36",
            "COX5B")
}else{
	gene_set=args$genes
}
print(paste0("Caculate score used : ",length(gene_set)," genes"))

ImmunoScore<-function(seurat,gene_set,sample_id="1"){
  metadata=seurat@meta.data
  if(!sampel_id%in%unique(metadata$idents)){
    stop("Invalid sample id!!!")
  }
  Idents(seurat)<-seurat$ident
  cells<-WhichCells(seurat,idents=sample_id)
  DATA<-GetAssayData(seurat,"data")[gene_set,]
  #DATA<-as.matrix(DATA)
  gs_data<-as.matrix(DATA[,cells])
  
  sd_j<-apply(DATA,1,sd)
  m_j<-apply(DATA,1,mean)
  
  Nes_ij<-(gs_data-m_j)/sd_j
  
}

GSscore<-function(seurat,gene_set,idents="celltype",cell_id="pDC"){
  metadata=seurat@meta.data
  Idents(seurat)<-metadata[[idents]]
  cells<-WhichCells(seurat,idents=cell_id)
  gs_seurat<-subset(seurat,cells = cells,features = gene_set)
  #gsij_umi<-sum(gs_seurat@meta.data$nCount_RNA)
  gsij_umi<-gs_seurat@meta.data$nCount_RNA
  #######################
  cj_seurat<-subset(seurat,idents = cell_id)
  cij_umi<-cj_seurat@meta.data$nCount_RNA
  #cij_umi<-sum(gs_seurat@meta.data$nCount_RNA)
  n_cells<-length(cells)
  Cj<-gsij_umi/cij_umi*100 #score of Cj cell for GSx gene set
  DATA<-data.frame("GS"=Cj,"Cluster"=cell_id)
  return(DATA)
}

metadata=seurat@meta.data
DATA<-lapply(unique(metadata[[args$ident]]),GSscore,seurat=seurat,idents=args$ident,gene_set=gene_set)
df<-do.call(rbind,DATA)
#df$Cluster=factor(df$Cluster,levels=c("CD4","CD8","NK","BC","Mono","DC")) # idents
#df$Cluster=factor(df$Cluster,levels=c(9:16,1:8))
#level=c("CD4 Naive","CD4 Tcm","CD4 Tem","CD4 Treg","CD4 Tex","CD8 Naive","CD8 Tem","CD8 CTL","CD8 Tex","CD4-CD8-","CD4+CD8+","T-mito","NK1","NK2","NK3","Naive","Memory","ASC","ABC","CD14","CD16","Intermed","cDC1","cDC2","pDC","pre-DC")
#df$Cluster=factor(df$Cluster,levels=level)
print(table(df$Cluster))
saveRDS(df,file.path(args$outdir,"agingScore.rds"))
p<-ggplot(df,aes(x=Cluster,y=GS,fill=Cluster))+
  #geom_violin()+geom_jitter()+theme_bw()+
  geom_violin()+theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15,angle=90))
pdf(file.path(args$outdir,"aging_Score.pdf"),width = 16,height = 16)
print(p)
dev.off()

