library(argparse)
library(stringr)


root="E:/Data/Su/20200609/clontype/TCR_BCR"
path="AA_BCR"
clusters<-read.csv("E:/Data/Su/20200605/BC/cluster.csv",stringsAsFactors = F)

upaired_file=list.files(file.path(root,path),"vloupe-unpaired-clonotypes.csv",full.names = T)
paired_file=list.files(file.path(root,path),"vloupe-clonotypes.csv",full.names = T)
df_upaired=read.csv(upaired_file,stringsAsFactors = F)
df_paired<-read.csv(paired_file,stringsAsFactors = F)

celltypes<-unique(clusters$celltype)
for(cell in celltypes){
  target_barcode=subset(clusters,celltype==cell)$barcode
  idx<-c()
  freq<-c()
  for(i in 1:nrow(aa_upaired)){
    upaired_barcode=str_split(aa_upaired$barcodes[i],";")[[1]]
    if(all(upaired_barcode%in%target_barcode)){
      idx<-c(idx,i)
      freq=c(freq,length(upaired_barcode))
    }
  }
  df=df_upaired[idx,]
  df$n_freq<-freq
  clonotype_ids=df$clonotype_ids
  v_df=df_paired[df_paired$clonotype_id%in%clonotype_ids,]
  write.table(v_df,file.path("E:/Data/Su/20200605/BC",paste0(cell,"_AA_BCR_vloupe_clonetype.csv")),sep=",",quote = F,row.names = F)
  write.table(df,file.path("E:/Data/Su/20200605/BC",paste0(cell,"_AA_BCR.csv")),sep=",",quote = F,row.names = F)
}
