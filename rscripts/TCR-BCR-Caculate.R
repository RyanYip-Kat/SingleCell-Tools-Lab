library(stringr)
outdir="E:/Data/Su/20200609/clontype/result" 
root="E:/Data/Su/20200609/clontype/TCR_BCR"

#path=paste0(ss,"_TCR")

cluster_file<-"E:/Data/Su/20200609/clontype/clusters"
celltype<-"Total"

clusters<-read.csv(file.path(cluster_file,celltype,"cluster.csv"),stringsAsFactors = F)
clusters$ident=unlist(lapply(clusters$barcode,function(b){return(str_split(b,"-")[[1]][2])}))
status=with(clusters,ifelse(as.integer(ident)%in%c(1:8),"AA","YA"))
clusters$status=status



##########
library(data.table)
DATA=fread("aging-11-YA-TCR/filtered_contig_annotations.csv")
ss="YA"
DATA=DATA[!duplicated(DATA$barcode),]
barcode<-unlist(lapply(DATA$barcode,function(b){
  x<-str_split(b,"-")[[1]]
  xx=paste0(x[1],"-",as.integer(x[2])+8)
  return(xx)
}))
if(ss=="YA"){
  DATA$barcode=barcode
}


# Total

clusters$celltype="Total"
celltypes<-unique(clusters$celltype)
idents<-unique(clusters$ident)
for(id in idents){
  target_barcode=subset(clusters,ident==id & status==ss)$barcode
  idx<-DATA$barcode%in%target_barcode
  df=DATA[idx,]
  
  clonotype_ids<-unique(df$raw_clonotype_id)
  freqs=c()
  clonotypes<-c()
  for(cid in clonotype_ids ){
    freq=sum(DATA$raw_clonotype_id%in%cid)
    freqs<-c(freqs,freq)
    clonotypes<-c(clonotypes,cid)
  }
  m<-data.frame("freqs"=freqs,"clonotypes"=clonotypes)
  res<-file.path(outdir,celltype,ss)
  if(!file.exists(res)){
    dir.create(res,recursive = T)
  }
  
  write.table(m,file.path(res,paste0(id,"_freqs_clonetype_total.csv")),
              sep=",",quote = F,row.names = F)
}



# subset
celltypes<-unique(clusters$celltype)
idents<-unique(clusters$ident)
for(id in idents){
  mat=subset(clusters,ident==id)
  for(cell in celltypes){
    target_barcode=subset(mat,celltype==cell & status==ss)$barcode
    idx<-DATA$barcode%in%target_barcode
    df=DATA[idx,]
    
    clonotype_ids<-unique(df$raw_clonotype_id)
    freqs=c()
    clonotypes<-c()
    for(cid in clonotype_ids ){
      freq=sum(DATA$raw_clonotype_id%in%cid)
      freqs<-c(freqs,freq)
      clonotypes<-c(clonotypes,cid)
    }
    m<-data.frame("freqs"=freqs,"clonotypes"=clonotypes)
    res<-file.path(outdir,celltype,"TCR",cell,ss,id)
    if(!file.exists(res)){
      dir.create(res,recursive = T)
    }
    
    write.table(m,file.path(res,"freqs_clonetype.csv"),
                sep=",",quote = F,row.names = F)
  }
}







