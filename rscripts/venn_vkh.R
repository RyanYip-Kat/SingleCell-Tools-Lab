library(stringr)
split_status_DE<-function(DATA,outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE)
  }
  stopifnot("cluster"%in%colnames(DATA))
  clusters=unique(DATA$cluster)
  for(c in clusters){
    df=subset(DATA,cluster==c)
    filename=paste0(c,"_markers.csv")
    write.table(df,file.path(outdir,filename),sep=",",quote = FALSE,row.names = F)
  }
}


library(VennDiagram)
library(UpSetR)
library(stringr)

path="Status-Markers/Mono"
outdir<-"Venn/Mono"
cl="HC"

color=ifelse(cl=="HC","#1F77B4","red")
if(!dir.exists(outdir)){
  dir.create(outdir,recursive=TRUE)
}

#files=list.files(path,pattern = "YA_markers.csv",recursive = TRUE,full.names = T)
#files=files[str_detect(files,"Total")]
name=c("CD14","CD16","Intermed")
gene_list=list()
for(i in 1:length(name)){
  file=file.path(path,name[i],paste0(cl,"_markers.csv"))
  df=read.csv(file)
  genes=unique(df$name)
  gene_list[[name[i]]]<-genes
}


genes<-unique(as.character(unlist(gene_list)))
gene_table<-lapply(names(gene_list),function(name){
  gl=gene_list[[name]]
  gc<-c()
  for(gene in genes){
    if(gene%in%gl){
      gc<-c(gc,1)
    }else{
      gc<-c(gc,0)
    }
  }
  df<-data.frame(gc,row.names = genes)
  colnames(df)=name
  return(df)
})

DATA=do.call(cbind,gene_table)
#DATA=fromList(gene_list)
pdf(file.path(outdir,paste0(cl,"_Upset.pdf")),width = 10,height = 12)
upset(DATA, nsets = 7, nintersects = 30,
      sets=rev(name),
      keep.order = T,mb.ratio = c(0.5, 0.5),
      order.by = c("freq"), decreasing = c(TRUE,FALSE),
      point.size=2.5,sets.bar.color="blue",matrix.color=color,text.scale=2.0,
      main.bar.color="gray23")

dev.off()

write.table(DATA,file.path(outdir,paste0(cl,"_UpSet_table.csv")),sep=",",quote = F)
