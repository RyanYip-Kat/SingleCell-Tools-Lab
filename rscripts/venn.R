library(VennDiagram)
library(UpSetR)
library(stringr)


df=read.csv("DC/DC_YA.csv",stringsAsFactors = F)
colnames(df)<-c("gene","celltype")
celltype<-unique(df$celltype)
gene_list=list()
for(cell in celltype){
  x<-subset(df,celltype==cell)
  print(dim(x))
  genes=unique(x$gene)
  gene_list[[cell]]<-genes
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
#pdf("CD4_Venn/YA_Upset.pdf",width = 16,height = 16)
pdf("DC/DC_YA_Upset.pdf",width = 6,height = 4.5)
upset(DATA, nsets = 7, nintersects = 30,
      #sets=rev(c("CD8 Naive","CD8 Tem","CD8 CTL","CD8 Tex")),
      sets=rev(c("cDC1","cDC2","pDC","pre-DC")),
      keep.order = T, mb.ratio = c(0.7, 0.3),
      order.by = c("freq"), decreasing = c(TRUE,FALSE),
      point.size=2.2,sets.bar.color="blue",matrix.color="#1F77B4",text.scale=1.2,
      main.bar.color="gray23")  # AA : red ,YA:#1F77B4

dev.off()

write.table(DATA,"DC_YA_table.csv",sep=",",quote = F)




