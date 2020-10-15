library(ggplot2)
library(stringr)

total=read.csv("label_fine_enrich.csv",stringsAsFactors = F)
DATA=read.csv("monocoyte-go-select.csv",stringsAsFactors = F)
DATA=total[total$native%in%DATA$native,]

DATA=subset(DATA,select = c(native,p_value,name,intersection_size,cluster))
DATA$LogP=log10(DATA$p_value)
colnames(DATA)=c("GO","p.value","Description","intersection_size","CellType","LogP")

DATA$CellType<-factor(DATA$CellType,levels = c("CD14","CD16","Intermed"))
pdf("Mono_Go.pdf",width = 10,height = 16)
ggplot(DATA,aes(x=CellType,y=Description,colour=-LogP,size=intersection_size))+
  geom_point()+scale_size(range = c(1,10))+
  #scale_colour_gradient(low ="#F0FFFF",high  = "#00008B")+   # blue
  scale_colour_gradient(low ="#FFE4B5",high  = "#FF4500")+  # red
  theme_bw()+guides( size = guide_legend(order = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.title = element_blank(),
        axis.text.y =element_text(size=12,face = "bold"),
        axis.text.x = element_text(size=12,face = "bold",angle = 45),
        axis.title = element_blank())


dev.off()





##################
GO_AllLists=read.csv("label_fine_enrich.csv",stringsAsFactors = F)
SelectedGO=read.csv("monocoyte-go-select.csv",stringsAsFactors = F)
#DATA=total[total$native%in%DATA$native,]

SelectedGO=subset(SelectedGO,select = c(native,p_value,name,intersection_size,cluster))
SelectedGO$LogP=log10(SelectedGO$p_value)
colnames(SelectedGO)=c("GO","p.value","Description","intersection_size","CellType","LogP")

SelectedGO<-SelectedGO[!duplicated(SelectedGO$GO),]
rownames(SelectedGO)=SelectedGO$GO
GO_Term=as.character(SelectedGO$GO)
for(cell in unique(GO_AllLists$cluster)){
  df<-subset(GO_AllLists,cluster==cell)
  df<-subset(df,native%in%GO_Term)
  df<-subset(df,select = c(native,intersection_size))
  not_in=GO_Term[!GO_Term%in%df$native]
  if(length(not_in)!=0){
    df_not=data.frame("native"=GO_Term[!GO_Term%in%df$native],"intersection_size"=0)
    df<-rbind(df,df_not)
  }
  rownames(df)<-df$native
  colnames(df)=c("GO",paste0(cell,"_","GeneInGOAndHitList"))
  df<-df[rownames(SelectedGO),]
  SelectedGO<-cbind(SelectedGO,df[,2,drop=F])
}

celltype=unique(GO_AllLists$cluster)
DATA<-list()
for(i in 1:length(celltype)){
  Name=colnames(SelectedGO)[str_detect(colnames(SelectedGO),celltype[i])]
  df<-SelectedGO[,c("GO","Description","LogP",Name)]
  colnames(df)<-c("GO","Description","LogP","GeneInGOAndHitList")
  df$CellType=celltype[i]
  DATA[[i]]<-df
}

DATA=do.call(rbind,DATA)
DATA$CellType<-factor(DATA$CellType,levels = c("CD14","CD16","Intermed"))
pdf("Mono_Go.pdf",width = 10,height = 16)
ggplot(DATA,aes(x=CellType,y=Description,colour=-LogP,size=GeneInGOAndHitList))+
  geom_point()+scale_size(range = c(1,10))+
  #scale_colour_gradient(low ="#F0FFFF",high  = "#00008B")+   # blue
  scale_colour_gradient(low ="#FFE4B5",high  = "#FF4500")+  # red
  theme_bw()+guides( size = guide_legend(order = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.title = element_blank(),
        axis.text.y =element_text(size=12,face = "bold"),
        axis.text.x = element_text(size=12,face = "bold",angle = 45),
        axis.title = element_blank())


dev.off()
