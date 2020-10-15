library(ggplot2)
library(stringr)

path="BCAA"
GO_AllLists=read.csv(file.path("GO_DotPlot",path,"GO_AllLists.csv"),stringsAsFactors = F,sep = ",")
SelectedGO=read.csv(file.path("GO_DotPlot",path,"HeatmapSelectedGO.csv"),sep = ",",stringsAsFactors = F)
SelectedGO<-SelectedGO[!duplicated(SelectedGO$GO),]
rownames(SelectedGO)=SelectedGO$GO
GO_Term=as.character(SelectedGO$GO)
for(cell in unique(GO_AllLists$GeneList)){
  df<-subset(GO_AllLists,GeneList==cell)
  df<-subset(df,GO%in%GO_Term)
  df<-subset(df,select = c(GO,X.GeneInGOAndHitList))
  not_in=GO_Term[!GO_Term%in%df$GO]
  if(length(not_in)!=0){
    df_not=data.frame("GO"=GO_Term[!GO_Term%in%df$GO],"X.GeneInGOAndHitList"=0)
    df<-rbind(df,df_not)
  }
  rownames(df)<-df$GO
  colnames(df)=c("GO",paste0("GeneInGOAndHitList_",str_replace(cell," ","\\.")))
  df<-df[rownames(SelectedGO),]
  SelectedGO<-cbind(SelectedGO,df[,2,drop=F])
}

celltype=unique(GO_AllLists$GeneList)
celltype<-str_replace(celltype," ","\\.")
#celltype<-str_replace(celltype,"-","\\.") #cDC2
DATA<-list()
for(i in 1:length(celltype)){
  Name=colnames(SelectedGO)[str_detect(colnames(SelectedGO),celltype[i])]
  logname=Name[str_detect(Name,"Log")]
  GOHit=Name[str_detect(Name,"GeneInGOAndHitList")]
  df<-SelectedGO[,c("GO","Description",logname,GOHit)]
  colnames(df)<-c("GO","Description","LogP","GeneInGOAndHitList")
  df$CellType=celltype[i]
  DATA[[i]]<-df
}

DATA=do.call(rbind,DATA)
# cDC2
DATA$CellType<-str_replace(DATA$CellType,"\\.","-")
DATA$CellType<-factor(DATA$CellType,levels = c("Naive-BC","Memory-BC","ASC"))
#DATA$CellType<-factor(DATA$CellType,levels = c("CD8-Naive","CD8-Tem","CD8-CTL","CD8-Tex"))
#DATA$CellType<-factor(DATA$CellType,levels = c("CD4-Naive","CD4-Tcm","CD4-Tem","CD4-Treg","CD4-Tex"))

pdf(file.path("GO_DotPlot",paste0(path,"_Go.pdf")),width = 8,height = 8)
#jpeg("GO-DotPlot.jpeg",width = 1960,height = 1960)
#DATA$Description<-factor(DATA$Description,levels = rev(DATA$Description[1:10]))
ggplot(DATA,aes(x=CellType,y=Description,colour=-LogP,size=GeneInGOAndHitList))+
  geom_point()+scale_size(range = c(1,8))+
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

