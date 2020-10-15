library(ggplot2)
library(stringr)
library(reshape2)
library(cowplot)

path="E:/Data/Su/20200617/Interaction/V3"
ACR=read.csv(file.path(path,"ACR","scores.csv"),row.names = 1)
YCR=read.csv(file.path(path,"YCR","scores.csv"),row.names = 1)
cell2cell=colnames(ACR)
for(cell in cell2cell){
  x<-ACR[,cell,drop=F]
  x$status="ACR"
  x$interaction=rownames(ACR)
  y<-YCR[,cell,drop=F]
  y$status="YCR"
  y$interaction=rownames(YCR)
  df<-rbind(x,y)
  #rownames(df)<-NULL
  colnames(df)<-c("Score","Status","LR")
  print(colnames(df))
  filename=file.path(path,paste0(cell,"_interaction.pdf"))#,width = 8,height = 64)
  ggplot(df,aes(x=Status,y=LR,colour=Score,size=Score))+
    geom_point(alpha=2)+scale_radius(range = c(0,15))+
    scale_color_gradient2(low = "#0000EE",mid = "#F5F5F5",high = "#FA0000")+
    theme_bw()+guides( size = guide_legend(order = 1))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.text.y =element_text(size=20,face = "bold"),
          axis.text.x = element_text(size=20,face = "bold"),
          axis.title = element_blank())
  ggsave(filename = filename,width = 8,height = 64,dpi = 250,limitsize = F,device = "pdf")
}



