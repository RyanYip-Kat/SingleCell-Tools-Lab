library(ggplot2)
library(stringr)
library(reshape2)
library(cowplot)

DATA<-read.csv("E:/Data/Su/20200423/1-InteractionScores.csv",stringsAsFactors = FALSE,row.names = 1)
DATA<-as.data.frame(t(DATA))
status<-unlist(lapply(rownames(DATA),function(cell){
  return(str_split(cell,"\\.")[[1]][2])
}))

cell2cell<-unlist(lapply(rownames(DATA),function(cell){
  return(str_split(cell,"\\.")[[1]][1])
}))


DATA$status<-factor(status,levels = c("YH","AH","YCR","ACR"))
DATA$cell2cell<-cell2cell


scale.func <- switch(EXPR = scale.by, size = scale_size, 
                     radius = scale_radius, 
                     stop("'scale.by' must be either 'size' or 'radius'"))

n<-ncol(DATA)
lapply(unique(DATA$cell2cell),function(cell){
  x<-subset(DATA,cell2cell==cell)
  xx<-x
  x<-x[,-c(n-1,n)]
  x<-as.data.frame(x)
  idx<-as.character(read.table(file.path("InteractionPlot",paste0(cell,".txt")))$V1)
  x<-x[,idx]
  print(dim(x))
  x$status<-xx$status
  df<-melt(x)
  colnames(df)<-c("Status","LR","Score")
  df<-subset(df,!Status%in%c("YH","AH"))
  
  pdf(file.path("InteractionPlot",paste0(cell,"_InteractionScore.pdf")),width=8,height = 24)
  print(ggplot(df,aes(x=Status,y=LR,colour=Score,size=Score))+
          geom_point(alpha=2)+scale_radius(range = c(0,15))+
          scale_color_gradient2(low = "#0000EE",mid = "#F5F5F5",high = "#FA0000")+
          theme_bw()+guides( size = guide_legend(order = 1))+ # colour is the first legend
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.title = element_blank(),
                axis.text.y =element_text(size=20,face = "bold"),
                axis.text.x = element_text(size=20,face = "bold"),
                axis.title = element_blank()))
  dev.off()
})




#########
lapply(unique(DATA$cell2cell),function(cell){
  x<-subset(DATA,cell2cell==cell)
  xx<-x
  x<-x[,-c(106,107)]
  x[2,]<-x[2,]-x[1,]
  x[3,]<-x[3,]-x[1,]
  x[4,]<-x[4,]-x[1,]
  x[1,]<-0
  x<-as.data.frame(x)
  
  x$status<-xx$status
  df<-melt(x)
  colnames(df)<-c("Status","LR","Score")
  df<-subset(df,Status!="YH")
  df$normal.score<-scale(df$Score)
  
  pdf(paste0(cell,"_InteractionScore.pdf"),width=8,height = 24)
  print(ggplot(df,aes(x=Status,y=LR,size=Score,color=normal.score))+
          geom_point()+
          scale_color_gradient2(low = "#0000EE",mid = "#F5F5F5",high = "#FA0000")+
          theme_classic()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y =element_text(size=10,face = "bold"),
                axis.text.x = element_text(size=18),
                axis.title = element_blank()))
  dev.off()
})





