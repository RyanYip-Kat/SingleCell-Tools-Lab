library(ggplot2)
library(stringr)
library(Cairo)

files<-list.files("E:/Data/Su/20200502",pattern = ".csv",full.names = TRUE)

lapply(files,function(file){
  DATA<-read.csv(file,sep = ",",stringsAsFactors = F,header = TRUE)
  DATA<-subset(DATA,select = c(GO,Description,LogP,X.GeneInGOAndHitList))
  colnames(DATA)<-c("GO","Description","LogP","GeneInGOAndHitList")
  DATA<-DATA[order(DATA$GeneInGOAndHitList,decreasing = TRUE),]
  DATA<-DATA[!duplicated(DATA$Description),]
  
  p<-ggplot(data=DATA,aes(y=reorder(Description,GeneInGOAndHitList),x=GeneInGOAndHitList, fill=-LogP))+
    geom_bar( stat='identity') + 
    geom_text(data=DATA,aes(y=reorder(Description,GeneInGOAndHitList),label=GO),hjust=0.25,size=4)+
    scale_x_reverse()+
    scale_y_discrete(position = "right")+
    #coord_flip() +
    scale_fill_gradient("-logP",low="grey", high ="red")+ #"#00BFFF") +
    ylab("") +
    xlab("Gene count") +
    theme_bw()+
    theme(
      axis.text.x=element_text(color="black",size=17),
      axis.text.y=element_text(color="black", size=17),
      axis.title.x = element_text(color="black", size=rel(1.6)),
      legend.text=element_text(color="black",size=rel(1.0)),
      legend.title = element_text(color="black",size=rel(1.1)),
      # legend.position=c(0,1),legend.justification=c(-1,0)
      legend.position="left",
    )
  filename<-paste0(str_split(basename(file),"\\.")[[1]][1],".pdf")
  pdf(filename,width = 16,height =8)
  print(p)
  dev.off()
})

