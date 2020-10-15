library(ggplot2)
library(RColorBrewer)
library(Cairo)
library(showtext)

mydata<- read.csv("Bar_data.csv",stringsAsFactors=FALSE)
mydata=mydata[!duplicated(mydata$Description),]
up_mydata<-subset(mydata,cluster=="CD8")
down_mydata=subset(mydata,cluster=="CD4")
down_mydata$Enrichment= -down_mydata$Enrichment
mydata<-rbind(up_mydata,down_mydata)

mydata$Description<-as.character(mydata$Description)
mydata<-transform(mydata, label1=ifelse(Enrichment>=0,Description, NA),
                  
                  label2=ifelse(Enrichment>0, NA,Description))
mydata$Description <- factor(mydata$Description, levels = mydata$Description[order(mydata$Enrichment)])
ggplot(data = mydata, aes(x =Description, y =Enrichment ,fill = Enrichment)) +
  
  geom_bar(stat = "identity", width = 0.8,colour="black",size=0.25)+
  
  scale_fill_gradient2(low=brewer.pal(7,"Set1")[2],mid="grey90",high=brewer.pal(7,"Set1")[1],midpoint=0)+
  
  geom_text(aes(y = 0,     label=label2),size=3,hjust=-0.1)+ #添加负值部分的数据标签
  
  geom_text(aes(y = -0.001,label=label1),size=3,hjust= 1.1)+ #添加正值部分的数据标签
  
  coord_flip() +   #坐标轴翻转
  
  ylim(-5,5)+
  
  theme_minimal() + #图表主题设定
  
  theme(
    
    panel.grid.major=element_blank(),
    
    panel.grid.minor=element_blank(),
    
    panel.grid.major.x = element_line(colour = "grey80",size=.25),
    
    panel.grid.minor.x = element_line(colour = "grey80",size=.25),
    
    plot.title=element_text(size=15,hjust=.5),
    
    axis.text.x = element_text(face="plain", color="black",
                               
                               size=11, angle=0),
    
    axis.text.y = element_blank(),
    legend.position="right",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10))

#showtext.end()

#dev.off()

