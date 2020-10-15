library(ggplot2)
library(stringr)
#DATA=read.csv("AA_GO.csv",stringsAsFactors =F,sep=",")


name=c("TC","NK","BC","Mono","DC")
files=c("E:/Data/Su/20200519/YA_GO/TC_FINAL_GO.csv",
        "E:/Data/Su/20200519/YA_GO/NK_FINAL_GO.csv",
        "E:/Data/Su/20200519/YA_GO/BC_FINAL_GO.csv",
        "E:/Data/Su/20200519/YA_GO/MONO_FINAL_GO.csv",
        "E:/Data/Su/20200519/YA_GO/DC_FINAL_GO.csv")

GOTerm=c("GO:0009611",
     "GO:0002683",
     "GO:0007162",
     "GO:0071559",
     "GO:0032733",
     "GO:0008380",
     "R-HSA-1169408",
     "GO:0046677",
     "GO:0006302",
     "GO:0002366"
)

df=lapply(1:length(name),function(i){
  d=read.csv(files[i],sep = ",",stringsAsFactors = F)
  d$CellType=name[i]
  rownames(d)=d$GO
  d<-d[,c("GO","Description","LogP","X.GeneInGOAndHitList","CellType")]
  d<-d[GOTerm,]
  return(d)
})

desc=c("response to wounding",
        "negative regulation of immune system process",
        "negative regulation of cell adhesion",
        "response to transforming growth factor beta",
        "positive regulation of interleukin-10 production",
        "RNA splicing",
        "ISG15 antiviral mechanism",
        "response to antibiotic",
        "double-strand break repair",
        "leukocyte activation involved in immune response"
)

df=do.call(rbind,df)
df<-as.data.frame(na.omit(df))
colnames(df)=c("GO","Description","LogP","GeneInGOAndHitList","CellType")
DATA=df
#colnames(DATA)<-c("GO","Description","LogP","GeneInGOAndHitList","CellType")
DATA$CellType<-factor(DATA$CellType,levels = c("TC","NK","BC","Mono","DC"))

pdf("YA-GO-DotPlot.pdf",width = 16,height = 16,pointsize = 15)
#jpeg("GO-DotPlot.jpeg",width = 1960,height = 1960)
ggplot(DATA,aes(x=CellType,y=Description,colour=LogP,size=GeneInGOAndHitList))+
  geom_point()+scale_size(range = c(1,15))+
  scale_colour_gradient2(low="blue",mid="white",high = "red")+
  theme_bw()+guides( size = guide_legend(order = 1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.title = element_blank(),
        axis.text.y =element_text(size=20,face = "bold"),
        axis.text.x = element_text(size=20,face = "bold"),
        axis.title = element_blank())


dev.off()
