df<-data.frame("CD4-Naive"=431,
               "CD4-Tcm"=321,
               "CD4-Tem"=421,
               "CD4-Treg"=169,
               "CD4-Tex"=52)

df<-as.data.frame(t(df))
df$CellType=rownames(df)
colnames(df)<-c("nGene","CellType")
df$Bubble="Bubble"
#df$nGene<-as.integer(df$nGene)
df<-df[order(df$nGene,decreasing = T),]
df$CellType<-factor(df$CellType,levels = 
                      rev(c("CD4.Naive","CD4.Tcm","CD4.Tem","CD4.Treg","CD4.Tex")))

pdf("nGene_Bubble.pdf",width = 12,height = 12)
ggplot(data=df,aes(y=CellType,x=Bubble,size=nGene,color=CellType))+
  geom_point(alpha=0.8)+scale_size(range = c(20,30))+
  #geom_text(data=df,aes(x=Bubble,y=CellType,label=nGene))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y =element_text(size=28,face = "bold"),
        axis.text.x = element_blank(),
        axis.title = element_blank())

dev.off()
