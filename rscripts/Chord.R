#加载包
library(statnet)
library(circlize)
#导入数据----
data<-read.csv("AA_BCR-vloupe-vj_gene_heatmap.csv",header=T,row.names =1)
mydata<-as.matrix(t(data))
#这里新建一个颜色，颜色都是提前选好的
grid.col = NULL
grid.col[rownames(mydata)]=sample(c("lavender","#85F29B","#FFA500","#9370DB", "#FFC0CB", "#EED5D2", 
                                    "#8B6969"),size = nrow(mydata),replace = T)

grid.col[colnames(mydata)]=rep("grey",ncol(mydata))

#开始画图----
circos.par(gap.degree=c(rep(2,nrow(mydata)-1),10, rep(2,ncol(mydata)-1),10),start.degree=180)
#定义画图环境
pdf("AA_BCR_Chord.pdf",height = 12,width = 12)
chordDiagram(mydata,directional=TRUE,
                            diffHeight = 0.06,
                            grid.col = grid.col, 
                            transparency = 0.5)
#circos.axis(labels = F,major.tick = FALSE)
#添加图例----
# legend("right",pch=20,legend=colnames(mydata),
#         col=grid.col[colnames(mydata)],bty="n",
#         cex=1,pt.cex=3,border="black")
dev.off()
circos.clear()
