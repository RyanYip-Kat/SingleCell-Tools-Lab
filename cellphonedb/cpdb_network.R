library(psych)
library(purrr)
library(qgraph)
library(igraph)
library(argparse)
library(stringr)

#############   reference to  :  https://github.com/zktuong/ktplots/
parser <- ArgumentParser(description='A Program to plot cellphonebd count network...')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--network",
                    type="character",
                    default=NULL,
                    help="count network file")


args <- parser$parse_args()
###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
makedir(outDir)

###############################
message("INFO : Loading network ...")
mynet <- read.delim(args$network, check.names = FALSE)
mynet=subset(mynet,count!=0)
net<- graph_from_data_frame(mynet)

############################### Plot 1
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局

E(net)$width  <- E(net)$count/10  # 边点权重（粗细）


message("INFO : Plots ...")
pdf(file.path(outDir,"Plot1.pdf"),width=16,height=12)
plot(net, edge.arrow.size=.1,
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7)
dev.off()

############################### Plot 2
net2 <- net  # 复制一份备用

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
}  # 这波操作谁有更好的解决方案？

pdf(file.path(outDir,"Plot2.pdf"),width=16,height=12)
plot(net, edge.arrow.size=.1,
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7)
dev.off()

pdf(file.path(outDir,"Plot3.pdf"),width=16,height=12)
plot(net, edge.arrow.size=.1,
     edge.curved=0.2, # 只是调了这个参数
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7)

dev.off()
####################################### Plot 3
pdf(file.path(outDir,"Plot4.pdf"),width=24,height=12)
par(mfrow=c(2,as.integer(length(unique(mynet$SOURCE))/2)+1), mar=c(.3,.3,.3,.3))
for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2

  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  # 故技重施

  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]

  plot(net1, edge.arrow.size=.1,
       edge.curved=0.4,
       edge.label = E(net1)$count, # 绘制边的权重
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1
  )

}
dev.off()
