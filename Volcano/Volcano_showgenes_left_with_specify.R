library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library(argparse)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

args<-parser$parse_args()


if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}



####
show_genes<-c("MALAT1","TGFBR2","FOS“,‘IL6R","JAK2","MALT1","CREB1","CX3CR1","ITGA4","CCAR1","RIF1","TXNIP","KLF6","TNFAIP8","PTPRC","CBL","FNTA","APPL1","ARF6")
###
markers<-read.csv("../20200420/output/aging-xxx/Status-Markers/DC/AS_YS_markers.csv",stringsAsFactors=FALSE)
markers$logFC<-log2(exp(markers$avg_logFC))
down_markers<-subset(markers,cluster=="YS" & logFC > 0)
down_markers$logFC<- -down_markers$logFC
up_markers<-subset(markers,cluster=="AS" & logFC > 0)

markers<-rbind(up_markers,down_markers)

direction<-with(markers,ifelse(logFC>0,"Up","Down"))  
markers$direction<-direction

label<-with(markers,ifelse(logFC > 0.5,
                           gene,ifelse(logFC < 0 & gene%in%show_genes,
                                       gene,"")))



#label<-with(markers,ifelse(abs(logFC)>0.5,gene,""))
markers$label<-label
label_color<-with(markers,ifelse(label!="",ifelse(direction=="Up" ,"red","blue"),"grey"))
markers$label_color<-label_color


p <- ggplot(data =markers,
            aes(x =logFC,
                y = -log10(p_val),
                colour=label_color)) +
  geom_point(size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-4.5, 4.5))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential metabolites")+
  geom_text_repel(data = markers, aes(x = logFC,
                                      y = -log10(p_val),
                                      label = label,
                                      colour=label_color),
                  size = 5,box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = "grey50",
                  show.legend = TRUE,
                  direction = "both")+theme_bw()

jpeg(file.path(args$outdir,"volcano_plot.jpeg"),width = 1024,height=1024)
p+ theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
         legend.position="none",
         axis.title = element_text(size = 22,face = "bold"),
         axis.text = element_text(size=20,face = "bold"),
         legend.title = element_blank(),
         legend.text=element_text(size=20,face="bold"))

dev.off()
