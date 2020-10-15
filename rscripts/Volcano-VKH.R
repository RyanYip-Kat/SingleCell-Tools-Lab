library(argparse)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggstatsplot)
library(ggsci)
library(Cairo)



outdir<-"Volcano/NK/Total"
cutoff=0.5

x_limit=5
y_limit=1e-500

if(!dir.exists(outdir)){
  dir.create(outdir,recursive=TRUE)
}

file<-"Status-Markers/NK/Total/_markers.csv"
markers<-read.csv(file,stringsAsFactors=FALSE)
##################
markers$logfoldchanges[markers$logfoldchanges > x_limit]=x_limit
markers$logfoldchanges[markers$logfoldchanges< -x_limit]= -x_limit

#### limit ayis
markers$pvals[markers$pvals < y_limit]= y_limit
#################
down_markers<-subset(markers,cluster=="HC" & logfoldchanges > 0)
down_markers$logfoldchanges<- -down_markers$logfoldchanges
up_markers<-subset(markers,cluster=="VKH" & logfoldchanges > 0)

markers<-rbind(up_markers,down_markers)
direction<-with(markers,ifelse(logfoldchanges>0,"Up","Down"))  # Up(ER),Down(HC)
markers$direction<-direction

label<-with(markers,ifelse(abs(logfoldchanges)>cutoff,names,""))

markers$label<-label
label_color<-with(markers,ifelse(label!="",ifelse(direction=="Up" ,"red","blue"),"grey"))
markers$label_color<-label_color

p <- ggplot(data =markers,
            aes(x =logfoldchanges,
                y = -log10(pvals),
                colour=label_color)) +
  geom_point(size=3.5) +
  scale_color_manual(values=c("#71B6A1", "grey","#1F77B4"))+
  xlim(c(-x_limit,x_limit))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential metabolites")+
  theme_half_open()+
  geom_text_repel(data = markers, aes(x = logfoldchanges,
                                      y = -log10(pvals),
                                      label = label),
                  size =2.5,box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = "grey50",
                  show.legend = TRUE,
                  colour = "black")+theme_bw()

#jpeg("volcano_plot.jpeg",width = 1024,height=1024)
#pdf(file.path(outdir,"volcano_plot.pdf"),width = 12,height = 16)
p+ theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
         legend.position="none",
         axis.title = element_text(size = 22,face = "bold"),
         axis.text = element_text(size=20,face = "bold"),
         legend.title = element_blank(),
         legend.text=element_text(size=20,face="bold"))

ggsave(file.path(outdir,"volcano_plot.pdf"),width = 16,height = 12,device =cairo_pdf)
dev.off()



