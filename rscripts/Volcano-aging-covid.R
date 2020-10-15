library(argparse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)


dbs<-c("KEGG_2016",
       "GO_Biological_Process_2018",
       "GO_Molecular_Function_2018",
       "GO_Cellular_Component_2018")

markers<-read.csv("Mono_ACR-YCR_markers.csv",
                  stringsAsFactors=FALSE)

down_markers<-subset(markers,cluster=="YCR" & logfoldchanges > 0)
down_markers$logfoldchanges<- -down_markers$logfoldchanges
up_markers<-subset(markers,cluster=="ACR" & logfoldchanges > 0)

markers<-rbind(up_markers,down_markers)
direction<-with(markers,ifelse(logfoldchanges>0,"Up","Down"))  # Up(ER),Down(HC)
markers$direction<-direction

show_genes<-c("IER2","ZFP36","CDKN1C","NAMPT",
              "CXCL2","IL1B","ATP2B1-AS1","NR4A1","JUNB",
              "DUSP2","CCL3","CXCL8","EGR1","JUND","HES1","HCAR2",
              "DUSP1","FOSB","RGCC","ARL5B","TMEM107","FOS",
              "CHMP1B","TRIB1","DUSP2","JUN","RGS1","TGFB1",
              "STAT1","TNFAIP2","IGF2R","TMEM176B","GBP5","RETN",
              "NCF1","IRF1","CD63","PDXK","MCEMP1","SULT1A2","LACTB","RAB5IF","GADD45B")

#label<-with(markers,ifelse(abs(logfoldchanges)>0.75,names,""))  # mono 1.25

label<-with(markers,ifelse(names%in%show_genes,names,
                           ifelse(abs(logfoldchanges)>0.75,names,"")))
#label<-with(markers,ifelse(abs(logfoldchanges)>0.5,gene,""))
markers$label<-label
label_color<-with(markers,ifelse(label!="",ifelse(direction=="Up" ,"red","blue"),"grey"))
markers$label_color<-label_color

p <- ggplot(data =markers,
            aes(x =logfoldchanges,
                y = -log10(pvals),
                colour=label_color))+
  geom_point(size=3.5) +
  scale_color_manual(values=c("#71B6A1", "grey","#1F77B4"))+
  xlim(c(-4.5, 4.5))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential metabolites")+
  geom_text_repel(data = markers, aes(x = logfoldchanges,
                                      y = -log10(pvals),
                                      label = label),
                  size = 5,box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = "grey50",
                  show.legend = TRUE,
                  colour = "black")+theme_bw()

#jpeg("volcano_plot.jpeg",width = 1024,height=1024)
pdf("mono_volcano_plot_1.pdf",width = 12,height = 16)
pdf("mono_volcano_plot_2.pdf",width = 16,height = 16)
p+ theme(plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
         legend.position="none",
         axis.title = element_text(size = 25,face = "bold"),
         axis.text = element_text(size=25,face = "bold"),
         legend.title = element_blank(),
         legend.text=element_text(size=25,face="bold"))

dev.off()
