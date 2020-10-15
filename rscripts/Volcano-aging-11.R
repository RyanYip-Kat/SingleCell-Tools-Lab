library(argparse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--markers",
                    type="character",
                    #nargs="+",
                    default="")

parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--cutoff",
                    type="double",
                    default="0.75",
                    help="")
args<-parser$parse_args()


if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

dbs<-c("KEGG_2016",
       "GO_Biological_Process_2018",
       "GO_Molecular_Function_2018",
       "GO_Cellular_Component_2018")

markers<-read.csv(args$marker,stringsAsFactors=FALSE)
down_markers<-subset(markers,cluster=="YA" & logfoldchanges > 0)
down_markers$logfoldchanges<- -down_markers$logfoldchanges
up_markers<-subset(markers,cluster=="AA" & logfoldchanges > 0)

markers<-rbind(up_markers,down_markers)
direction<-with(markers,ifelse(logfoldchanges>0,"Up","Down"))  # Up(ER),Down(HC)
markers$direction<-direction

label<-with(markers,ifelse(abs(logfoldchanges)> args$cutoff,names,""))
#label<-with(markers,ifelse(abs(logfoldchanges)>0.5,gene,""))
markers$label<-label
label_color<-with(markers,ifelse(label!="",ifelse(direction=="Up" ,"red","blue"),"grey"))
markers$label_color<-label_color

p <- ggplot(data =markers,
            aes(x =logfoldchanges,
                y = -log10(pvals),
                colour=label_color)) +
  geom_point(size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
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
pdf(file.path(args$outdir,"volcano_plot.pdf"),width = 12,height = 16)
p+ theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
         legend.position="none",
         axis.title = element_text(size = 22,face = "bold"),
         axis.text = element_text(size=20,face = "bold"),
         legend.title = element_blank(),
         legend.text=element_text(size=20,face="bold"))

dev.off()

jpeg(file.path(args$outdir,"volcano_plot.jpeg"),width = 1024,height=1024)
p+ theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
         legend.position="none",
         axis.title = element_text(size = 22,face = "bold"),
         axis.text = element_text(size=20,face = "bold"),
         legend.title = element_blank(),
         legend.text=element_text(size=20,face="bold"))

dev.off()
