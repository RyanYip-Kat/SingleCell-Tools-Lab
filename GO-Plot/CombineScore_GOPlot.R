Skip to content
Search or jump to…

Pull requests
Issues
Marketplace
Explore
 
@RyanYip-Kat 
yejg2017
/
CovID2019
1
00
Code
Issues
Pull requests
Actions
Projects
Wiki
Security
Insights
CovID2019/20200317/GO-Volcano.R
@yejg2017
yejg2017 CovID 2019 Paper
Latest commit c930282 on 23 Mar
 History
 1 contributor
169 lines (140 sloc)  5.73 KB
  
library(clusterProfiler)
library(pathview)
library(topGO)
library(AnnotationHub)
library(biomaRt)
library(Rgraphviz)
library(DOSE)
library(org.Hs.eg.db)
library(argparse)
library(dplyr)
library(enrichR)
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

args<-parser$parse_args()


if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

dbs<-c("KEGG_2016",
       "GO_Biological_Process_2018",
       "GO_Molecular_Function_2018",
       "GO_Cellular_Component_2018")

markers<-readRDS(args$marker)
#markers<-readRDS("CD4_HER_markers.rds")
markers$logFC<-log2(exp(markers$avg_logFC))
down_markers<-subset(markers,cluster=="HC" & logFC > 0)
down_markers$logFC<- -down_markers$logFC
up_markers<-subset(markers,cluster=="ER" & logFC > 0)

markers<-rbind(up_markers,down_markers)
direction<-with(markers,ifelse(logFC>0,"Up","Down"))  # Up(ER),Down(HC)
markers$direction<-direction

label<-with(markers,ifelse(abs(logFC)>0.75,gene,""))
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
  theme(plot.title = element_text(hjust = 0.5,size=20,face = "bold"),
        legend.position="right",
        axis.title = element_text(size = 20,face = "bold"),
        axis.text = element_text(size=18,face = "bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=18,face="bold"))

jpeg(file.path(args$outdir,"volcano_plot.jpeg"),width = 1024,height=1024)
p+geom_text_repel(data = markers, aes(x = logFC,
                                         y = -log10(p_val),
                                         label = label),
                     size = 3,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",
                     show.legend = TRUE,
                     colour = "black")+theme_bw()+theme(legend.position = "none")


dev.off()



print("### Pathway analysis")
up_genes<-up_markers$gene
up_enriched <- enrichr(up_genes, dbs)
down_genes<-down_markers$gene
down_enriched<-enrichr(down_genes, dbs)
saveRDS(up_enriched,file.path(args$outdir,"ER_enriched.rds"))
saveRDS(down_enriched,file.path(args$outdir,"HC_enriched.rds"))

up_result<-up_enriched$KEGG_2016
up_result$direction<-"ER"

write.table(up_result,file.path(args$outdir,"ER_KEGG_table.csv"),sep=",",quote = FALSE,row.names = FALSE)
up_result<-subset(up_result,P.value<0.05)
up_result<-up_result%>%top_n(30,wt=Combined.Score)

down_result<-down_enriched$KEGG_2016
down_result$direction<-"HC"
write.table(down_result,file.path(args$outdir,"HC_KEGG_table.csv"),sep=",",quote = FALSE,row.names = FALSE)
down_result<-subset(down_result,P.value<0.05)
down_result<-down_result%>%top_n(30,wt=Combined.Score)

df<-rbind(up_result,down_result)
df$Term<-str_replace_all(df$Term,"Homo sapiens hsa\\d+","")
df$qvalue<-factor(df$P.value,levels = unique(df$P.value))
df<-df[!duplicated(df$qvalue),]


p<-ggplot(data=df,aes(x=qvalue,y=Combined.Score,fill=direction))+
  geom_bar(stat="identity")+
  geom_text(data=df,aes(x=qvalue,label=Term),hjust=-0.1,size=2.5)+
  facet_wrap(~direction,scales = "free",ncol = 1)+ylim(0,500)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,size=12),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()

pdf(file.path(args$outdir,"HER_KEGG.pdf"),width=12,height=8)
print(p)
dev.off()


up_result<-up_enriched$GO_Biological_Process_2018
up_result$direction<-"ER"
write.table(up_result,file.path(args$outdir,"ER_BP_table.csv"),sep=",",quote = FALSE,row.names = FALSE)
up_result<-subset(up_result,P.value<0.05)
up_result<-up_result%>%top_n(30,wt=Combined.Score)

down_result<-down_enriched$GO_Biological_Process_2018
down_result$direction<-"HC"
write.table(down_result,file.path(args$outdir,"HC_BP_table.csv"),sep=",",quote = FALSE,row.names = FALSE)
down_result<-subset(down_result,P.value<0.05)
down_result<-down_result%>%top_n(30,wt=Combined.Score)

df<-rbind(up_result,down_result)
df$Term<-str_replace_all(df$Term,"Homo sapiens hsa\\d+","")
df$qvalue<-factor(df$P.value,levels = unique(df$P.value))
df<-df[!duplicated(df$qvalue),]

p<-ggplot(data=df,aes(x=qvalue,y=Combined.Score,fill=direction))+
  geom_bar(stat="identity")+
  geom_text(data=df,aes(x=qvalue,label=Term),hjust=-0.1,size=2.5)+
  facet_wrap(~direction,scales = "free",ncol = 1)+ylim(0,800)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,size=12),
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()

pdf(file.path(args$outdir,"HER_BP.pdf"),width=12,height=8)
print(p)
dev.off()
© 2020 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
