library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggstatsplot)

cutoff=2.0
x_limit=5.5
y_limit=1e-200

markers=read.csv("_markers.csv",stringsAsFactors = F)
markers<-subset(markers,pvals<0.05)

#### limit axis
markers$logfoldchanges[markers$logfoldchanges > x_limit]=x_limit
markers$logfoldchanges[markers$logfoldchanges< -x_limit]= -x_limit

#### limit ayis
markers$pvals[markers$pvals < y_limit]= y_limit

down_markers<-subset(markers,cluster=="HC" & logfoldchanges > 0)
down_markers$logfoldchanges<- -down_markers$logfoldchanges
up_markers<-subset(markers,cluster=="DR" & logfoldchanges > 0)


markers<-rbind(up_markers,down_markers)
#direction<-with(markers,ifelse(logfoldchanges>0,"Up","Down"))  # Up(ER),Down(HC)
direction<-with(markers,ifelse(logfoldchanges>cutoff,"Up",
                               ifelse(logfoldchanges < -cutoff,"Down","NS")))
markers$direction<-direction

###########
p1<-ggplot(data = markers,
       aes(x = logfoldchanges, 
           y = -log10(pvals),color=direction))+
  geom_point(aes(x=logfoldchanges,y=-log10(pvals)),size=2,shape=19)+
  scale_color_manual(values = c('skyblue','grey','red'))+
  theme_ggstatsplot()

##############

p2 <- ggplot(data = markers,
             aes(x = logfoldchanges, 
                 y = -log10(pvals))) +
  geom_point(size =1.5,aes(color = direction),show.legend = T) +
  scale_color_manual(values = c('skyblue', 'gray', 'pink')) +
  labs(x = 'Log2(fold change)', y = '-log10(p-value)') +  
  theme_ggstatsplot() +
  theme_half_open() +
  labs(x = 'Log2(fold change)', 
       y = '-log10(p-value)')+
  geom_hline(yintercept = -log10(0.05),  
             linetype = 'dotdash',
             color = 'grey30') + 
  geom_vline(xintercept = c(-cutoff, cutoff), 
             linetype = 'dotdash', 
             color = 'grey30')  +
  theme_half_open()

############## show
p3 <- p2+ geom_text_repel(data = subset(markers,abs(logfoldchanges)>3.5),
                          size = 3.0,
                          aes(label = names))


############# add specify genes
p4 <- ggplot(data = markers,
             aes(x = logfoldchanges, 
                 y = -log10(pvals))) +
  geom_point(size = 3,aes(color = direction),show.legend = T) +
  scale_color_manual(values = c('skyblue', 'gray', 'pink')) +
  labs(x = 'Log2(fold change)', y = '-log10(p-value)')+
  theme_ggstatsplot() +
  theme_half_open() +
  labs(x = 'Log2(fold change)', 
       y = '-log10(p-value)')+
  geom_hline(yintercept = -log10(0.05),  
             linetype = 'dotdash',
             color = 'grey30') + 
  geom_vline(xintercept = c(-cutoff, cutoff), 
             linetype = 'dotdash', 
             color = 'grey30')  +
  theme_half_open()+
  geom_text_repel(data = markers%>% 
                    filter(names%in% c('GSTM1',"TENM4","NACA","MYL6","CXCL2","MTX3","SOS1","CHMP7")),
                  size = 4.5,
                  aes(label = names))


#############
library(ggplot2)
library(GOplot)
library(reshape2)
library(stringr)

data("EC")
david <- EC$david
genelist <- EC$genelist
circ <- circle_dat(EC$david, EC$genelist)
chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

#######
DATA=read.csv("NPDR-DME-CHOOSE-METASCPE-ENRICH-INT.csv",stringsAsFactors = F)
df1<-subset(DATA,select = c(Category,GO,Description,Hits,Log.q.value.))
colnames(df1)<-c("Category","ID","Term","Genes","adj_pval")
df1$adj_pval<-10^(df1$adj_pval)
df1$Genes<-str_replace_all(df1$Genes,"\\|",",")

markers=read.csv("_markers.csv",stringsAsFactors = F)
df2<-subset(markers,select = c(names,logfoldchanges,scores,pvals,pvals_adj))
colnames(df2)<-c("ID","logFC","AveExpr","P.Value","adj.P.Val")

circ <- circle_dat(df1, df2)
GOBar(subset(circ, category == 'GO Biological Processes'))
GOBar(circ, display = 'multiple')
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))
GOBubble(circ, labels = 3)

GOBubble(circ, title = 'Bubble plot', 
         colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  

GOCircle(circ,label.size = 2.5)
##################
tmp=do.call(rbind,
            apply(df1, 1,function(x){
              data.frame(go=x[3],
                         gene=strsplit(x[4],'\\|')[[1]])
            })
)
tmp2=dcast(tmp,go~gene)
tmp2[is.na(tmp2)]=0
rownames(tmp2)=tmp2[,1]
tmp2=tmp2[,-1]
tmp2=t(tmp2)
tmp2[tmp2!=0]=1
tmp2=as.data.frame(tmp2)
tmp2$logFC=0
cg=rownames(tmp2)
tmp2=apply(tmp2,2,as.numeric)
rownames(tmp2)=cg

GOChord(tmp2, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 2.5)
GOBubble(tmp2, labels = 3)

