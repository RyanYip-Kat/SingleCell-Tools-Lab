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
                    filter(names%in% c('GSTM1',"TENM4")),
                  size = 4.5,
                  aes(label = names))

