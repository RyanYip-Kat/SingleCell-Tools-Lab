library(ggplot2)
library(stringr)
library(ggpubr)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 12,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
 
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)',colors=my_palette,limits=c(-6,0)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}

cellphonedb_data_input<-
  function(selected_rows = NULL,
           selected_columns = NULL,
           filename = 'plot.pdf',
           width = 8,
           height = 10,
           means_path = './means.txt',
           pvalues_path = './pvalues.txt',
           means_separator = '\t',
           pvalues_separator = '\t',
           output_extension = '.pdf'
  ){
    
    all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
    all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
    
    intr_pairs = all_pval$interacting_pair
    all_pval = all_pval[,-c(1:11)]
    all_means = all_means[,-c(1:11)]
    
    if(is.null(selected_rows)){
      selected_rows = intr_pairs
    }
    
    if(is.null(selected_columns)){
      selected_columns = colnames(all_pval)
    }
    
    sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
    sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
    
    df_names = expand.grid(selected_rows, selected_columns)
    pval = unlist(sel_pval)
    pval[pval==0] = 0.0009
    plot.data = cbind(df_names,pval)
    pr = unlist(as.data.frame(sel_means))
    pr[pr==0] = 1
    plot.data = cbind(plot.data,log2(pr))
    colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
    return(plot.data)
  
    }




my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

# write.table(hc,"hc_cellphone.csv",row.names = T,sep=",",quote = F)
# write.table(vkh,"vkh_cellphone.csv",row.names = T,sep=",",quote = F)


cols1=c("Mono|TC","Mono|BC","Mono|NovelBC")
cols2=c("DC|TC","DC|BC","DC|NovelBC")
cols3=c("TC|Mono","TC|DC")
cols4=c("BC|Mono","BC|DC","NovelBC|Mono","NovelBC|DC")
cols5=c("TC|TC","Mono|Mono")

cols6=c("TC|Mono","Mono|TC")
rows=as.character(read.csv("cellphoedb-vkh/TC-MONO.csv",header = TRUE,stringsAsFactors = F)[,1])

h=cellphonedb_data_input(means_path = "cellphone_vkh/HC/out/means.txt",
                         pvalues_path ="cellphone_vkh/HC/out/pvalues.txt")

v=cellphonedb_data_input(means_path = "cellphone_vkh/VKh/out/means.txt",
                         pvalues_path ="cellphone_vkh/VKh/out/pvalues.txt")

h=subset(h,pvalue<1)
v=subset(v,pvalue<1)

write.table(h,"cellphone_vkh/hc_cellphone.csv",row.names = T,sep=",",quote = F)
write.table(v,"cellphone_vkh/vkh_cellphone.csv",row.names = T,sep=",",quote = F)
p1=ggplot(h,aes(x=clusters,y=pair,size=-log10(pvalue),color=mean)) +
  geom_point()+scale_size(range = c(0,2))+
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text( hjust = 1,angle = 45),
        axis.text.y = element_text(size=5, colour = "black"),
        axis.title=element_blank(),
        #legend.background = element_rect( size = 0.5),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

p2=ggplot(v,aes(x=clusters,y=pair,size=-log10(pvalue),color=mean)) +
  geom_point()+scale_size(range = c(0,2))+
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text( hjust = 1,angle = 45),
        axis.text.y = element_text(size=5, colour = "black"),
        axis.title=element_blank(),
        #legend.background = element_rect( size = 0.5),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

ggarrange(p1,p2,common.legend = T,legend="right")
ggsave("cellphone_vkh/plot.pdf", width = 8, height = 24, device = cairo_pdf, limitsize=F)
dev.off()

p2
ggsave("cellphone_vkh/vkh.pdf", width = 8, height = 24, device = cairo_pdf, limitsize=F)
dev.off()

df=read.csv("cellphoedb-vkh/TC-MONO_cellphone3.csv",stringsAsFactors = F)
colnames(df)<-c("pair","mean","cluster")
p=ggplot(df,aes(x=cluster,y=pair,size=abs(mean),color=mean)) +
  geom_point()+scale_size(range = c(0,5))+
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text( hjust = 1,angle = 45),
        axis.text.y = element_text(size=5, colour = "black"),
        axis.title=element_blank(),
        #legend.background = element_rect( size = 0.5),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

pdf("cellphoedb-vkh/TC-MONO.pdf", width = 6, height = 24)
print(p)
dev.off()


files=list.files("cellphone_vkh","*.csv",full.names = T,recursive = T)
lapply(files,function(file){
  Data=read.csv(file,stringsAsFactors = F)
  colnames(Data)=c("pair","mean","cluster")
  p=ggplot(Data,aes(x=cluster,y=pair,size=mean,color=mean)) +
    geom_point()+scale_size(range = c(0,5))+
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=10, colour = "black"),
          axis.text.x = element_text( hjust = 1,angle = 45),
          axis.text.y = element_text(size=5, colour = "black"),
          axis.title=element_blank(),
          #legend.background = element_rect( size = 0.5),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  filename=str_split(basename(file),"\\.")[[1]][1]
  filename=file.path(dirname(file),paste0(filename,".pdf"))
  pdf(filename, width = 6, height = 18)
  print(p)
  dev.off()
  
})

############ heatmap
library(pheatmap)
heatmaps_plot_data = function(meta_file, pvalues_file, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = F,border_color='white', cluster_rows = F, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    # print(head(count_matrix))
    # 
    # all_sum = rowSums(count_matrix)
    # all_sum = cbind(names(all_sum), all_sum)
    # 
    # col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    # 
    # p<-pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
    #          border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
    #          main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col,
    #          )
    # return(p)
    return(count_matrix)
  }
  else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}

library(ComplexHeatmap)
library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)

col1 = "dodgerblue4";col2 = 'peachpuff';col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
pmat1=heatmaps_plot_data(meta_file = "HC/cell_meta.txt",
                pvalues_file = "HC/out/pvalues.txt")

pmat2=heatmaps_plot_data(meta_file = "VKH/cell_meta.txt",
                 pvalues_file = "VKH/out/pvalues.txt")


write.table(pmat1,"heatmap_count_hc.csv",sep=",",quote = F)
write.table(pmat2,"heatmap_count_vkh.csv",sep=",",quote = F)
p1=Heatmap(pmat1,col = col.heatmap,name = "hc",cluster_rows = F,cluster_columns = F,)
p2=Heatmap(pmat2,col = col.heatmap,name = "vkh",cluster_rows = F,cluster_columns = F)

p1 <- as.ggplot(p1)+scale_fill_continuous()
p2 <- as.ggplot(p2)
pdf("heatmap.pdf", width = 32, height = 16)
#p1+p2
ggarrange(p1,p2,common.legend = T,labels = )
dev.off()


hc=read.csv("heatmap_count_hc.csv",stringsAsFactors = F)
colnames(hc)=rownames(hc)
colnames(hc)=paste("HC",colnames(hc),sep = "_")
vkh=read.csv("heatmap_count_vkh.csv",stringsAsFactors = F)
colnames(vkh)=rownames(vkh)
colnames(vkh)=paste("VKH",colnames(vkh),sep = "_")
mat=cbind(hc,vkh)

mat=read.csv("heatmap_count.csv",stringsAsFactors = F,row.names = 1)
pdf("heatmap.pdf", width = 16, height = 8)
cols=c()
for(i in 1:length(colnames(vkh))){
  cols<-c(cols,c(colnames(hc)[i],colnames(vkh)[i]))
}
pheatmap(mat,
         cluster_rows = F,
         cluster_cols = F,
         color = col.heatmap)
         #gaps_col =ncol(hc))

dev.off()
