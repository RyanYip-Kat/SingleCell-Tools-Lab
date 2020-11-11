library(pheatmap)
library(Seurat)
color=colorRampPalette(c("navy","#40E0D0","#B0171F"))(50)
pheatmap.idents<-function(seurat,features,slot="scale.data"){
        cts <- GetAssayData(seurat, slot = slot)
        if(slot=="counts"){
                cts <- log10(cts + 1)
        }
        features<-features[features%in%rownames(cts)]
        new_cluster <- sort(seurat$celltype)
        cts <- as.matrix(cts[features, names(new_cluster)])
        ac=data.frame(cluster=new_cluster)
        #color<-c("#40E0D0","#B0171F")
	n=length(unique(new_cluster))-1
	gap_cols<-cumsum(as.integer(table(new_cluster)))[1:n]
        p<-pheatmap(cts,color=color,show_colnames =F,show_rownames = T,
		    gaps_col=gap_cols,
                    cluster_rows = F,
                    cluster_cols = F,
		    fontsize_row=10,
                    fontsize=8,
                    annotation_col=ac)
        return(p)
}

#genes=topn$gene
#jpeg(file.path(args$outdir,"heatmap1.jpeg"),width=1024,height=1024)
#pheatmap.idents(seurat,genes)
#dev.off()

#jpeg(file.path(args$outdir,"heatmap2.jpeg"),width=1024,height=1024)
#pheatmap.idents(seurat,genes,"data")
#dev.off()

