library(ggplot2)
library(monocle)

cds=readRDS("MSCE-Cluster1/c245-vst/cds.rds")
BEAM_res <- BEAM(cds, branch_point = 1, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
write.table(BEAM_res,"BEAM_res.csv",sep=",",quote=FALSE)
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf("plot_genes_branched_heatmap.pdf",width=8,height=24)
res=plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T,
			    return_heatmap=T)
ph_res=res[["ph_res"]]
labels=ph_res$tree_row$labels
annotation_row=as.integer(res[["annotation_row"]][["Cluster"]])

mat=res[["heatmap_matrix"]]
mat=mat[ph_res$tree_row$order,]
annotation_row=annotation_row[ph_res$tree_row$order]
df=data.frame("genes"=rownames(mat),"cluster"=annotation_row)

write.table(df,"heatmap_labels.csv",sep=",",quote=FALSE,row.names=F)
#saveRDS(res,"plot_genes_branched_heatmap.rds")
dev.off()
