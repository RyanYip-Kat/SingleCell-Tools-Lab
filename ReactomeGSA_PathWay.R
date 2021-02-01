library(argparse)
library(stringr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ReactomeGSA.data)
library(ReactomeGSA)

parser <- ArgumentParser(description='ReactomeGSA Pathway Analysis')
parser$add_argument("--seurat",
                    type="character",
		    default=NULL,
                    help="the path of project saved")

parser$add_argument("--outdir",
                    type="character",
                    default="ReactomeGSA_Result")


parser$add_argument("--useIdent",
		    type="character",
		    default="seurat_clusters",
		    help="which group as idents in seurat object")
args <- parser$parse_args()

############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

###########################
outDir=args$outdir
useIdent=args$useIdent
makedir(outDir)

email="ryanyip_@hotmail.com"

##########################
message("INFO : Loading dataset ...")
seurat=readRDS(args$seurat)
metadata=seurat@meta.data
Idents(seurat)=metadata[[useIdent]]

message("INFO : ReactomeGSA analyse clusters ...")
gsva_result <- analyse_sc_clusters(seurat, 
				   verbose = TRUE,
				   assay = "RNA",
				   slot="counts",
				   report_email=email,
				   create_reports=TRUE)


saveRDS(gsva_result,file.path(outDir,"gsva_result.rds"))
message("INFO : Pathway analyse ...")
pathway_expression <- pathways(gsva_result)
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
saveRDS(pathway_expression,file.path(outDir,"pathway_expression.rds"))

message("INFO : find the maximum differently expressed pathway ...")
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
   values <- as.numeric(row[2:length(row)])
   return(data.frame(name = row[1], min = min(values), max = max(values)))
 }))

max_difference$diff <- max_difference$max - max_difference$min
#sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

message("INFO : Visualize Pathway ...")
p1 <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1]) 
p1$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> df

p1=ggplot(data=df,aes(cluster_id,   expr ,fill=absmy    ))+
  geom_bar(stat='identity') + theme_bw()+
  theme(axis.text.x = element_text(angle =45,hjust = .9,size = 10,vjust = 0.9))+
  ggtitle(max_difference$name[1]) + theme(legend.position="none")

ggsave(file.path(outDir,"PathwayIdentBar.pdf"),plot=p1,width=12,height=10)

pdf(file.path(outDir,"gsva_heatmap.pdf"),width=8,height=16)
plot_gsva_heatmap(gsva_result,
		  pathway_ids=NULL,  #show rownames,A vector of pathway ids if not NULL 
		  truncate_names=TRUE, #If set, long pathway names are truncated
		  max_pathways = 15, 
		  margins = c(6,20))
dev.off()

p=plot_gsva_pca(gsva_result)
ggsave(file.path(outDir,"gsva_pca.pdf"),plot=p,width=12,height=10)
