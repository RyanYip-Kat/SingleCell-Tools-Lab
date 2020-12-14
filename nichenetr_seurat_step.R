library(stringr)
library(ggplot2)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(argparse)
library(future)
#############################
parser <- ArgumentParser(description='A Program to Perform NicheNet analysis starting from a Seurat object ')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat object")


parser$add_argument("--outdir",
                    type="character",
                    default="Result")


parser$add_argument("--condition",
		    type="character",
		    default=NULL,
		    help="condition colname in Seurat metadata to  explain differential expression between two conditions,length must be two!!!")

parser$add_argument("--condition_oi",
                    type="character",
                    default=NULL,
		    help="condition of interest name")

parser$add_argument("--condition_reference",
                    type="character",
                    default=NULL,
                    help="condition of reference name")


parser$add_argument("--receiver",
                    type="character",
                    default=NULL,
		    help="receiver celltype name,like : CD8 T,CD4 T,...")

parser$add_argument("--sender",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="sender celltype name,like : CD8 T,CD4 T,...")


args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

get_lfc_celltype=function (celltype_oi, seurat_obj, condition_colname, condition_oi,
    condition_reference, expression_pct = 0.1)
{
    requireNamespace("Seurat")
    requireNamespace("dplyr")
    seurat_obj_celltype = SetIdent(seurat_obj, value = seurat_obj[["celltype"]])
    seuratObj_sender = subset(seurat_obj_celltype, idents = celltype_oi)
    seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname]])
    DE_table_sender = FindMarkers(object = seuratObj_sender,
        ident.1 = condition_oi, ident.2 = condition_reference,
        min.pct = expression_pct, logfc.threshold = 0.05) %>%
        rownames_to_column("gene")
    DE_table_sender = DE_table_sender %>% as_tibble() %>% select(-p_val) %>%
        select(gene, avg_log2FC)
    colnames(DE_table_sender) = c("gene", celltype_oi)
    return(DE_table_sender)
}


##############################  Paramenters
outDir=args$outdir
condition=args$condition
condition_oi=args$condition_oi
condition_reference=args$condition_reference
receiver=args$receiver
sender_celltypes=args$sender

makedir(outDir)

pct=0.10
lfc_cutoff=0.25
############################## Configure files
message("INFO : Loading Nichenetr Configure ...")
ligand_target_matrix=readRDS("/home/ye/Work/R/SingleCell/Nichenetr/ligand_target_matrix.rds")
weighted_networks=readRDS("/home/ye/Work/R/SingleCell/Nichenetr/weighted_networks.rds")
lr_network=readRDS("/home/ye/Work/R/SingleCell/Nichenetr/lr_network.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
##############################
message("INFO : Loading dataset ...")
seuratObj=readRDS(args$seurat)
seuratObj$celltype=seuratObj$label_fine
Idents(seuratObj)=seuratObj$label_fine

##############################
message("INFO : Perform the NicheNet analysis ...")
message("INFO : Step 1 ...")
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj, pct = pct)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, pct) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()


##############################
message("INFO : Step 2 ...")
seurat_obj_receiver= subset(seuratObj, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[[condition]])

message("INFO : get Different markers ...")
DE_file=file.path(outDir,"DE_table_receiver.rds")
if(!file.exists(DE_file)){
	plan("multiprocess", workers = 16)
        DE_table_receiver = FindMarkers(object = seurat_obj_receiver, 
				ident.1 = condition_oi, ident.2 = condition_reference, min.pct = pct) %>% rownames_to_column("gene")

        message("INFO : Save Result ...")
        saveRDS(DE_table_receiver,file.path(outDir,"DE_table_receiver.rds"))
}else{
	DE_table_receiver=readRDS(DE_file)
}


geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= lfc_cutoff) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

message("INFO : Step 3 ...")
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

message("INFO : Step 4 ...")
message("INFO : Perform NicheNet ligand activity analysis ...")
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
saveRDS(ligand_activities,file.path(outDir,"ligand_activities.rds"))

message("INFO : Step 5 ...")
message("INFO : Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis")
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
saveRDS(vis_ligand_target,file.path(outDir,"vis_ligand_target.rds"))

###################################
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
	theme(axis.text.x = element_text(face = "italic")) + 
	scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
ggsave(file.path(outDir,"ligand_target_network.pdf"),plot=p_ligand_target_network,width=16,height=16)

###################################
message("INFO : Receptors of top-ranked ligands")
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
saveRDS(vis_ligand_receptor_network,file.path(outDir,"vis_ligand_receptor_network.rds"))

#####################################
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

ggsave(file.path(outDir,"ligand_receptor_network.pdf"),plot=p_ligand_receptor_network,width=16,height=16)


#####################################
message("INFO : Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases")
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

saveRDS(vis_ligand_receptor_network_strict,file.path(outDir,"vis_ligand_receptor_network_strict.rds"))

####################################
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
ggsave(file.path(outDir,"ligand_receptor_network_strict.pdf"),plot=p_ligand_receptor_network_strict,width=16,height=16)


###################################
message("INFO : Step 6 ...")
message("INFO : Add log fold change information of ligands from sender cells ..")
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = condition, condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = pct) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]
colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
saveRDS(vis_ligand_lfc,file.path(outDir,"vis_ligand_lfc.rds"))

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
ggsave(file.path(outDir,"ligand_lfc.pdf"),plot=p_ligand_lfc,width=16,height=16)

####################################
message("INFO : Step 7 ...")
message("INFO : Summary visualizations of the NicheNet analysis...")
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
saveRDS(vis_ligand_pearson,file.path(outDir,"vis_ligand_pearson.rds"))

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
ggsave(file.path(outDir,"ligand_pearson.pdf"),plot=p_ligand_pearson,width=16,height=16)


#####################################
# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" 
order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" 
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype %in% sender_celltypes), features = order_ligands_adapted %>% rev(), cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
    ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
pdf(file.path(outDir,"combined_plot.pdf"),width=24,height=24)
print(combined_plot)
dev.off()

message("INFO : Done!")
