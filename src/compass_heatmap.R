library(compassR)
library(ggrepel)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(stringr)
library(argparse)
library(dplyr)
library(pheatmap)
library(Cairo)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default=NULL,
                    help="directory should include files named :cell_metadata.csv, reactions.tsv, and linear_gene_expression_matrix.tsv")


parser$add_argument("--symbol",
                    type="character",
                    default="HGNC.symbol")

parser$add_argument("--outdir",
                    type="character",
                    default="compass_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


compass_data=readRDS(args$data)
wilcoxon_results=read.csv(file.path(dirname(args$data),"wilcoxon_results.csv"),stringsAsFactors =F)

facets <- c( "Glycolysis", "TCA cycle", "Fatty acid oxidation", "Amino acid metabolism")

compass_scores_by_cell_type <-
    wilcoxon_results %>%
    left_join(
        select(compass_data$reaction_partitions, "reaction_id", "reaction_no_direction"),
        by = "reaction_id"
    ) %>%
    left_join(
        compass_data$reaction_metadata,
        by = "reaction_no_direction"
    ) %>%
    # Keep only "confident reactions", as defined in our paper.
    filter(!is.na(EC_number)) %>%
    filter(confidence == "0" | confidence == "4") %>%
    # Exclude non-mitochondrially localized reactions from TCA.
    mutate(subsystem = case_when(
        reaction_id == "SPMDOX_pos" ~ "Arginine and Proline Metabolism",
        subsystem == "Citric acid cycle" & !grepl("[m]", formula, fixed = TRUE) ~ "Other",
        TRUE ~ subsystem
    )) %>%
    # Assign reactions to the appropriate subsystem.
    mutate(
        subsystem_priority = factor(subsystem) %>%
        fct_recode(
            "Glycolysis" = "Glycolysis/gluconeogenesis",
            "TCA cycle" = "Citric acid cycle"
        ) %>%
        fct_collapse("Amino acid metabolism" = c(
            "Alanine and aspartate metabolism",
            "Arginine and Proline Metabolism",
            "beta-Alanine metabolism",
            "Cysteine Metabolism",
            "D-alanine metabolism",
            "Folate metabolism",
            "Glutamate metabolism",
            "Glycine, serine, alanine and threonine metabolism",
            "Histidine metabolism",
            "Lysine metabolism",
            "Methionine and cysteine metabolism",
            "Taurine and hypotaurine metabolism",
            "Tryptophan metabolism",
            "Tyrosine metabolism",
            "Urea cycle",
            "Valine, leucine, and isoleucine metabolism"
        )) %>%
        fct_other(keep = facets) %>%
        fct_relevel(facets)
    ) %>%
    # Keep only the subsystems for which we want to plot a facet.
    filter(subsystem_priority != "Other") %>%
    # Lower-bound the adjusted p-value.
    mutate(adjusted_p_value = if_else(
        subsystem_priority == "Amino acid metabolism" & adjusted_p_value <= 1e-12,
        1e-12,
        adjusted_p_value
    )) %>%
    # Assign descriptive labels to various reactions.
    mutate(label = case_when(
        reaction_id == "PGM_neg" ~ "phosphoglycerate mutase (PGAM)",
        reaction_id == "LDH_L_neg" ~ "lactate dehydrogenase",
        reaction_id == "PDHm_pos" ~ "pyruvate dehydrogenase (PDH)",
        reaction_id == "TPI_neg" ~ "triosephosphate isomerase (DHAP forming)",
        reaction_id == "FACOAL1821_neg" ~ "long-chain fatty-acid-CoA ligase",
        reaction_id == "r1257_pos" ~ "long-chain fatty-acid-CoA ligase",
        reaction_id == "FACOAL1831_neg" ~ "long-chain fatty-acid-CoA ligase",
        reaction_id == "CSNATr_neg" ~ "carnitine O-acetyltransferase",
        reaction_id == "C160CPT1_pos" ~ "carnitine O-palmitoyltransferase",
        reaction_id == "ACONTm_pos" ~ "aconitate hydratase",
        reaction_id == "SUCOASm_pos" ~ "succinate-CoA ligase",
        reaction_id == "AKGDm_pos" ~ "alpha-ketoglutarate dehydrogenase",
        reaction_id == "SUCD1m_pos" ~ "succinate dehydrogenase",
        reaction_id == "ICDHyrm_pos" ~ "isocitrate dehydrogenase",
        reaction_id == "CK_pos" ~ "creatine\nkinase",
        reaction_id == "PGCD_pos" ~ "phosphoglycerate dehydrogenase",
        reaction_id == "ARGSS_pos" ~ "arginosuccinate synthase",
        reaction_id == "r0281_neg" ~ "putrescine diamine oxidase",
        reaction_id == "SPMDOX_pos" ~ "spermidine dehydrogenase (spermidine -> GABA)",
        reaction_id == "ARGDCm_pos" ~ "arginine decarboxylase",
        reaction_id == "AGMTm_pos" ~ "agmatinase",
        reaction_id == "GHMT2r_pos" ~ "serine hydroxymethyltransferase",
        reaction_id == "AHC_pos" ~ "adenosylhomocysteinase",
        reaction_id == "METAT_pos" ~ "methionine adenosyltransferase",
        reaction_id == "METS_pos" ~ "methionine\nsynthase",
        reaction_id == "ARGN_pos" ~ "arginase",
        TRUE ~ ""
    ))

#write.table(compass_scores_by_cell_type,file.path(args$outdir,"compass_scores_by_cell.csv"),sep=",",quote=F,row.names=F)
mat=compass_data$reaction_consistencies[compass_scores_by_cell_type$reaction_id,]
cell_meta=as.data.frame(compass_data$cell_metadata)

cell_meta=cell_meta[order(cell_meta$status),]
mat=mat[,cell_meta$cell_id]

ac=data.frame(Status=cell_meta$status,Sample=cell_meta$idents)
rownames(ac)=colnames(mat)

pdf("pheatmap_status.pdf",width=8,height=32)
print(pheatmap(as.matrix(mat),
	       show_colnames =F,
	       show_rownames = T,
	       cluster_rows = F,
	       cluster_cols=F,
	       fontsize=5,
	       scale="column",
	       annotation_col=ac))
dev.off()
