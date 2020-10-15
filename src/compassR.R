library(compassR)
library(ggrepel)
library(tidyverse)
library(magrittr)
library(stringr)
library(argparse)
library(dplyr)
#source("compass_function.R")

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

print("# Compass setting")
files=list.files(args$data)
stopifnot("cell_metadata.csv"%in%files)
stopifnot("reactions.tsv"%in%files)
stopifnot("linear_gene_expression_matrix.tsv"%in%files)


compass_settings <- CompassSettings$new(
    user_data_directory = args$data,
    cell_id_col_name = "cell_id",
    gene_id_col_name = args$symbol
)

print("# Loading data")
compass_data <- CompassData$new(compass_settings)

print("# Exploring the statistical analysis suite")
compass_analyzer <- CompassAnalyzer$new(compass_settings)

saveRDS(compass_data,file.path(args$outdir,"compass_data.rds"))
saveRDS(compass_analyzer,file.path(args$outdir,"compass_analyzer.rds"))

group_A_cell_ids <-
    compass_data$cell_metadata %>%
    #filter(label_fine == "IL8+CM") %>%
    filter(status=="PDR")%>%
    pull(cell_id)
group_B_cell_ids <-
    compass_data$cell_metadata %>%
    #filter(label_fine == "S100A8+CM") %>%
    filter(status=="HC")%>%
    pull(cell_id)

#wilcoxon_results=conduct_wilcoxon_test_1(compass_data,group.by="status",for_metareactions=FALSE)
wilcoxon_results <- compass_analyzer$conduct_wilcoxon_test(
    compass_data$reaction_consistencies,
    group_A_cell_ids,
    group_B_cell_ids,
    for_metareactions = FALSE
)

print(head(wilcoxon_results))
write.table(wilcoxon_results,file.path(args$outdir,"wilcoxon_results.csv"),sep=",",quote=F,row.names=F)

reaction_metadata=compass_data$reaction_metadata
reaction_metadata$reaction_name=str_replace(reaction_metadata$reaction_name,",","|")
reaction_metadata=subset(reaction_metadata,!is.na(EC_number))
write.table(reaction_metadata,file.path(args$outdir,"reaction_metadata.csv"),sep=",",quote=F,row.names=F)

reaction_partitions=compass_data$reaction_partitions
write.table(reaction_partitions,file.path(args$outdir,"reaction_partitions.csv"),sep=",",quote=F,row.names=F)
