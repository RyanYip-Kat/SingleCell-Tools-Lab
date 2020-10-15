# Reproduces Fig. 2(e) from our paper. (This script works out-of-the-box with the Th17 cell data set included in the package.)

library(compassR)
library(ggrepel)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(stringr)
library(argparse)
library(dplyr)
library(Cairo)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default=NULL,
                    help="directory should include files named :cell_metadata.csv, reactions.tsv, and linear_gene_expression_matrix.tsv")


parser$add_argument("--outdir",
                    type="character",
                    default="compass_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


compass_data=readRDS(args$data)
wilcoxon_results=read.csv(file.path(dirname(args$data),"wilcoxon_results.csv"),stringsAsFactors =F)

cohens_d_by_subsystem <-
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
    # Keep only "interesting subsystems", as defined in our paper.
    filter(!(subsystem == "Miscellaneous" | subsystem == "Unassigned")) %>%
    filter(!(startsWith(subsystem, "Transport") | startsWith(subsystem, "Exchange"))) %>%
    # Keep only subsystems of non-negligible size.
    group_by(subsystem) %>%
    filter(n() > 5) %>%
    ungroup() %>%
    # Order subsystems in a manner that will lend itself to a visually aesthetic plot.
    mutate(
        subsystem_priority = factor(subsystem) %>%
        fct_reorder2(
            cohens_d,
            adjusted_p_value,
            .fun = function(cohens_d, adjusted_p_value) {
                abs(median(cohens_d[adjusted_p_value < 0.1]))
            },
            .desc = FALSE
        )
    )

write.table(cohens_d_by_subsystem,file.path(args$outdir,"cohens_d_by_subsystem.csv"),sep=",",quote=F,row.names=F)
saveRDS(cohens_d_by_subsystem,file.path(args$outdir,"cohens_d_by_subsystem.rds"))
cols=c("associated_genes","subsystem","EC_number","confidence","subsystem_priority")
cohens_d_by_subsystem_df=cohens_d_by_subsystem[,cols]
write.table(cohens_d_by_subsystem_df,file.path(args$outdir,"cohens_d_by_subsystem_2.csv"),sep=",",quote=F,row.names=F)

cohens_d_by_subsystem_df=cohens_d_by_subsystem[,setdiff(colnames(cohens_d_by_subsystem),cols)]
write.table(cohens_d_by_subsystem_df,file.path(args$outdir,"cohens_d_by_subsystem_1.csv"),sep=",",quote=F,row.names=F)

cols=c("cohens_d","p_value","adjusted_p_value","reaction_no_direction","subsystem_priority","reaction_id")
cohens_d_by_subsystem_df=cohens_d_by_subsystem[,cols]
write.table(cohens_d_by_subsystem_df,file.path(args$outdir,"cohens_d_by_subsystem_3.csv"),sep=",",quote=F,row.names=F)
ggplot(
    cohens_d_by_subsystem,
    aes(
        x = subsystem_priority,
        y = cohens_d,
        color = if_else(cohens_d > 0, "up_regulated", "down_regulated"),
        alpha = if_else(adjusted_p_value < 0.1, "significant", "insignificant")
    )
) +
ggtitle("Up- and Down-Regulated Reactions Cross Pathway Boundaries") +
xlab("") + ylab("Cohen's d") +
scale_color_manual(
    values = c(up_regulated = "#ca0020", down_regulated = "#0571b0"),
    guide = FALSE
) +
scale_alpha_manual(
    name = "",
    values = c(significant = 1, insignificant = 0.25),
    labels = c(significant = "BH-adjusted p-value < 0.1", insignificant = "insignificant")
) +
coord_flip() +
geom_point() +
geom_hline(yintercept = 0, linetype = "dashed") +
theme_bw() +
theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(file.path(args$outdir,"compass_up_down_plot.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()
