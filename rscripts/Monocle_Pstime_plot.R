library(ggplot2)
library(monocle)
library(argparse)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--cds",
                    type="character",
                    default="",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--genes",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="")
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

cds<-readRDS(args$cds)
#pdf(file.path(args$outdir,"mono_celltype.pdf"),width=16,height=12)
#plot_cell_trajectory(cds, color_by = "label_fine")
#dev.off()

pdf(file.path(args$outdir,"mono_states.pdf"),width=16,height=12)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf(file.path(args$outdir,"mono_states_facet.pdf"),width=16,height=12)
plot_cell_trajectory(cds, color_by = "State") +
    facet_wrap(~State, nrow = 1)
dev.off()

pdf(file.path(args$outdir,"mono_states_Pseudotime.pdf"),width=16,height=12)
plot_cell_trajectory(cds, color_by = "Pseudotime")
dev.off()

pdf(file.path(args$outdir,"mono_status_facet.pdf"),width=16,height=12)
plot_cell_trajectory(cds, color_by = "Status") +
    facet_wrap(~Status, nrow = 1)
dev.off()

pdf(file.path(args$outdir,"mono_status.pdf"),width=16,height=12)
plot_cell_trajectory(cds, color_by = "Status") 
dev.off()

if(!is.null(args$genes)){
	print(paste0("plot :",length(args$genes)," genes in pseudotime"))
	my_genes <- row.names(subset(fData(cds), gene_short_name %in% args$genes))
        cds_subset <- cds[my_genes,]

	pdf(file.path(args$outdir,"mono_states_ps.pdf"),width=16,height=12)
        print(plot_genes_in_pseudotime(cds_subset, color_by ="State"))
	dev.off()

	pdf(file.path(args$outdir,"mono_status_ps.pdf"),width=16,height=12)
	print(plot_genes_in_pseudotime(cds_subset, color_by ="Status"))
	dev.off()

}

#plot_cell_trajectory(cds, markers="MYH3",color_by="State")

