suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(scRepertoire))
suppressMessages(library(Cairo))
suppressMessages(library(Seurat))

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--combined",
		    type="character",
		    default=NULL,
		    help="the combined object from contig list")


parser$add_argument("--groupby",
                    type="character",
                    default="sample",
                    help="the groupby in contig combined object")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the seurat object consistent with contigs")

parser$add_argument("--cloneCall",
		    type="character",
		    default="gene+nt",
		    choices=c("gene+nt","aa","gene","nt"))
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
###############################
height=12
width=12
###############################

print("### Loading combined contig data")
combined=readRDS(args$combined)
seurat=readRDS(args$seurat)

print("### combineExpression with seurat")
seurat <- combineExpression(combined, seurat, cloneCall=args$cloneCall, groupBy = args$groupby)

slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType,
                levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)",
                            "Medium (5 < X <= 20)", "Small (1 < X <= 5)",
                            "Single (0 < X <= 1)", NA))


filename=file.path(args$outdir,"seurat_cloneType.pdf")
DimPlot(seurat, group.by = "cloneType") +
    scale_color_manual(values = colorblind_vector(5), na.value="grey")
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

#seurat <- highlightClonotypes(seurat, cloneCall= "aa", sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF", "NA_CATSATLRVVAEKLFF"))
#DimPlot(seurat, group.by = "highlight")
occupiedscRepertoire(seurat, x.axis = "sample")
filename=file.path(args$outdir,"occupiedscRepertoire_sample.pdf")
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()


