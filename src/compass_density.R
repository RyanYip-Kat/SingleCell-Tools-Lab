library(compassR)
library(ggrepel)
library(ggplot2)
library(stringr)
library(argparse)
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

gene_expression=as.data.frame(compass_data$gene_expression_statistics)
rownames(gene_expression)=gene_expression$cell_id
cell_meta=as.data.frame(compass_data$cell_metadata)

rownames(cell_meta)=cell_meta$cell_id
gene_expression=gene_expression[cell_meta$cell_id,]
#gene_expression$status=cell_meta$label_fine
gene_expression$status=cell_meta$status
gene_expression$idents=cell_meta$sample

ggplot(gene_expression,aes(x=metabolic_expression, fill=status)) +
	geom_density(alpha=.3)+theme_bw()+
	theme(axis.text=element_blank())


ggsave(file.path(args$outdir,"metabolic_expression.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()

ggplot(gene_expression,aes(x=metabolic_activity, fill=status)) +
        geom_density(alpha=.3)+theme_bw()+
        theme(axis.text=element_blank())


ggsave(file.path(args$outdir,"metabolic_activity.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()


##################################
ggplot(gene_expression,aes(x=metabolic_expression, fill=idents)) +
        geom_density(alpha=.3)+theme_bw()+
        theme(axis.text=element_blank())


ggsave(file.path(args$outdir,"metabolic_expression_idents.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()

ggplot(gene_expression,aes(x=metabolic_activity, fill=idents)) +
        geom_density(alpha=.3)+theme_bw()+
        theme(axis.text=element_blank())


ggsave(file.path(args$outdir,"metabolic_activity_idents.pdf"),width = 16,height = 16,device =cairo_pdf)
dev.off()

