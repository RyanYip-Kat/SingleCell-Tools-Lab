suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(scRepertoire))
suppressMessages(library(Cairo))

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--combined",
		    type="character",
		    default=NULL,
		    help="the combined object from contig list")


parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="subset with column for visual")

parser$add_argument("--subset",
		    nargs="+",
		    type="character",
		    default=NULL,
		    help="the subset in column factors,(eg :A,B,C,D)")

parser$add_argument("--sample",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="plot by group:sample,ID,batch")

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
if(!is.null(args$column) & !is.null(args$subset)){
	variables=args$subset
	combined=subsetContig(combined, name =args$column, variables =variables)
}

print("### Visualizing")
#################################
print(paste0("#### compareClonotypes plot with ",args$cloneCall))
compareClonotypes(combined, numbers = 10, cloneCall=args$cloneCall,samples=args$sample,graph = "alluvial")
filename=file.path(args$outdir,paste0(args$cloneCall,"_alluvial_compareClonotypes.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

compareClonotypes(combined, numbers = 10, cloneCall=args$cloneCall,samples=args$sample, graph = "area")
filename=file.path(args$outdir,paste0(args$cloneCall,"_area_compareClonotypes.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

print(paste0("#### clonalHomeostasis  plot with ",args$cloneCall))
clonalHomeostasis(combined, cloneCall =args$cloneCall)
filename=file.path(args$outdir,paste0(args$cloneCall,"_clonalHomeostasis.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

output=clonalHomeostasis(combined, cloneCall =args$cloneCall,exportTable=TRUE)
filename=file.path(args$outdir,paste0(args$cloneCall,"_clonalHomeostasis.csv"))
write.table(output,filename,sep=",",quote=F,row.names=F)

print(paste0("#### clonalProportion plot with ",args$cloneCall))
clonalProportion(combined, cloneCall =args$cloneCall)
filename=file.path(args$outdir,paste0(args$cloneCall,"_clonalProportion.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

output=clonalProportion(combined, cloneCall =args$cloneCall,exportTable=TRUE)
filename=file.path(args$outdir,paste0(args$cloneCall,"_clonalProportion.csv"))
write.table(output,filename,sep=",",quote=F,row.names=F)

##################################
print("### Overlap Analysis")
print(paste0("#### clonalOverlap plot with ",args$cloneCall))
clonalOverlap(combined, cloneCall =args$cloneCall, method = "morisita")
filename=file.path(args$outdir,paste0(args$cloneCall,"_morisita_clonalOverlap.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

clonalOverlap(combined, cloneCall =args$cloneCall, method = "overlap")
filename=file.path(args$outdir,paste0(args$cloneCall,"_overlap_clonalOverlap.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()


print(paste0("#### clonesizeDistribution plot with ",args$cloneCall))
clonesizeDistribution(combined, cloneCall =args$cloneCall, method="ward.D2")
filename=file.path(args$outdir,paste0(args$cloneCall,"_clonesizeDistribution.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

#########################################
print("### Diversity Analysis")
print(paste0("#### clonalDiversity plot with ",args$cloneCall))
clonalDiversity(combined, cloneCall = args$cloneCall, group = "samples")
filename=file.path(args$outdir,paste0(args$cloneCall,"_sample_clonalDiversity.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

clonalDiversity(combined, cloneCall = args$cloneCall, group = "batch")
filename=file.path(args$outdir,paste0(args$cloneCall,"_batch_clonalDiversity.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

clonalDiversity(combined, cloneCall = args$cloneCall, group = "ID")
filename=file.path(args$outdir,paste0(args$cloneCall,"_ID_clonalDiversity.pdf"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()


