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

parser$add_argument("--groupby",
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
print(paste0("#### quantContig plot with ",args$cloneCall))
quantContig(combined, cloneCall=args$cloneCall, scale = T,group=args$groupby)
filename=file.path(args$outdir,paste0(args$cloneCall,"_quantContig.pdf"))
if(!is.null(args$groupby)){
        filename=file.path(args$outdir,paste0(args$cloneCall,"_",args$groupby,"_quantContig.pdf"))
}

ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()


quantContig_output <- quantContig(combined, cloneCall=args$cloneCall, scale = F, exportTable = T,group=args$groupby)
filename=file.path(args$outdir,paste0(args$cloneCall,"_quantContig.csv"))
if(!is.null(args$groupby)){
        filename=file.path(args$outdir,paste0(args$cloneCall,"_",args$groupby,"_quantContig.csv"))
}

write.table(quantContig_output,filename,sep=",",quote=F,row.names=F)

#################################
print(paste0("#### abundanceContig plot with ",args$cloneCall))
abundanceContig(combined, cloneCall = args$cloneCall, scale = T,group=args$groupby)
filename=file.path(args$outdir,paste0(args$cloneCall,"_abundanceContig.pdf"))
if(!is.null(args$groupby)){
	filename=file.path(args$outdir,paste0(args$cloneCall,"_",args$groupby,"_abundanceContig.pdf"))
}

ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()


abundanceContig_output=abundanceContig(combined, cloneCall =args$cloneCall, exportTable = T,group=args$groupby)
filename=file.path(args$outdir,paste0(args$cloneCall,"_abundanceContig.csv"))
if(!is.null(args$groupby)){
        filename=file.path(args$outdir,paste0(args$cloneCall,"_",args$groupby,"_abundanceContig.csv"))
}

write.table(abundanceContig_output,filename,sep=",",quote=F,row.names=F)

###############################
print(paste0("#### compareClonotypes plot with ",args$cloneCall))
compareClonotypes(combined, numbers = 10, cloneCall=args$cloneCall, graph = "alluvial")
filename=file.path(args$outdir,paste0(args$cloneCall,"_alluvial_compareClonotypes.csv"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

compareClonotypes(combined, numbers = 10, cloneCall=args$cloneCall, graph = "area")
filename=file.path(args$outdir,paste0(args$cloneCall,"_area_compareClonotypes.csv"))
ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
dev.off()

