library(argparse)
library(stringr)
library(Matrix)
library(edgeR)
library(matrixStats)
library(SummarizedExperiment)
library(BiocParallel)
library(chromVAR)
library(motifmatchr)
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Mmusculus.UCSC.mm10)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="BulkATAC_result")


parser$add_argument("--jobs",
                    type="integer",
                    default=16,
                    help="number of threads for command")


parser$add_argument("--species",
                    type="character",
                    default="human",
		    choices=c("human","hs","mouse","mm"))


parser$add_argument("--bamfiles",
		    nargs="+",
                    type="character",
                    default=NULL,
		    help="bam files from Bulk ATAC Alignment pipeline")


parser$add_argument("--peakfile",
                    type="character",
                    default=NULL,
		    help="the peak file(of *.bed)")

parser$add_argument("--names",
		    nargs="+",
                    type="character",
                    default=NULL,
		    help="the name of each bam file")


#parser$add_argument("--binarize",
#                    action="store_true")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

jobs=args$jobs
cat(sprintf("### Use %d cores to run the scripts\n",jobs))
register(MulticoreParam(jobs)) 

peakfile=args$peakfile
cat(sprintf("### Loading the peakfile from : %s\n",peakfile))
peaks <- getPeaks(peakfile, sort_peaks = TRUE)

bamfiles=args$bamfiles
cat(sprintf("### Loading  %d bamfiles\n",length(bamfiles))
cat(sprintf("--- : %s\n",bamfiles))

if(!is.null(args$names)){
	stopifnot(length(bamfiles)==length(args$names))
}

se = getCounts(bamfiles, peaks,
	       paired = FALSE, 
               by_rg = TRUE, 
               format = "bam", 
               colData = DataFrame("Status" =args$names))


############################
if(str_to_lower(args$species)%in%c("human","hs")){
	library(BSgenome.Hsapiens.UCSC.hg38)
	genome=BSgenome.Hsapiens.UCSC.hg38
}else if(str_to_lower(args$species)%in%c("mouse","mm")){
	library(BSgenome.Mmusculus.UCSC.mm10)
	genome=BSgenome.Mmusculus.UCSC.mm10
}else{
	stop("### Invalid species provived!!!")
}


ser<- addGCBias(se,genome = genome)
#ser <- filterSamples(ser, min_depth = 1500,
#                                  min_in_peaks = 0.15)
#ser <- filterPeaks(ser)

print("### Match Motif")
motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, ser,
                         genome =genome)

# computing deviations
print("### Compute Deviations")
dev <- computeDeviations(object = ser,annotations = motif_ix)

print("### Save Objects")
dev_file=file.path(args$outdir,"SummarizedExperiment-Deviations.rds")
se_file=file.path(args$outdir,"SummarizedExperiment-SE.rds")
cat(sprintf("### Save Se Object into : %s\n",se_file))
cat(sprintf("### Save Deviations Object into : %s\n",dev_file))
saveRDS(se,se_file)
saveRDS(dev,dev_file)
print("### Done!")

