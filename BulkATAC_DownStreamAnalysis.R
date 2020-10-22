library(argparse)
library(stringr)
library(Matrix)
library(edgeR)
library(pheatmap)
library(matrixStats)
library(SummarizedExperiment)
library(BiocParallel)
library(chromVAR)
library(motifmatchr)
#library(BSgenome.Hsapiens.UCSC.hg38)
#library(BSgenome.Mmusculus.UCSC.mm10)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--featureCounts",
                    type="character",
                    default="the path of the featureCounts.txt")


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


parser$add_argument("--excludeChr",
		    nargs="+",
                    type="character",
                    default=NULL)


#parser$add_argument("--binarize",
#                    action="store_true")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

jobs=args$jobs
cat(sprintf("### Use %d cores to run the scripts\n",jobs))
register(MulticoreParam(jobs)) 


cat(sprintf("### Loading the featureCounts from : %s\n",args$featureCounts))
counts=read.table(args$featureCounts,skip=1,header=TRUE)

if(!is.null(args$excludeChr)){
	if(length(args$excludeChr)>0){
		exclude=paste(args$excludeChr,collapse="|")
	}
	cat(sprintf("### Exclude Chr Pattern : %s\n",exclude))
	chroms_idx=str_detect(counts$Chr,exclude)
	counts=counts[!chroms_idx,]
}

print("### Create SummarizedExperiment Object")
counts_matrix=counts[,7:ncol(counts)]
meta=counts[,c(1,6)]
rowranges=counts[,c(2,3,4,5)]
colnames(rowranges)=c("chr","start","end","strand")
gs=makeGRangesFromDataFrame(rowranges)
mcols(gs)=meta

rows=paste(rowranges$chr,rowranges$start,rowranges$end,sep="-")
rownames(counts_matrix)=rows
cols=unlist(lapply(colnames(counts_matrix),function(x)return(str_split(x,"\\.")[[1]][1])))

se=SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts_matrix)),
                          rowRanges=gs,
                          colData=DataFrame("Status"=cols))

rownames(se)=rows

############################
if(str_to_lower(args$species)%in%c("human","hs")){
	library(BSgenome.Hsapiens.UCSC.hg38)
	genome=BSgenome.Hsapiens.UCSC.hg38
	species="Homo sapiens"
}else if(str_to_lower(args$species)%in%c("mouse","mm")){
	library(BSgenome.Mmusculus.UCSC.mm10)
	genome=BSgenome.Mmusculus.UCSC.mm10
	species="Mus musculus"
}else{
	stop("### Invalid species provived!!!")
}


ser<- addGCBias(se,genome = genome)
#ser <- filterSamples(ser, min_depth = 1500,
#                                  min_in_peaks = 0.15)
#ser <- filterPeaks(ser)

print("### Match Motif")

motifs <- getJasparMotifs(species=species)
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

print("### Plots")
variability <- computeVariability(dev)
pdf(file.path(args$outdir,"Plots.pdf"),width=12,height=12)
p1=plotVariability(variability, n=5,use_plotly = FALSE)
print(p1)

sample_cor <- getSampleCorrelation(dev)

pheatmap(as.dist(sample_cor), 
         annotation_row = as.data.frame(colData(dev)), 
         clustering_distance_rows = as.dist(1-sample_cor), 
	 clustering_distance_cols = as.dist(1-sample_cor))
dev.off()


print("### differentialVariability")
diff_var <- differentialVariability(dev, "Status")
saveRDS(diff_var,file.path(args$outdir,"diff_var.rds"))
print(head(diff_var))




