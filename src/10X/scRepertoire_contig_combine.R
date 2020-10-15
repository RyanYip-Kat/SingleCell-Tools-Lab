suppressMessages(library(argparse))
suppressMessages(library(stringr))
suppressMessages(library(scRepertoire))

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--meta",
		    type="character",
		    default=NULL,
		    help="meta csv file : Path,Sample,Batch(ID)")

parser$add_argument("--class",
		    type="character",
		    default="BCR",
		    choices=c("BCR","TCR"))
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

meta=read.csv(args$meta,stringsAsFactors=FALSE)
meta$ID=1:nrow(meta)
stopifnot("Path"%in%colnames(meta))
stopifnot("Sample"%in%colnames(meta))
stopifnot("Batch"%in%colnames(meta))


contig_list=lapply(1:length(meta$Path),function(i){
			   file=meta$Path[i]
			   contig=read.csv(file,stringsAsFactors=FALSE)
			   barcode=unlist(lapply(contig$barcode,function(b)return(str_split(b,"-")[[1]][1])))
			   barcode=paste(barcode,i,sep="-")
			   contig$barcode=barcode

			   contig_id=unlist(lapply(contig$contig_id,function(b)return(str_split(b,"_")[[1]][2])))
			   contig_id=paste(barcode,contig_id,sep="_")
			   contig$contig_id=contig_id
			   return(contig)
		    })
if(args$class=="TCR"){
	combined <- combineTCR(contig_list, samples = meta$Sample, ID =meta$Batch , cells ="T-AB")
}else{
	combined=combineBCR(contig_list,samples = meta$Sample, ID =meta$Batch)
}

combined<- addVariable(combined, name = "batch", variables = meta$Batch)
saveRDS(combined,file.path(args$outdir,"scRepertoire_combine.rds"))
