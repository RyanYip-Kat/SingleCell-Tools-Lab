library(argparse)
library(harmony)
library(stringr)
library(flowCore)
library(Rtsne)
library(ggplot2)
#library(sva)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--path",
                    type="character",
                    default="")

parser$add_argument("--batch_correct",
                    dest="batch_correct",
                    action="store_true")

parser$add_argument("--number",
                    type="integer",
                    default="20000")
args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
cols_file<-file.path(args$path,"cols.csv")
if(file.exists(cols_file)){
	cols<-read.csv(file.path(args$path,"cols.csv"),stringsAsFactors=F)
}else{
	cols<-NULL
}

### Run ComBat batch correction from the SVA package
batch.normalise.comBat <- function(counts, batch.groups, max.val=6)
{
  batch.groups = factor(batch.groups) ## drop zero levels
  batch.id = 1:length(unique(batch.groups))
  names(batch.id) = unique(batch.groups)
  batch.ids = batch.id[batch.groups]
  correct.data = ComBat(counts,batch.ids, prior.plots=FALSE, par.prior=T)
  correct.data[correct.data > max.val] = max.val
  as.data.frame(correct.data)
}

#if(!is.null(cols)){
#        column<-as.character(cols$name)
#        column=ifelse(str_detect(column,"\\d+"),str_extract(column,"\\d+"),column)
#        antigen<-as.character(cols$marker)
#        data<-data[,column]
#        colnames(data)<-antigen
#}


print("### Loading data")
files<-list.files(args$path,pattern=".fcs",recursive=TRUE,full.names=TRUE)
print(files)
counts_list<-lapply(files,function(file){
	            p1<-basename(dirname(file))
                    p2<-str_split(basename(file),"\\.")[[1]][1]
		    prefix<-paste(p1,p2,sep="-")
		    ff<-flowCore::read.FCS(file,transformation = FALSE, truncate_max_range = FALSE)
		    data <- flowCore::exprs(ff)
		    Names=colnames(data)
		    col_names=ifelse(str_detect(Names,"\\d+"),str_extract(Names,"\\d+"),Names)
		    colnames(data)=col_names
		    row_names<-paste(prefix,1:nrow(data),sep="_")
		    rownames(data)<-row_names
		    if(!is.null(cols)){
		    	    column<-as.character(cols$name)
                            column=ifelse(str_detect(column,"\\d+"),str_extract(column,"\\d+"),column)
                            antigen<-as.character(cols$marker)
                            data<-data[,column]
                            colnames(data)<-antigen
		    }

		    if(args$number!=0){
			    idx<-sample(1:nrow(data),as.integer(args$number),replace=FALSE)
			    data<-data[idx,]
		    }
		    print(dim(data))
		    return(data)
		    })

print("# Merge data")
data<-do.call(rbind,counts_list)
print(paste0("number of cells : ",nrow(data)))
print(tail(rownames(data)))
print(dim(data))

Names<-rownames(data)
saveRDS(data,file.path(args$outdir,"counts.rds"))
status=unlist(lapply(Names,function(x){return(str_split(x,"-")[[1]][1])}))
sample=unlist(lapply(Names,function(name){return(str_split(name,"_")[[1]][1])}))

metadata<-data.frame("cell_id"=Names,"sample"=sample,"status"=status)
rownames(metadata)=Names
print("# Write data into csv")
write.table(data,file.path(args$outdir,"counts.csv"),sep=",",quote=FALSE)
write.table(metadata,file.path(args$outdir,"cell_meta.csv"),sep=",",quote=FALSE)


print(table(sample))
print(table(status))
if(args$batch_correct){
        metadata<-data.frame("cell_id"=Names,"sample"=sample,"status"=status)
        harmony_embeddings<-HarmonyMatrix(as.matrix(data), metadata, c("sample")) #v3
        rownames(harmony_embeddings)<-Names
        colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")
        saveRDS(harmony_embeddings,file.path(args$outdir,"harmony.rds"))
}


