library(stringr)
library(argparse)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--data",
                    type="character",
                    default="",
                    help="")
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

split_status_DE<-function(DATA,outdir){
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = TRUE)
  }
  stopifnot("cluster"%in%colnames(DATA))
  clusters=unique(DATA$cluster)
  for(c in clusters){
    df=subset(DATA,cluster==c )
    filename=paste0(c,"_markers.csv")
    write.table(df,file.path(outdir,filename),sep=",",quote = FALSE,row.names = F)

  }
}

DATA=read.csv(args$data,sep=",")
split_status_DE(DATA,args$outdir)
