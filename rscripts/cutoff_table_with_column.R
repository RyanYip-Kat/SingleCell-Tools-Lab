library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--table",
                    type="character",
                    default="")

parser$add_argument("--column",
                    type="character",
                    default="pvals")

parser$add_argument("--cutoff",
                    type="double",
                    default="0.05")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

df=read.csv(args$table,stringsAsFactors=F)
if(!args$column%in%colnames(df)){
	stop("Invalid columns input!!!")
}

if(class(df[[args$column]])=="numeric"){
	m<-df[df[[args$column]]<args$cutoff,]
	write.table(m,file.path(args$outdir,"filtered_table.csv"),sep=",",quote=F,row.names=F)
}

