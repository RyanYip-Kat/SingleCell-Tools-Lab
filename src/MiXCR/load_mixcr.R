suppressMessages(library(argparse))
suppressMessages(library(stringr))
suppressMessages(library(immunarch))

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--path",
		    type="character",
		    default=NULL,
		    help="the mixcr result path")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Load MiXCR data with repLoad")
immdata_mixcr <- repLoad(args$path)
saveRDS(immdata_mixcr,file.path(args$outdir,"immdata_mixcr.rds"))

