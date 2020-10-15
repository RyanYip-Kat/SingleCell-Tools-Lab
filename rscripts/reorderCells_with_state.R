library(monocle)
library(argparse)
print("### configure parameters ###")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--cds",
                    type="character",
                    default="",
                    help="which slot data to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="",
                    help="the dataset  to be used")

parser$add_argument("--state",
                    type="character",
                    default="1",
                    help="the dataset  to be used")

parser$add_argument("--reverse",
                    action="store_true",
                    default=FALSE)
args<-parser$parse_args()

if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

cds<-readRDS(args$cds)
cds <- orderCells(cds,root_state=args$state,reverse=args$reverse)
print("### Saving")
saveRDS(cds,file.path(args$outdir,"cds_reorder.rds"))
