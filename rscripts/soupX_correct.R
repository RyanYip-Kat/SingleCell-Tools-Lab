library(stringr)
library(argparse)
library(SoupX)
library(DropletUtils)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--crg",
                    type="character",
                    default=NULL,
                    help="the output from cellranger counts or aggr" )
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="output path")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
sc = load10X(args$crg)

print("### Correct")
sc = autoEstCont(sc)
out = adjustCounts(sc)

print("### Write into mtx matrix")
write10xCounts(x =out, path=file.path(args$outdir,"strainedCounts"),overwrite=TRUE)
