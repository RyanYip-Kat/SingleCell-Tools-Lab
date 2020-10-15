library(CATALYST)
library(stringr)
library(argparse)
library(flowCore)
library(Seurat)
library(harmony)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--samples",
                    type="character",
                    default=NULL,
                    help="the path of fcs files")


parser$add_argument("--outdir",
                      type="character",
                      default="output",
                      help="save path")


parser$add_argument("--number",
                    type="integer",
                    default="2000")

parser$add_argument("--config",
                    type="character",
                    default=NULL,
		    help="config file : etc,channel,markers")

parser$add_argument("--batch_correct",
                    dest="batch_correct",
                    action="store_true")

parser$add_argument("--batch_key",
                    type="character",
                    default="sample")

parser$add_argument("--pt_size",
                    type="double",
                    default=0.75)

parser$add_argument("--transformation",
                    action="store_true",
		    default=TRUE)

args<-parser$parse_args()
if(!dir.exists(args$outdir)){
            dir.create(args$outdir,recursive=TRUE)
}

samples_DF=read.csv(args$samples,stringsAsFactors=FALSE,sep=",")
stopifnot("Condition"%in%colnames(samples_DF))
stopifnot("Path"%in%colnames(samples_DF))  # fcs file
stopifnot("Sample"%in%colnames(samples_DF)) # fcs sample name

print(paste0("There are : ",nrow(samples_DF)," sample fcs file"))

print("### Loading config file")
config=read.csv(args$config,stringsAsFactors=FALSE,sep=",")
stopifnot("marker"%in%colnames(config))
stopifnot("channel"%in%colnames(config)) 

channel=config$channel
markers=config$marker
pattern=str_extract(channel,"\\d+")  # extract channel number(unique)
markers_df=data.frame("markers"=markers,"pattern"=pattern)
rownames(markers_df)=markers_df$pattern
pattern=paste(pattern,collapse="|")

names(markers)=channel
print("### Read flowSet")
fset=read.flowSet(samples_DF$Path,which.lines=args$number,column.pattern=pattern,transformation =args$transformation)
column=colnames(fset)
beat=str_extract(column,"\\d+")
marker=markers_df[beat,]$markers
colnames(fset)=marker

#markernames(fset)=markers
colnames(fset)=markers
sampleNames(fset)=unlist(lapply(pData(fset)$name,function(name){return(str_split(name,"\\.")[[1]][1])}))
pData(fset)$Sample=samples_DF$Sample
pData(fset)$Ids=samples_DF$Path
pData(fset)$Condition=samples_DF$Condition

print("### Data preparation and convert into SingleCellExperiment")
fcs_panel=data.frame(fcs_colname=colnames(fset),antigen=colnames(fset),marker_class=sample(c("type","state"),size=length(colnames(fset)),replace=TRUE))
fcs_md=data.frame(file_name=pData(fset)$Ids,sample_id=pData(fset)$Sample,condition=pData(fset)$Condition,patient_id=pData(fset)$Sample)

sce=prepData(fset,fcs_panel,fcs_md)
print(n_cells(sce))

saveRDS(sce,file.path(args$outdir,"sceCytof.rds"))
