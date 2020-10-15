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

fset=list()
for(i in 1:length(samples_DF$Path)){
	fr=read.FCS(samples_DF$Path[i],which.lines=args$number,column.pattern=pattern,transformation =args$transformation)
	column=colnames(fr)
	beat=str_extract(column,"\\d+")
	marker=markers_df[beat,]$markers
	colnames(fr)=marker
	fset[[samples_DF$Sample[i]]]=fr
}



print("### Rename fset exprs")
cofactor=5
exprs_list=list()
for(i in 1:nrow(samples_DF)){
	sample=samples_DF$Sample[i]
	condition=samples_DF$Condition[i]
	rowName=paste0(condition,"_",paste(sample,1:nrow(exprs(fset[[sample]])),sep="_"))
	rownames(exprs(fset[[sample]]))=rowName
	expr=exprs(fset[[sample]])
	if(!args$transformation){
		expr <- asinh(expr / cofactor) # transform data
	}
	exprs(fset[[sample]])=expr
	exprs_list[[i]]=expr
}
print("### Save flowSet")
saveRDS(fset,file.path(args$outdir,"flowSet.rds"))


print("### Concat dataset")
exprs_matrix=do.call(rbind,exprs_list)
print(paste0("There are : ",nrow(exprs_matrix)," cells"))

Names=rownames(exprs_matrix)
status=unlist(lapply(Names,function(x){return(str_split(x,"_")[[1]][1])}))
samples=unlist(lapply(Names,function(name){return(str_split(name,"_")[[1]][2])}))
metadata=data.frame("cell_id"=Names,"sample"=samples,"status"=status,"condition"=status)
print(table(metadata$sample))
print(table(metadata$condition))

use_rep="pca"
ndims=20
print("### Create Seurat object")
seurat<-CreateSeuratObject(counts=t(exprs_matrix),
                       project ="cytof",
                       names.delim="_",
                       min.cells=0,
                       min.features=1)

print(dim(seurat))
rownames(metadata)=Names
seurat<-AddMetaData(seurat,metadata)

if(args$transformation){
	print("### Normalize Data")
        seurat=NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)
}

if(args$batch_correct){
        print("### Run Harmony")
        harmony_embeddings=HarmonyMatrix(exprs_matrix, metadata, args$batch_key) #v3
        rownames(harmony_embeddings)=Names
        colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")
        use_rep="harmony"
        ndims=20
        #saveRDS(harmony_embeddings,file.path(args$outdir,"harmony.rds"))
	print("### Add harmony")
	mat=as.matrix(harmony_embeddings)
        colnames(mat)<-paste("Harmony_",1:ncol(mat),sep = "")
        seurat[["harmony"]]<-CreateDimReducObject(embeddings =mat,
						  key = "Harmony_",
                                                  assay = DefaultAssay(seurat))
}

if(!args$batch_correct){
	print("### Run PCA")
        seurat<-ScaleData(seurat,features=rownames(seurat))
        seurat<-RunPCA(seurat,npcs=30,features=rownames(seurat))
}
print("### Run UMAP")
seurat <- RunUMAP(object = seurat, reduction =use_rep, dims = 1:ndims)
print("### Run TSNE")
seurat <- RunTSNE(object = seurat, reduction = use_rep, dims = 1:ndims)
print("Run FindNeighbors")
seurat <- FindNeighbors(object = seurat,reduction =use_rep, dims = 1:ndims)

print("Run Find Clusters")
seurat <- FindClusters(object = seurat,resolution=0.8, verbose =TRUE)


saveRDS(seurat,file.path(args$outdir,"seurat.rds"))

jpeg(file.path(args$outdir,"cluster1.jpeg"),width=1024,height=1024)
p<-DimPlot(seurat,reduction="tsne",label=TRUE,label.size=7.5)
print(p)
dev.off()

jpeg(file.path(args$outdir,"cluster2.jpeg"),width=1024,height=1024)
p<-DimPlot(seurat,reduction="umap",label=TRUE,label.size=7.5)
print(p)
dev.off()

print("### Antigen plot")
library(ggplot2)
library(Cairo)
plot_jpeg=file.path(args$outdir,"antigen_jpeg")
plot_pdf=file.path(args$outdir,"antigen_pdf")
if(!dir.exists(plot_jpeg)){
            dir.create(plot_jpeg,recursive=TRUE)
}

if(!dir.exists(plot_pdf)){
            dir.create(plot_pdf,recursive=TRUE)
}


my_theme<-theme(axis.title.x = element_text(size=25),
                  axis.text.x = element_text(size=18),
                  axis.text.y = element_text(size=18),
                  axis.title.y = element_text(size=25),
                  plot.title=element_text(size=25,face="bold"),
                  legend.text = element_text(size=25),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())

antigens<-rownames(seurat)
for(antigen in antigens){
        jpeg(file.path(plot_jpeg,paste0(antigen,"_tsne_condition_seurat.jpeg")),width=1024,height=1024)
        p=FeaturePlot(seurat,reduction="tsne",pt.size=args$pt_size,features=antigen,split.by="condition")+
                my_theme
        print(p)
        dev.off()

        jpeg(file.path(plot_jpeg,paste0(antigen,"_tsne_seurat.jpeg")),width=1024,height=1024)
        p=FeaturePlot(seurat,reduction="tsne",pt.size=args$pt_size,features=antigen)+my_theme
        print(p)
        dev.off()


}

cols=c("#1A237E","#F8FCB7","#DD2C00")
for(antigen in antigens){
        p=FeaturePlot(seurat,reduction="tsne",pt.size=args$pt_size,features=antigen,split.by="condition",cols=cols)+
                my_theme
        print(p)
        ggsave(file.path(plot_pdf,paste0(antigen,"_tsne_condition_seurat.pdf")),width=12,height=12,device =cairo_pdf)
        dev.off()

        p=FeaturePlot(seurat,reduction="tsne",pt.size=args$pt_size,features=antigen,cols=cols)+my_theme
        print(p)
        ggsave(file.path(plot_pdf,paste0(antigen,"_tsne_seurat.pdf")),width=12,height=12,device =cairo_pdf)
        dev.off()

        p=FeaturePlot(seurat,reduction="umap",pt.size=args$pt_size,features=antigen,split.by="condition",cols=cols)+
                my_theme
        print(p)
        ggsave(file.path(plot_pdf,paste0(antigen,"_umap_condition_seurat.pdf")),width=12,height=12,device =cairo_pdf)
        dev.off()

        p=FeaturePlot(seurat,reduction="umap",pt.size=args$pt_size,features=antigen,cols=cols)+my_theme
        print(p)
        ggsave(file.path(plot_pdf,paste0(antigen,"_umap_seurat.pdf")),width=12,height=12,device =cairo_pdf)
        dev.off()


}
