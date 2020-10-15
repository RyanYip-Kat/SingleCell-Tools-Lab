library(SingleCellSignalR)
library(Seurat)

# Define your working directory
setwd("~/example/")

# Pre-processing using Seurat (https://satijalab.org/seurat/)
pbmc.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc1k")

# Data filtering and normalization
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 40)
pbmc <- NormalizeData(pbmc,scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Data clustering
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.1)

# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
data = data.frame(pbmc[["RNA"]]@data)

# Ligand/Receptor analysis using SingleCellSignalR
signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)

# Visualization
visualize(signal)
intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 3",signal = signal)

