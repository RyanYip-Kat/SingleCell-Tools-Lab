import argparse
import os
#import stream as st
import scanpy as sc
import numpy as np
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-s","--seurat",type=str,default="./seurat_export",help="seurat export path:should include cell_label.txt,cell.txt,X_tsne,X_umap.txt")
parser.add_argument("-d","--data",type=str,default=None,help="path of cellranger count or aggr")

args=parser.parse_args()


sc._settings.ScanpyConfig.n_jobs=4
sc._settings.ScanpyConfig.autoshow=False
sc.set_figure_params(figsize=[36,24])
sc.settings.verbosity = 3

matrix_path=os.path.join(args.data,"outs/filtered_feature_bc_matrix")
graphclust_path=os.path.join(args.data,"outs/analysis/clustering/graphclust/clusters.csv")

print("-----------")
print("Loading counts matrix from : {}".format(matrix_path))
print("Loading original cluster  from : {}".format(graphclust_path))

adata=sc.read_10x_mtx(matrix_path,cache=True)
print("origin adata shape is [{},{}]".format(adata.shape[0],adata.shape[1]))

print("*** add meta data ***")
orig_cluster=pd.read_csv(graphclust_path)
adata.obs["orig_cluster"]=orig_cluster.Cluster.to_list()


cells=adata.obs_names.to_list()
idents=[cell.split("-")[1] for cell in cells]
adata.obs["idents"]=idents

cell_subset=pd.read_table(os.path.join(args.seurat,"cell.txt"),sep="\t",header=None)
cell_subset=cell_subset.iloc[:,0].values
cell_label_subset=pd.read_table(os.path.join(args.seurat,"cell_label.txt"),sep="\t",header=None)
cell_label_subset=cell_label_subset.iloc[:,0].values

tsne_mat=pd.read_table(os.path.join(args.seurat,"X_tsne.txt"),sep="\t")
tsne_mat=np.array(tsne_mat)

umap_mat=pd.read_table(os.path.join(args.seurat,"X_umap.txt"),sep="\t")
umap_mat=np.array(umap_mat)

def subset(adata,cells=None,features=None):
    Cells=adata.obs_names
    Features=adata.var_names
    if cells is not None:
       cell_idx=[True if cell in cells else False for cell in Cells]
       adata=adata[cell_idx,:]
    if features is not None:
       feature_idx=[True if feature in features else False for feature in Features]
       adata=adata[:,feature_idx]
    
    return adata
    

adata_subset=subset(adata,cell_subset)
adata_subset.obs["cell_label"]=cell_label_subset

adata_subset.obsm["X_tsne"]=tsne_mat
adata_subset.obsm["X_umap"]=umap_mat

adata_subset.write(os.path.join(args.outdir,"adata_seurat.h5ad"), compression='gzip')

       
        
    



