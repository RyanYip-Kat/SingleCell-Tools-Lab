import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys


parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--seurat_dir",type=str,default=None)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("# Loading seurat result from : {}".format(args.seurat_dir))

#counts_dir=os.path.join(args.seurat_dir,"counts.csv")
counts_dir=os.path.join(args.seurat_dir,"matrix")
meta_dir=os.path.join(args.seurat_dir,"cell_meta.csv")
projections=[file for file in  os.listdir(args.seurat_dir) if "projection" in file]

print("# Create adata with scanpy")
#adata=sc.read_csv(counts_dir)
adata=sc.read_10x_mtx(counts_dir,cache=True)
adata.raw=adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2000)

meta=pd.read_csv(meta_dir)
meta=meta.loc[adata.obs_names,:]
adata.obs=meta

print("# Add Embedding")
for proj in projections:
    proj_dir=os.path.join(args.seurat_dir,proj)
    data=pd.read_csv(proj_dir,index_col=0)
    data=np.array(data)
    assert data.shape[0]==adata.shape[0]
    try:
        method=proj.split("_")[0]
        emb="X_"+method
        print("## Add method : {} into obsm : {}".format(method,emb))
        adata.obsm[emb]=data
    except:
        print("Please check your adata and embedding  data!!!")

print("# Save data")
adata.write(os.path.join(args.outdir,"seurat_adata.h5ad"))
