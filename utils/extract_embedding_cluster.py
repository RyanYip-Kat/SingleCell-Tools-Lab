import scanpy as sc
import numpy as np
import pandas as pd
import scipy.io
import argparse
import os
from collections import Counter

import sys
sys.path.append("../")
from core.utils import subset_by_column,subset_by_cell_feature


parser=argparse.ArgumentParser()
parser.add_argument("-s","--subset",nargs="+",type=str,default=None)
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None,help="adata object")
parser.add_argument("-c","--column",type=str,default="leiden")
parser.add_argument("-v","--invert",action="store_true",default=False)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)


print("### Laoding data")
adata=sc.read_h5ad(args.data)

if args.subset is not None:
    print("*** subset")
    print("Subset data by {}".format(args.column))
    adata=subset_by_column(adata,args.subset,args.column,invert=args.invert)

metadata=adata.obs.copy()
metadata["barcode"]=metadata.index
print(metadata.head())

assert args.column in metadata.columns
print("Export cluster")
df_cluster=metadata[["barcode",args.column]]
df_cluster.to_csv(os.path.join(args.outdir,"cluster.csv"),sep=",",index=False)

print("Export embeding")
############### umap
umap_mat=adata.obsm["X_umap"]
umap_mat=pd.DataFrame(umap_mat,columns=["UMAP_1","UMAP_2"],index=adata.obs_names)
umap_mat["Barcode"]=umap_mat.index

df_umap=umap_mat.iloc[:,[2,0,1]]
df_umap.to_csv(os.path.join(args.outdir,"umap_projection.csv"),sep=",",index=False)
################ tsne
tsne_mat=adata.obsm["X_tsne"]
tsne_mat=pd.DataFrame(tsne_mat,columns=["tSNE_1","tSNE_2"],index=adata.obs_names)
tsne_mat["Barcode"]=tsne_mat.index

df_tsne=tsne_mat.iloc[:,[2,0,1]]
df_tsne.to_csv(os.path.join(args.outdir,"tsne_projection.csv"),sep=",",index=False)

#df_ident_cluster=metadata[["idents",args.column]]
#df_ident_cluster=metadata[["Sample",args.column]]
#df_ident_cluster.to_csv(os.path.join(args.outdir,"sample_cluster.csv"),sep=",",index=False)

#df_ident_cluster=metadata[["Status",args.column]]
#df_ident_cluster.to_csv(os.path.join(args.outdir,"status_cluster.csv"),sep=",",index=False)
