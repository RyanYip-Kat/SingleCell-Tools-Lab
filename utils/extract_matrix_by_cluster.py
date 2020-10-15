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
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None,help="adata object")
parser.add_argument("--column",type=str,default="leiden")
parser.add_argument("--use_raw",action="store_true",default=False)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)


print("### Laoding data")
adata=sc.read_h5ad(args.data)


metadata=adata.obs.copy()
metadata["barcode"]=metadata.index

assert args.column in metadata.columns
cluster=np.unique(metadata[[args.column]].values.tolist())

if args.use_raw:
    X=adata.raw.X.toarray()
else:
    X=adata.X

for c in cluster:
    idx=metadata[args.column].isin([c])
    m=X[idx,:]
    print("The size of matrix from : {} is : nrow {},ncol {}".format(c,m.shape[0],m.shape[1]))
    filename=os.path.join(args.outdir,"cluster_"+str(c)+'_data_array.npy')
    print("Save {} matrix array into {}".format(c,filename))
    np.save(filename,m)

print("Done!")

