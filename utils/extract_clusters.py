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
parser.add_argument("-n","--name",type=str,default=None)
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
if args.name is not None:
    df_cluster[[args.column]]=args.name

df_cluster.to_csv(os.path.join(args.outdir,"cluster.csv"),sep=",",index=False)


