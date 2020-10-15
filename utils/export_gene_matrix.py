import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys

sys.path.append("../")

from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--n_top",type=int,default=10000)
parser.add_argument("--subset",nargs="+",type=str,default=None)
parser.add_argument("--column",type=str,default="leiden")
parser.add_argument("--invert", action='store_true', default=False)


args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
print("*** Loading data")
adata=sc.read_h5ad(args.data)

print("Subset data by {}".format(args.column))
if args.column is not None and args.subset is not None:
    print("*** subset")
    print("Subset data by {}".format(args.column))
    adata=subset_by_column(adata,args.subset,args.column,invert=args.invert)

if args.n_top is not None:
    adata=adata.raw.to_adata()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=args.n_top)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]


X=adata.X.toarray()
print(X.shape)
matrix=pd.DataFrame(X.transpose(),columns=adata.obs_names,index=adata.var_names)

filename=os.path.join(args.outdir,"expression_matrix.txt")
print("*** Expression matrix output in {}".format(filename))
matrix.to_csv(filename,sep="\t")

print("*** Done! ***")
