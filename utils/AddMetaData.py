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
#parser.add_argument("--subset",nargs="+",type=str,default=None)
parser.add_argument("--column",type=str,default="new_cluster")
parser.add_argument("--metadata",type=str,default=None)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("*** Loading data")
adata=sc.read_h5ad(args.data)
DATA=pd.read_csv(args.metadata)

assert DATA.shape[1]==2
DATA.columns=["barcode","label_fine"]
DATA.index=DATA.barcode.values
DATA.drop("barcode",axis=1,inplace=True)

metadata=DATA.loc[adata.obs_names,:]
meta=metadata.iloc[:,0].values.tolist()

if args.column is None:
    column=metadata.columns.tolist()[0]
else:
    column=args.column

print("Add column : {} into adata metadata".format(column))
adata.obs[column]=meta

print("Save")
adata.write(os.path.join(args.outdir,"adata.h5ad"))
