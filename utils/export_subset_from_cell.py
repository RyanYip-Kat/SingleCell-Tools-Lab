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
parser.add_argument("--barcode",type=str,default=None,help="barcode csv")
parser.add_argument("--invert", action='store_true', default=False)
parser.add_argument("--use_raw",action="store_true",default=False)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
print("*** Loading data")
adata=sc.read_h5ad(args.data)
Data=pd.read_csv(args.barcode)
Data.columns=["barcode","label"]
cells=Data["barcode"].values.tolist()
print("*** subset")
adata=subset_by_cell_feature(adata,cells=cells,invert=args.invert)
if args.use_raw:
    adata=adata.raw.to_adata()

data=adata.copy()
print("*** save")
data.write(os.path.join(args.outdir,"adata.h5ad"),compression='gzip')

