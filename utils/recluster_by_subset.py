import argparse
import os
import scanpy as sc
import scanpy.external as sce
import numpy as np
import sys

sys.path.append("../")

from core.functions import reduction,subset_by_column,subset_by_cell_feature,batch_correct
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--subset",nargs="+",type=str,default=None)
parser.add_argument("--column",type=str,default="leiden")
parser.add_argument("--batch_key",type=str,default=None)
parser.add_argument("--invert", action='store_true', default=False)
parser.add_argument("--batch_correct", action='store_true', default=False,help="batch correct")
parser.add_argument("--batch_method",type=str,default="combat",choices=["combat","mnn","bbknn","harmony"])
parser.add_argument("--resolution",type=float,default=1.2)
parser.add_argument("--use_rep",type=str,default="X_pca")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

########################
sc.settings.figdir=args.outdir
sc.settings.n_jobs=8
print("*** Loading data")
adata=sc.read_h5ad(args.data)

if args.column is not None and args.subset is not None:
    print("*** subset")
    print("Subset data by {}".format(args.column))
    adata=subset_by_column(adata,args.subset,args.column,invert=args.invert)
    #adata=adata.raw.to_adata()

if args.batch_correct and args.batch_key is not None:
    adata=batch_correct(adata,args.batch_key,args.batch_method)

print("*** tranform and reduction")
#adata=process(adata,args.batch_key)
adata=reduction(adata,resolution=args.resolution,use_rep=args.use_rep,batch_correct=args.batch_correct)
print(Counter(adata.obs["leiden"]))

print("*** find markers")
sc.tl.rank_genes_groups(adata,groupby="leiden", method='t-test',n_genes=500,corr_method="bonferroni")
adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')
