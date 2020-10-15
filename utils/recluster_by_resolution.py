import argparse
import os
import scanpy as sc
import numpy as np
import sys

sys.path.append("../")

from core.tl import reduction,process
from core.core import Model
from core.utils import subset_by_column,subset_by_cell_feature

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None)
parser.add_argument("--use_rep",type=str,default="X_pca")
parser.add_argument("-r","--resolution",type=float,default=1.5)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
sc.settings.n_jobs=8
print("*** Loading data")
adata=sc.read_h5ad(args.data)

print("**** clustering ****")
sc.pp.neighbors(adata, n_neighbors=20,use_rep=args.use_rep)
sc.tl.leiden(adata,resolution=args.resolution)
#sc.tl.louvain(adata,resolution=args.resolution)

print("*** find markers")
sc.tl.rank_genes_groups(adata,groupby="leiden", method='t-test',n_genes=500,corr_method="bonferroni")
adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')
