import argparse
import os
import scanpy as sc
import scanpy.external as sce
import scvelo as scv
import numpy as np
import sys

sys.path.append("../")

from core.functions import reduction,subset_by_column,subset_by_cell_feature,batch_correct
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)

parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default=None,help="data")
parser.add_argument("--column2",type=str,default=None,help="data")

parser.add_argument("--dynamics",action="store_true",default=False)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

########################
sc.settings.figdir=args.outdir
sc.settings.n_jobs=8
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

print("### Loading data")
adata=sc.read_h5ad(args.data)

if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)


print("### Run monments")
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
mode="stochastic"
if args.dynamics:
    print(" Run recover_dynamics")
    scv.tl.recover_dynamics(adata)
    mode="dynamical"
scv.tl.velocity(adata, mode=mode)

print("###  Run velocity")
scv.tl.velocity_graph(adata)
#scv.pl.velocity_embedding_stream(adata, basis='umap',save=)

print("### Run latent time")
scv.tl.latent_time(adata)

print("### Run pseudotime")
scv.tl.velocity_pseudotime(adata)

print("### Save adata")
adata.write(os.path.join(args.outdir,"adata_velocity.h5ad"),compression='gzip')
