import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import argparse
import sys
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column
parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column",type=str,default="leiden",help="data")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
sc.settings.set_figure_params(dpi_save=200)
sc.settings.figdir=args.outdir

print("### Loading dataset")
adata=sc.read_h5ad(args.data)
metadata=adata.obs

if args.subset is not None and args.column in metadata.columns:
    adata=subset_by_column(adata,args.subset,args.column)

if hasattr(adata,"raw"):
    nGene=np.array( np.sum(adata.raw.X.transpose()>0 , axis=0)).flatten()
    nUMI=np.array( np.sum(adata.raw.X.transpose() , axis=0)).flatten()
    X=adata.raw.X.transpose()
    row_attrs = {
    "Gene": np.array(adata.raw.var_names)
    }
    col_attrs = {
            "CellID":  np.array(adata.obs_names) ,
            "nGene": nGene ,
            "nUMI": nUMI 
            }

else:
    nGene=np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten()
    nUMI=np.array( np.sum(adata.X.transpose() , axis=0)).flatten()
    X=adata.X.transpose()
#################
    row_attrs = { 
    "Gene": np.array(adata.var_names)
    }
    col_attrs = { 
            "CellID":  np.array(adata.obs_names) ,
            "nGene": nGene ,
            "nUMI": nUMI}

loom_path=os.path.join(args.outdir,"scenic.loom")
print("### Save loom data into {}".format(loom_path))
lp.create(loom_path,X, row_attrs, col_attrs )
