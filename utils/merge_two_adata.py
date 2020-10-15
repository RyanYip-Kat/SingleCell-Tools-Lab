import argparse
import os
import scanpy as sc
import numpy as np
import sys

sys.path.append("../")

from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--d1",type=str,default=None)
parser.add_argument("--c1",type=str,default=None)
parser.add_argument("--s1",nargs="+",type=str,default=None)

parser.add_argument("--d2",type=str,default=None)
parser.add_argument("--c2",type=str,default=None)
parser.add_argument("--s2",nargs="+",type=str,default=None)

parser.add_argument("--subset",nargs="+",type=str,default=None)
parser.add_argument("--level",nargs="+",type=str,default=None)
parser.add_argument("--genes",type=str,default=None,help="show genes")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.set_figure_params(dpi_save=200)
sc.settings.figdir=args.outdir

print("*** Loading data")
d1=sc.read_h5ad(args.d1)
d1=d1.raw.to_adata()

d2=sc.read_h5ad(args.d2)
d2=d2.raw.to_adata()

print("*** subset")
if args.c1 is not None and args.s1 is not None:
    d1=subset_by_column(d1,args.s1,args.c1,invert=False)

if args.c2 is not None and args.s2 is not None:
    d2=subset_by_column(d2,args.s2,args.c2,invert=False)

print("*** merge")
d1.obsm=None
d2.obsm=None

adata=d1.concatenate(d2)
sc.pp.normalize_total(adata, target_sum=1e4)
adata.obs["status"]=adata.obs["status"].astype("category")
#################

df=pd.read_csv(args.genes,sep="\t",header=None)
marker_genes=df.loc[:,0].values.tolist()

if args.subset is not None:
    adata=subset_by_column(adata,args.subset,"status")

if args.level is not None:
    adata.obs["status"]=adata.obs["status"].cat.reorder_categories(args.level)

sc.pl.dotplot(adata,var_names=marker_genes,groupby="status",color_map="Oranges",dendrogram=False,
        figsize=[16,8],show=False,save=True,smallest_dot=40,standard_scale='var')
