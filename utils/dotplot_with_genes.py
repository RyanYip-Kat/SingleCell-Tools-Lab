import scanpy as sc
import argparse
import sys
import os
import scanyuan as scy
import pandas as pd

sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--genes",type=str,default=None,help="show genes")

parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default=None,help="data")
parser.add_argument("--column2",type=str,default=None,help="data")

parser.add_argument("--groupby",type=str,default="status")
parser.add_argument("--color",type=str,default="Oranges")
parser.add_argument("--level",nargs="+",type=str,default=None)


args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
sc.settings.figdir=args.outdir

df=pd.read_csv(args.genes,sep="\t",header=None)
marker_genes=df.loc[:,0].values.tolist()
adata=sc.read_h5ad(args.data)
if args.subset1 is not None and args.column1 is not None: 
    adata=subset_by_column(adata,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)

if args.level is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.level)

# "Oranges"
sc.pl.dotplot(adata,var_names=marker_genes,groupby=args.groupby,color_map=args.color,dendrogram=False,
        figsize=[12,4],show=False,save=True,dot_min=0,dot_max=1,smallest_dot=40,standard_scale='var') # smallest_dot=40,[16,8]

scy.stacked_violin_t(adata, marker_genes, figsize=[12,6], groupby=args.groupby,show=False,save="_stacked_violin_t",stripplot=False,jitter=False)
sc.pl.stacked_violin(adata,var_names=marker_genes,groupby=args.groupby,figsize=[12,6],show=False,save=True,stripplot=True,jitter=False)
