import scanpy as sc
import argparse
import sys
import os
import scanyuan as scy
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--order",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column",type=str,default="leiden",help="data")
parser.add_argument("--genes",nargs="+",type=str,default=None,help="genes")
parser.add_argument("--groupby",type=str,default="label_fine")
parser.add_argument("--theme",type=str,default=None)
#parser.add_argument("--size",type=int,default=15)

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
sc.settings.set_figure_params(frameon=False,dpi_save=600,figsize=[16,12])
sc.settings.figdir=args.outdir


adata=sc.read_h5ad(args.data)
if args.order is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.order)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)

sc.pl.umap(adata,color=args.groupby,palette=args.theme,size=28,show=False,save="_"+args.groupby+".pdf")
sc.pl.tsne(adata,color=args.groupby,palette=args.theme,size=28,show=False,save="_"+args.groupby+".pdf")

if args.genes is not None:
    for gene in args.genes:
        sc.pl.umap(adata,color=gene,size=24,show=False,save="_"+gene+".pdf")
        sc.pl.tsne(adata,color=gene,size=24,show=False,save="_"+gene+".pdf")

