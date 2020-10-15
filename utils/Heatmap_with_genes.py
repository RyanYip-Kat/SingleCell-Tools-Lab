import scanpy as sc
import argparse
import sys
import os
import scanyuan as scy
import numpy as np
import pandas as pd
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
parser.add_argument("--groupby",type=str,default="status")
parser.add_argument("--genes",type=str,default=None)

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
sc.settings.set_figure_params(dpi_save=200)
sc.settings.figdir=args.outdir

viridisBig = cm.get_cmap('viridis', 512)
newcmp = ListedColormap(viridisBig(np.linspace(0.25, 0.75, 256)))

adata=sc.read_h5ad(args.data)
if args.order is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.order)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.column)
    #adata=adata.raw.to_adata()
    #adata.raw=adata
    sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.highly_variable_genes(adata,n_top_genes=5000)
    #adata= adata[:, adata.var.highly_variable]
    #sc.pp.scale(adata, max_value=10)
    adata.obs[args.column]=adata.obs[args.column].astype("category")
    adata.obs[args.column]=adata.obs[args.column].cat.reorder_categories(args.subset)

df=pd.read_csv(args.genes,sep="\t",header=None)
genes=df.loc[:,0].values.tolist()

#genes=np.unique(genes)


sc.pl.heatmap(adata,var_names=genes,groupby=args.groupby,
        swap_axes=True,show=False,show_gene_labels=True,cmap="viridis",
        vmin=-2.5, vmax=2.5,save="_genes_groups_viridis",figsize=[24,24])

sc.pl.stacked_violin(adata,var_names=genes,groupby=args.groupby,
        swap_axes=True,show=False,dendrogram=False,save="_genes_groups",figsize=[24,24])

scy.stacked_violin_t(adata, genes, figsize=[12,6], groupby=args.groupby,show=False,save="_stacked_violin_t",stripplot=False,jitter=False)
