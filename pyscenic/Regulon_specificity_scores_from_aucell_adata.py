# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import argparse

import sys
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="aucell_adata.h5ad")
parser.add_argument("--cluster",type=str,default="CellType",choices=["Status","CellType","ClusterID","Leiden_clusters_Scanpy"],help="cellinfo")
parser.add_argument("--subset",nargs="+",type=str,default=None)
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("-"*16)
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=150)

print("-"*16)
# scenic output

adata=sc.read_h5ad(args.data)
if args.subset is not None:
    adata=subset_by_column(adata,args.subset,args.cluster)

auc_mtx=adata.to_df()
# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
    [
        pd.DataFrame( adata.obs.CellType, index=adata.obs.index),
        pd.DataFrame( adata.obs.ClusterID, index=adata.obs.index ),
        pd.DataFrame( adata.obs.Leiden_clusters_Scanpy, index=adata.obs.index ),
        pd.DataFrame( adata.obs.nGene, index=adata.obs.index),
        pd.DataFrame( adata.obs.nUMI, index=adata.obs.index ),
        pd.DataFrame( adata.obs.Status, index=adata.obs.index ),
    ],
    axis=1
)
cellAnnot.columns = [
 'CellType',
 'ClusterID',
 'Leiden_clusters_Scanpy',
 'nGene',
 'nUMI',
 "Status"]
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

rss_cellType = regulon_specificity_scores( auc_mtx, cellAnnot[args.cluster] )

print("#RSS panel plot with all cell types")
cats = sorted(list(set(cellAnnot[args.cluster])))
fig = plt.figure(figsize=(16, 10))
for c,num in zip(cats, range(1,len(cats)+1)):
    x=rss_cellType.T[c]
    ax = fig.add_subplot(2,2,num)
    plot_rss(rss_cellType, c, top_n=10, max_n=None, ax=ax)
    ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
        'figure.titlesize': 'large' ,
        'axes.labelsize': 'medium',
        'axes.titlesize':'large',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        })
plt.savefig(os.path.join(args.outdir,"cellType-RSS-top.pdf"), dpi=600, bbox_inches = "tight")

