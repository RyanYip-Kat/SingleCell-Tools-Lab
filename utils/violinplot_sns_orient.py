import scanpy as sc
import argparse
import sys
import os
import seaborn as sns
sys.path.append("../")
from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter
from core.plotting import violin_hue,violin_sns

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--genes",nargs="+",type=str,default=None)
parser.add_argument("--orient",type=str,default="h",help="")
parser.add_argument("--invert", action='store_true', default=False)
parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default=None,help="data")
parser.add_argument("--column2",type=str,default=None,help="data")
parser.add_argument("--groupby",type=str,default="status",help="")
parser.add_argument("--order",nargs="+",type=str,default=None,help="")
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir

genes=args.genes
adata=sc.read_h5ad(args.data)

if args.subset1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)

#if args.groupby=="status":
    #adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(["YA","AA"])
if args.order is not None:
    adata.obs[args.groupby]=adata.obs[args.groupby].cat.reorder_categories(args.order)
for gene in genes:
    #if args.groupby not in ["idents","status"]:
        sc.pl.violin(adata,gene,groupby=args.groupby,show=False,size=0.2,jitter=False,save="_"+gene+"_"+args.groupby+".pdf")
    #else:
        violin_sns(adata,gene_name=gene,save=True,groupby=args.groupby,figdir=args.outdir,orient=args.orient,strip=True)
