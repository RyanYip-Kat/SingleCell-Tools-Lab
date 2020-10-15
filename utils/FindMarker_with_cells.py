import scanpy as sc
import argparse
import sys
import os
import pandas as pd
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column,subset_by_cell_feature

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--cells",type=str,default=None,help="cells.csv")
parser.add_argument("--groupby",type=str,default="status",help="data")
parser.add_argument("--pval",type=float,default=None,help="pvalue")
parser.add_argument("--fc",type=float,default=None,help="logfc")
parser.add_argument("--method",type=str,default="t-test",help="t-test_overestim_var,wilcoxon")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
if args.cells is not None:
    df=pd.read_csv(args.cells,sep="\t",header=None)
    cells=df.loc[:,0].values.tolist()
    adata=subset_by_cell_feature(adata,cells=cells)

sc.tl.rank_genes_groups(adata,groupby=args.groupby,method=args.method,n_genes=500,corr_method="bonferroni")
get_rank_group_genes(adata,pval=args.pval,fc=args.fc,outdir=args.outdir,n_top=50)
get_rank_group_genes(adata,pval=args.pval,fc=args.fc,outdir=args.outdir)
