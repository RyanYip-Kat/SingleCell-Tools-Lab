import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default="leiden",help="data")
parser.add_argument("--column2",type=str,default=None,help="data")
parser.add_argument("--groupby",type=str,default="status",help="data")
parser.add_argument("--n_genes",type=int,default=500,help="the n genes to show")
parser.add_argument("--method",type=str,default="t-test",help="t-test_overestim_var,wilcoxon")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
if args.subset1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)

    #adata=adata.raw.to_adata()
    #adata.raw=adata
    #sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.highly_variable_genes(adata,n_top_genes=5000)
    #adata= adata[:, adata.var.highly_variable]
    #sc.pp.scale(adata, max_value=10)
sc.tl.rank_genes_groups(adata,groupby=args.groupby,method=args.method,n_genes=args.n_genes,corr_method="bonferroni")
get_rank_group_genes(adata,pval=None,fc=None,outdir=args.outdir,n_top=50)
get_rank_group_genes(adata,pval=None,fc=None,outdir=args.outdir)
