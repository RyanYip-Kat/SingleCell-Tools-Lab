import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--groupby",type=str,default="status",help="data")
parser.add_argument("--column",type=str,default="label_main",help="data")
parser.add_argument("--method",type=str,default="wilcoxon",help="t-test_overestim_var,wilcoxon")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)


labels=[cell for cell in adata.obs[args.column].unique()]
print(labels)
for celltype in labels:
    c=[]
    c.append(celltype)
    print("Get from : {}".format(celltype))
    data=subset_by_column(adata,c,args.column)
    sc.pp.normalize_total(data, target_sum=1e4)

    sc.tl.rank_genes_groups(data,groupby=args.groupby,method=args.method,n_genes=500,corr_method="bonferroni")
    sc.tl.filter_rank_genes_groups(data,groupby=args.groupby,key_added='rank_genes_groups_filtered',min_fold_change=0.25)
    out=os.path.join(args.outdir,celltype)

    if not os.path.exists(out):
        os.makedirs(out)
    get_rank_group_genes(data,pval=None,fc=0.25,outdir=out,n_top=50,key="rank_genes_groups_filtered")
    get_rank_group_genes(data,pval=None,fc=0.25,outdir=out,key="rank_genes_groups_filtered")

